 /* Include constant curvature drivers for all dimensions here */
 
/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives2d.h"
#include "qss_tvd_runge_kutta2d.h"
#include "qss_data_arrays.h"
#include "qss_util2d.h"
#include "qss_macros.h"
#include "qss_reinitialization2d.h"
#include "qss_grid.h"
#include "constCurvModel2d.h"
#include "qss_general_util.h"
#include "connectivity.h"
#include "reinitialize_top.h"

/* Main driver for constant curvature level set method model */

QSSLIB_REAL constCurvModel2d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
{

    QSSLIB_REAL     zero = 0.0;
    char fname[256];
  
    int flag = 0, OUTER_STEP = 0, INNER_STEP, idx;     
  
    QSSLIB_REAL dt = 0, dt_sub, t = 0;
    QSSLIB_REAL   mask_sign = -1;
                                                                                                                  
    QSSLIB_REAL eps, cur_max_H_over_dX = -1, cfl_number = 0.5;
    QSSLIB_REAL max_abs_err = 100, vol_phi = 100, vol_very_small;
    QSSLIB_REAL vol_pore_sp, satn, satn_old, max_satn;
    
    max_satn = options->vol_frac_max;
    
    vol_very_small = (g->dx)[0]*(g->dx)[1];    
    eps = (options->eps_coefficient)*(g->dx[0]);
    
    QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_pore_sp, p->mask, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);
		        
    QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);
		        
    satn = vol_phi/vol_pore_sp;
    
    while( (t < options->tmax) && (max_abs_err > options->eps_stop) && (vol_phi > vol_very_small) && (satn < max_satn) )
    { 
        /* outer loop */
        OUTER_STEP++;
        dt_sub = 0;
        COPY_DATA(p->phi_prev,p->phi,g)
        
        /* Begin parallel region */
        #pragma omp parallel default(none) shared(p, g, flag, cfl_number,\
            cur_max_H_over_dX, zero, dt_sub, vol_phi, max_abs_err, eps, options, fp_out) \
            private(INNER_STEP, dt)
        {

            INNER_STEP = 0;
            QSSLIB_REAL max_H_over_dX;
            int    bdry_location_idx = 9; /* all boundaries */
            QSSLIB_REAL disconn_overlap = 0;
    
            /* Set up variables for multi-threading */
            int cur_thread, cur_jlo_fb, cur_jhi_fb, num_threads, nslices, i;
            int cur_jlo_gb, cur_jhi_gb;
                       
            cur_thread = omp_get_thread_num();
            num_threads = omp_get_num_threads();
    
            if (num_threads > 1) flag = 1;

            nslices = g->jhi_fb - g->jlo_fb + 1;

            cur_jlo_fb = g->jlo_fb + nslices*cur_thread/num_threads;
            cur_jhi_fb = g->jlo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            double t1 = omp_get_wtime();
            
            /* Keeping track of thread-local ghost boundaries, mainly for imposing mask */
            if (cur_jhi_fb > (g->jhi_fb))
	            cur_jhi_fb = (g->jhi_fb);

            if (cur_thread == 0)
                cur_jlo_gb = cur_jlo_fb - 3;
            else
                cur_jlo_gb = cur_jlo_fb;
            
            if (cur_thread == (num_threads - 1))
                cur_jhi_gb = cur_jhi_fb + 3;
            else
                cur_jhi_gb = cur_jhi_fb;
            

            while( dt_sub < options->tplot )
            { 
                /* inner loop */
                INNER_STEP++;
                      
                (options->space_deriv_func)(p->phi, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));	 

	            #pragma omp barrier
	            
	            #pragma omp critical
	            {
	                if(max_H_over_dX > cur_max_H_over_dX)
	                    cur_max_H_over_dX = max_H_over_dX;
	            }
	
	            /* 
	               Barrier to ensure same cur_max_H_over_dX across all threads.
	            */
	            #pragma omp barrier	            
	            
	            /* get final correct dt due to parabolic (curvature) term */
	            dt =  cfl_number / (cur_max_H_over_dX + options->b_max_over_dx + options->max_U_over_dx);
		                
	            QSS2D_RK1_STEP(p->phi_next,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);

                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
                /* phi_next should now store the RK_stage1 output */
	            if (options->order_time_accur >= 2) {
	                (options->space_deriv_func)(p->phi_next, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));
	                
	                #pragma omp barrier
	                
	                QSS2D_TVD_RK2_STAGE2(p->scratch1,GB_DIMS_2D,p->phi_next,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);
		            
		            #pragma omp barrier
		            /* boundary conditions */
                     #pragma omp single
                        signedLinearExtrapolationBCqss(p->scratch1,g,bdry_location_idx);         
                     #pragma omp single
                        COPY_DATA(p->phi_next, p->scratch1, g);	
                    /* scratch1 should now have RK_stage2 output. If third order is not desired, phi_next 
                       stores final output */

	            }
	            
	            /* 3rd order is the highest time order accuracy handled */
                if (options->order_time_accur >= 3) {
	                (options->space_deriv_func)(p->phi_next, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));
	                
	                #pragma omp barrier
	                
	                QSS2D_TVD_RK3_STAGE3(p->scratch2,GB_DIMS_2D,p->scratch1,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);
		        
		            #pragma omp barrier
		            /* boundary conditions */
                     #pragma omp single
                        signedLinearExtrapolationBCqss(p->scratch2,g,bdry_location_idx);
                     #pragma omp single
                        COPY_DATA(p->phi_next, p->scratch2, g);	
                    /* scratch3 should now have RK_stage3 output. phi_next stores final output */
                 }
                 
		        #pragma omp barrier
	            IMPOSE_MASK_PAR_2D(p->phi, p->mask, p->phi_next, &(options->overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	            
	            #pragma omp barrier
	            if(options->check_connectivity) IMPOSE_MASK_PAR_2D(p->phi, p->mask_w, p->phi, &(disconn_overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	            #pragma omp barrier
	              
	            if(options->check_connectivity) IMPOSE_MASK_PAR_2D(p->phi, p->mask_nw, p->phi, &(disconn_overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	              
            	#pragma omp barrier

                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
	            /* 
	               Implicit barrier after omp single, so all threads should sync here
	            */
	             #pragma omp single
	            {
	                SET_DATA_TO_CONSTANT(p->lse_rhs,g,zero);
	                dt_sub += dt;
	            }
	            	    
            }   /* End inner loop */
            
            double t2 = omp_get_wtime();
            
            if (cur_thread == 0)
	            fprintf(fp_out, "%d threads: Level set time = %lf\n", num_threads, t2 - t1);
	            
        }   /* End parallel region */
            
            /* 
                Reinitialization of the level set function 
                - may want to parallelize the following functions later
            */
        
        t += dt_sub;
        
        double t3 = omp_get_wtime();
        
        
        IMPOSE_MASK(p->phi, p->mask, p->phi, g);
        
        if(options->check_connectivity) { //if (((int)t % 1) == 0) 
            trapComponents_mask(p, g, options);
            reinitializeSubcellFix2d(p->mask_w,g,options);
            reinitializeSubcellFix2d(p->mask_nw,g,options);
            
            /* update velocities based on trapped components */
            if (options->conserve_imbibe)
                IMPOSE_TRAP_VEL_2D(p->mask_nw, p->mask_w, p->normal_velocity, p->curvature_coeff, p->external_velocity_x,
                    p->external_velocity_y, GB_DIMS_2D, FB_DIMS_2D, &((g->dx)[0]));   
        }
        
        /* merge trapped nw phase for ensuring reinitialization, and correct saturations */
	    MERGE_SETS(p->phi, p->mask_nw, g);
	    
        fprintf(fp_out, "Reinitializing....");
        reinitializeSubcellFix2d(p->phi,g,options);
        fprintf(fp_out, "Reinitialized\n");

        /* compute stopping criteria */    
        /* max abs error */
        QSS2D_MAX_NORM_DIFF_LOCAL(&max_abs_err,p->phi,GB_DIMS_2D, p->phi_prev,
		        GB_DIMS_2D, FB_DIMS_2D, &(options->err_check_zone));
	    	    
	    satn_old = satn;
	    
	    
        QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);
        
        satn = vol_phi/vol_pore_sp;
                
        if ( (satn > max_satn) || (vol_phi <= vol_very_small) )
		    return -1;
		    
        if(options->checkpoint)
        {
            sprintf(fname,"checkpoint_phi");
            writeDataArrayQSS(p->phi,g,fname,GZIP);
	        sprintf(fname,"checkpoint_phi_prev");
            writeDataArrayQSS(p->phi_prev,g,fname,GZIP);
            sprintf(fname,"mask_w");
            writeDataArrayQSS(p->mask_w,g,fname,GZIP);
            sprintf(fname,"mask_nw");
            writeDataArrayQSS(p->mask_nw,g,fname,GZIP);
            sprintf(fname,"curvature_coeff");
            writeDataArrayQSS(p->curvature_coeff,g,fname,GZIP);
        }

	    double t4 = omp_get_wtime();
        fprintf(fp_out, "connectivity time = %lf\n", t4 - t3);
        
        fprintf(fp_out, "t = %f\t",t);
        fprintf(fp_out, "max_abs_err = %4.3f,\t",max_abs_err);
        fprintf(fp_out, "nw phase satn pct = %4.3f\n",100*satn);

        printf("t = %4.1f\t",t);
        printf("max_abs_err = %4.3f,\t",max_abs_err);
        printf("nw phase satn pct = %4.3f\n", 100*satn);
        
        /* If nw phase saturation isn't changing much, then continue on */
        if (options->use_satn_stop) {
            if (fabsf(100*(satn - satn_old)) < 0.001)
                break;
            }
    } /* End outer loop */
    
    /* Merge disconnected components for writing to saved data */
    //MERGE_SETS(p->phi, p->mask_nw, g);
   
    return t;
}

QSSLIB_REAL constCurvModel2dNoVar(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out, QSSLIB_REAL a) 
{

    QSSLIB_REAL     zero = 0.0;
    char fname[256];
  
    int flag = 0, OUTER_STEP = 0, INNER_STEP, idx;     
  
    QSSLIB_REAL dt = 0, dt_sub, t = 0;
    QSSLIB_REAL   mask_sign = -1;
                                                                                                                  
    QSSLIB_REAL eps, cur_max_H_over_dX = -1, cfl_number = 0.5;
    QSSLIB_REAL max_abs_err = 100, vol_phi = 100, vol_very_small;
    QSSLIB_REAL vol_pore_sp, satn, satn_old, max_satn;
    
    vol_very_small = (g->dx)[0]*(g->dx)[1];    
    eps = (options->eps_coefficient)*(g->dx[0]);
    
    max_satn = options->vol_frac_max;
    
    QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_pore_sp, p->mask, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);
		        
    QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);
		        
    satn = vol_phi/vol_pore_sp;
    
    while( (t < options->tmax) && (max_abs_err > options->eps_stop) && (vol_phi > vol_very_small) && (satn < max_satn) )
    { 
        /* outer loop */
        OUTER_STEP++;
        dt_sub = 0;
        COPY_DATA(p->phi_prev,p->phi,g)
        
        /* Begin parallel region */
        #pragma omp parallel default(none) shared(p, g, a, flag, cfl_number,\
            cur_max_H_over_dX, zero, dt_sub, vol_phi, max_abs_err, eps, options, fp_out) \
            private(INNER_STEP, dt)
        {

            INNER_STEP = 0;
            QSSLIB_REAL max_H_over_dX;
            int    bdry_location_idx = 9; /* all boundaries */
            QSSLIB_REAL disconn_overlap = 0;
    
            /* Set up variables for multi-threading */
            int cur_thread, cur_jlo_fb, cur_jhi_fb, num_threads, nslices, i;
            int cur_jlo_gb, cur_jhi_gb;
                       
            cur_thread = omp_get_thread_num();
            num_threads = omp_get_num_threads();
    
            if (num_threads > 1) flag = 1;

            nslices = g->jhi_fb - g->jlo_fb + 1;

            cur_jlo_fb = g->jlo_fb + nslices*cur_thread/num_threads;
            cur_jhi_fb = g->jlo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            double t1 = omp_get_wtime();
            
            /* Keeping track of thread-local ghost boundaries, mainly for imposing mask */
            if (cur_jhi_fb > (g->jhi_fb))
	            cur_jhi_fb = (g->jhi_fb);

            if (cur_thread == 0)
                cur_jlo_gb = cur_jlo_fb - 3;
            else
                cur_jlo_gb = cur_jlo_fb;
            
            if (cur_thread == (num_threads - 1))
                cur_jhi_gb = cur_jhi_fb + 3;
            else
                cur_jhi_gb = cur_jhi_fb;
            

            while( dt_sub < options->tplot )
            { 
                /* inner loop */
                INNER_STEP++;
                      
                (options->space_deriv_func)(p->phi, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));	 

	            #pragma omp barrier
	            
	            #pragma omp critical
	            {
	                if(max_H_over_dX > cur_max_H_over_dX)
	                    cur_max_H_over_dX = max_H_over_dX;
	            }
	
	            /* 
	               Barrier to ensure same cur_max_H_over_dX across all threads.
	            */
	            #pragma omp barrier	            
	            
	            /* get final correct dt due to parabolic (curvature) term */
	            dt =  cfl_number / (cur_max_H_over_dX + options->b_max_over_dx + options->max_U_over_dx);
		                
	            QSS2D_RK1_STEP(p->phi_next,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);

                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
                /* phi_next should now store the RK_stage1 output */
	            if (options->order_time_accur >= 2) {
	                (options->space_deriv_func)(p->phi_next, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));
	                
	                #pragma omp barrier
	                
	                QSS2D_TVD_RK2_STAGE2(p->scratch1,GB_DIMS_2D,p->phi_next,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);
		            
		            #pragma omp barrier
		            /* boundary conditions */
                     #pragma omp single
                        signedLinearExtrapolationBCqss(p->scratch1,g,bdry_location_idx);         
                     #pragma omp single
                        COPY_DATA(p->phi_next, p->scratch1, g);	
                    /* scratch1 should now have RK_stage2 output. If third order is not desired, phi_next 
                       stores final output */

	            }
	            
	            /* 3rd order is the highest time order accuracy handled */
                if (options->order_time_accur >= 3) {
	                (options->space_deriv_func)(p->phi_next, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS_2D, FB_DIMS_PAR_2D,
	                &((g->dx)[0]),&((g->dx)[1]));
	                
	                #pragma omp barrier
	                
	                QSS2D_TVD_RK3_STAGE3(p->scratch2,GB_DIMS_2D,p->scratch1,GB_DIMS_2D,p->phi,GB_DIMS_2D,p->lse_rhs,
		            GB_DIMS_2D, FB_DIMS_PAR_2D, &dt);
		        
		            #pragma omp barrier
		            /* boundary conditions */
                     #pragma omp single
                        signedLinearExtrapolationBCqss(p->scratch2,g,bdry_location_idx);
                     #pragma omp single
                        COPY_DATA(p->phi_next, p->scratch2, g);	
                    /* scratch3 should now have RK_stage3 output. phi_next stores final output */
                 }
                 
		        #pragma omp barrier
	            IMPOSE_MASK_PAR_2D(p->phi, p->mask, p->phi_next, &(options->overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	            
	            #pragma omp barrier
	            if(options->check_connectivity) IMPOSE_MASK_PAR_2D(p->phi, p->mask_w, p->phi, &(disconn_overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	            #pragma omp barrier
	              
	            if(options->check_connectivity) IMPOSE_MASK_PAR_2D(p->phi, p->mask_nw, p->phi, &(disconn_overlap),
	              GB_DIMS_2D, &(cur_jlo_gb), &(cur_jhi_gb));
	              
            	#pragma omp barrier

                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
	            /* 
	               Implicit barrier after omp single, so all threads should sync here
	            */
	             #pragma omp single
	            {
	                SET_DATA_TO_CONSTANT(p->lse_rhs,g,zero);
	                dt_sub += dt;
	            }
	            	    
            }   /* End inner loop */
            
            double t2 = omp_get_wtime();
            
            if (cur_thread == 0)
	            fprintf(fp_out, "%d threads: Level set time = %lf\n", num_threads, t2 - t1);
	            
        }   /* End parallel region */
            
            /* 
                Reinitialization of the level set function 
                - may want to parallelize the following functions later
            */
        
        t += dt_sub;
        
        double t3 = omp_get_wtime();
        
        if(options->check_connectivity) {//if (((int)t % 1) == 0) 
            trapComponents_mask(p, g, options);
            reinitializeSubcellFix2d(p->mask_w,g,options);
            reinitializeSubcellFix2d(p->mask_nw,g,options);
        }
        
        /* merge trapped nw phase for ensuring reinitialization */
	    MERGE_SETS(p->phi, p->mask_nw, g);
	    
        fprintf(fp_out, "Reinitializing....");
        //reinitialize2d_subcell_fix_qss(p,g,options);
        reinitializeSubcellFix2d(p->phi,g,options);
        fprintf(fp_out, "Reinitialized\n");

        /* compute stopping criteria */    
        /* max abs error */
        QSS2D_MAX_NORM_DIFF_LOCAL(&max_abs_err,p->phi,GB_DIMS_2D, p->phi_prev,
		        GB_DIMS_2D, FB_DIMS_2D, &(options->err_check_zone));
	    	    
	    satn_old = satn;
	    
        QSS2D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS_2D, FB_DIMS_2D, 
		        &(g->dx[0]),&(g->dx[1]), &eps);

        satn = vol_phi/vol_pore_sp;
                
        if ( (satn > max_satn) || (vol_phi <= vol_very_small) )
		    return -1;
		    
        if(options->checkpoint)
        {
            sprintf(fname,"checkpoint_phi");
            writeDataArrayQSS(p->phi,g,fname,GZIP);
	        sprintf(fname,"checkpoint_phi_prev");
            writeDataArrayQSS(p->phi_prev,g,fname,GZIP);
            sprintf(fname,"mask_w");
            writeDataArrayQSS(p->mask_w,g,fname,GZIP);
            sprintf(fname,"mask_nw");
            writeDataArrayQSS(p->mask_nw,g,fname,GZIP);
        }

	    double t4 = omp_get_wtime();
        fprintf(fp_out, "connectivity time = %lf\n", t4 - t3);
        
        fprintf(fp_out, "t = %f\t",t);
        fprintf(fp_out, "max_abs_err = %4.3f,\t",max_abs_err);
        fprintf(fp_out, "nw phase satn pct = %4.3f\n", 100*satn);

        printf("t = %4.1f\t",t);
        printf("max_abs_err = %4.3f,\t",max_abs_err);
        printf("nw phase satn pct = %4.3f\n", 100*satn);
        
        /* If nw phase saturation isn't changing much, then continue on */
        if (options->use_satn_stop) {
            if (fabsf(100*(satn - satn_old)) < 0.001)
                break;
            }
    } /* End outer loop */
    
    /* Merge disconnected components for writing to saved data */
    //MERGE_SETS(p->phi, p->mask_nw, g);
   
    return t;
}


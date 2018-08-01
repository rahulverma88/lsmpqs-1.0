/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file constCurvModel3d.c

    Function definitions for performing level set computations to reach constant curvature solutions
    in three dimensions.
    
    Both drainage and imbibition simulations call these functions.
    
    Sets up time-stepping parameters based on options structure, and then sets up parallelization.
    Domain decomposition for parallelization is done based on a simple bread-slicing 
    plan in the z-direction. 
    
    There are two functions - constCurvModel3d, and constCurvModel3dNoVar.
    The latter is for zero contact angle cases, where "a" and "b" parameters in the level
    set equation do not need to vary, and the convective "V" term does not exist. So not creating
    matrices for these terms results in significant memory and computational time savings,
    especially for large geometries in 3 dimensions.
             
*/
/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_reinitialization3d.h"
#include "qss_grid.h"
#include "qss_macros.h"
#include "constCurvModel3d.h"
#include "reinitialize_top.h"
#include "qss_general_util.h"
#include "connectivity.h"

#define DSZ  sizeof(QSSLIB_REAL)

/* Main driver for constant curvature level set method model */

QSSLIB_REAL constCurvModel3d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
{

    QSSLIB_REAL     zero = 0.0;
    char fname[256];
  
    int flag = 0, OUTER_STEP = 0, INNER_STEP, idx;     
  
    QSSLIB_REAL dt = 0, dt_sub, t = 0;
    QSSLIB_REAL   mask_sign = -1;
                                                                                                    
    QSSLIB_REAL eps, cur_max_H_over_dX = -1, cfl_number = 0.5;
    QSSLIB_REAL max_abs_err = 100, vol_phi = 100, vol_very_small;
    QSSLIB_REAL vol_pore_sp, satn, satn_old, max_satn;
    
    vol_very_small = (g->dx)[0]*(g->dx)[1]*(g->dx)[2];    
    eps = (options->eps_coefficient)*(g->dx[0]);

    max_satn = options->vol_frac_max;
    
    QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_pore_sp, p->mask, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]),&(g->dx[2]), &eps);
		        
    QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]), &(g->dx[2]), &eps);
	
	satn = vol_phi/vol_pore_sp;
	
    while( (t < options->tmax) && (max_abs_err > options->eps_stop) )
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
            int cur_thread, cur_klo_fb, cur_khi_fb, num_threads, nslices, i;
            int cur_klo_gb, cur_khi_gb;
                       
            cur_thread = omp_get_thread_num();
            num_threads = omp_get_num_threads();
            
            if (num_threads > 1) flag = 1;

            nslices = g->khi_fb - g->klo_fb + 1;

            cur_klo_fb = g->klo_fb + nslices*cur_thread/num_threads;
            cur_khi_fb = g->klo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            double t1 = omp_get_wtime();
            
            /* Keeping track of thread-local ghost boundaries, mainly for imposing mask */
            if (cur_khi_fb > (g->khi_fb))
	            cur_khi_fb = (g->khi_fb);

            if (cur_thread == 0)
                cur_klo_gb = cur_klo_fb - 3;
            else
                cur_klo_gb = cur_klo_fb;
            
            if (cur_thread == (num_threads - 1))
                cur_khi_gb = cur_khi_fb + 3;
            else
                cur_khi_gb = cur_khi_fb;
            

            while( dt_sub < options->tplot )
            { 
                /* inner loop */
                INNER_STEP++;
                      
                (options->space_deriv_func)(p->phi, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y, p->external_velocity_z,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	 
	                
	            #pragma omp barrier
	            
	            #pragma omp critical
	            {
	                if(max_H_over_dX > cur_max_H_over_dX)
	                    cur_max_H_over_dX = max_H_over_dX;
	            }
	
	            /* Barrier to ensure same cur_max_H_over_dX across all threads. */
	            #pragma omp barrier	            
	            
	            /* get final correct dt due to parabolic (curvature) term */
	            dt =  cfl_number / (cur_max_H_over_dX + options->b_max_over_dx + options->max_U_over_dx);

	            QSS3D_RK1_STEP(p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);

                #pragma omp barrier
                
                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
	            
	            /* phi_next should now store the RK_stage1 output */
	            if (options->order_time_accur >= 2) {
	                (options->space_deriv_func)(p->phi_next, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y, p->external_velocity_z,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	                
	                #pragma omp barrier
	                
	                QSS3D_TVD_RK2_STAGE2(p->scratch2,GB_DIMS,p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);
		            
		            #pragma omp barrier
		            #pragma omp single
                    {
                       COPY_DATA(p->phi_next, p->scratch2, g);	
                    }
		            /* boundary conditions */
                     #pragma omp single
                    {
                        signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);       
                    }  
                     
                    /* scratch1 should now have RK_stage2 output. If third order is not desired, phi_next 
                       stores final output */
                }
                
                /* 3rd order is the highest time order accuracy handled */
                if (options->order_time_accur >= 3) {
	                (options->space_deriv_func)(p->phi_next, p->normal_velocity, p->curvature_coeff,
                    p->external_velocity_x, p->external_velocity_y, p->external_velocity_z,
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	                
	                #pragma omp barrier
	                
	                QSS3D_TVD_RK3_STAGE3(p->scratch2,GB_DIMS,p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);
		        
		            #pragma omp barrier
		            #pragma omp single
		            {
                        COPY_DATA(p->phi_next, p->scratch2, g);	
                    }
		            /* boundary conditions */
                     #pragma omp single
                    {
                        signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);
                    }
                     
                    /* scratch3 should now have RK_stage3 output. phi_next stores final output */
                 }
		        
		         #pragma omp single
	             {
	                SET_DATA_TO_CONSTANT(p->lse_rhs,g,zero);
	                SET_DATA_TO_CONSTANT(p->scratch1,g,zero);
	                SET_DATA_TO_CONSTANT(p->scratch2,g,zero);
	                dt_sub += dt;
	             }  
		            
		        /* Since all final outputs are stored in phi_next, IMPOSE_MASK works for all time stepping methods */
	            #pragma omp barrier
	            if(~options->use_var_theta)
	                IMPOSE_MASK_PAR(p->phi, p->mask, p->phi_next, &(options->overlap),
	                        GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            else
	                IMPOSE_MASK_PAR_VAR(p->phi, p->mask, p->phi_next, p->overlap,
	                        GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	            if(options->check_connectivity) IMPOSE_MASK_PAR(p->phi, p->mask_w, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	            if(options->check_connectivity) IMPOSE_MASK_PAR(p->phi, p->mask_nw, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	              
            	#pragma omp barrier
	            	    
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
        
        if(options->check_connectivity) {
            trapComponents_mask(p, g, options);
            reinitializeSubcellFix3d(p->mask_w, g, options);
            reinitializeSubcellFix3d(p->mask_nw,g,options);
            
            /* update velocities based on trapped components */
            if (options->conserve_imbibe)
                IMPOSE_TRAP_VEL_3D(p->mask_nw, p->mask_w, p->normal_velocity, p->curvature_coeff, p->external_velocity_x,
                    p->external_velocity_y, p->external_velocity_z, GB_DIMS, FB_DIMS, &((g->dx)[0]));
            
	        MERGE_SETS(p->phi, p->mask_nw, g);
        }
        
	    
        fprintf(fp_out, "Reinitializing....");
        reinitializeSubcellFix3d(p->phi,g,options);
        fprintf(fp_out, "Reinitialized\n");
        
        /* compute stopping criteria */
        /* max abs error */
        QSS3D_MAX_NORM_DIFF_LOCAL(&max_abs_err,p->phi,GB_DIMS, p->phi_prev,
		        GB_DIMS, FB_DIMS, &(options->err_check_zone));

        QSSLIB_REAL vol_phi_old = vol_phi;
        
        satn_old = satn;
        
        QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]),&(g->dx[2]), &eps);
		
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

        fprintf(fp_out, "t = %4.1f\t",t);
        fprintf(fp_out, "max_abs_err = %4.3f,\t",max_abs_err);
        fprintf(fp_out, "nw phase satn pct = %4.3f\n", satn*100);
        
        printf("t = %4.1f\t",t);
        printf("max_abs_err = %4.3f,\t",max_abs_err);
        printf("nw phase satn = %4.3f\n", satn*100);
        fflush(fp_out);

        if (options->use_satn_stop) {
            if (fabsf(100*(satn - satn_old)) < 0.001)
                break;
            }
    } /* End outer loop */
    
    //MERGE_SETS(p->phi, p->mask_nw, g);
    
    return t;
}

QSSLIB_REAL constCurvModel3dNoVar(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out, QSSLIB_REAL a) 
{

    QSSLIB_REAL     zero = 0.0;
    char fname[256];
  
    int flag = 0, OUTER_STEP = 0, INNER_STEP, idx;     
  
    QSSLIB_REAL dt = 0, dt_sub, t = 0;
    QSSLIB_REAL   mask_sign = -1;
                                                                                                    
    QSSLIB_REAL eps, cur_max_H_over_dX = -1, cfl_number = 0.5;
    QSSLIB_REAL max_abs_err = 100, vol_phi = 100, vol_very_small;
    QSSLIB_REAL vol_pore_sp, satn, satn_old, max_satn;
    
    vol_very_small = (g->dx)[0]*(g->dx)[1]*(g->dx)[2];    
    eps = (options->eps_coefficient)*(g->dx[0]);

    max_satn = options->vol_frac_max;
    
    QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_pore_sp, p->mask, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]),&(g->dx[2]), &eps);
		        
    QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]), &(g->dx[2]), &eps);
	
	satn = vol_phi/vol_pore_sp;
	
    while( (t < options->tmax)  && (max_abs_err > options->eps_stop) )
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
            int cur_thread, cur_klo_fb, cur_khi_fb, num_threads, nslices, i;
            int cur_klo_gb, cur_khi_gb;
                       
            cur_thread = omp_get_thread_num();
            num_threads = omp_get_num_threads();
    
            if (num_threads > 1) flag = 1;

            nslices = g->khi_fb - g->klo_fb + 1;

            cur_klo_fb = g->klo_fb + nslices*cur_thread/num_threads;
            cur_khi_fb = g->klo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            double t1 = omp_get_wtime();
            /* Keeping track of thread-local ghost boundaries, mainly for imposing mask */
            if (cur_khi_fb > (g->khi_fb))
	            cur_khi_fb = (g->khi_fb);

            if (cur_thread == 0)
                cur_klo_gb = cur_klo_fb - 3;
            else
                cur_klo_gb = cur_klo_fb;
            
            if (cur_thread == (num_threads - 1))
                cur_khi_gb = cur_khi_fb + 3;
            else
                cur_khi_gb = cur_khi_fb;
            

            while( dt_sub < options->tplot )
            { 
                /* inner loop */
                INNER_STEP++;
                      
                (options->space_deriv_func)(p->phi, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));	 
	                
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

	            QSS3D_RK1_STEP(p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);

                #pragma omp barrier
                
                /* boundary conditions */
                #pragma omp single
                {
                    signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	
                    cur_max_H_over_dX = -1;
                    max_H_over_dX = -1;
                } 
	            
	            /* phi_next should now store the RK_stage1 output */
	            if (options->order_time_accur >= 2) {
	                (options->space_deriv_func)(p->phi_next, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	                
	                #pragma omp barrier
	                
	                QSS3D_TVD_RK2_STAGE2(p->scratch1,GB_DIMS,p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);
		            
		            #pragma omp barrier
		            #pragma omp single
		            {
                        COPY_DATA(p->phi_next, p->scratch1, g);
                    }
		            /* boundary conditions */
                     #pragma omp single
                     {
                        signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);     
                     }    
                     	
                    /* scratch1 should now have RK_stage2 output. If third order is not desired, phi_next 
                       stores final output */
                }
                
                /* 3rd order is the highest time order accuracy handled */
                if (options->order_time_accur >= 3) {
	                (options->space_deriv_func)(p->phi_next, &a, &(options->b),
		            &(max_H_over_dX), p->lse_rhs, GB_DIMS, FB_DIMS_PAR,
	                &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
	                
	                #pragma omp barrier
	                
	                QSS3D_TVD_RK3_STAGE3(p->scratch2,GB_DIMS,p->phi_next,GB_DIMS,p->phi,GB_DIMS,p->lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt);
		        
		            #pragma omp barrier
		            #pragma omp single
		            {
                        COPY_DATA(p->phi_next, p->scratch2, g);
                    }
		            /* boundary conditions */
                     #pragma omp single
                    {
                        signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);
                    }
                     	
                    /* scratch3 should now have RK_stage3 output. phi_next stores final output */
                 }
		        
		         #pragma omp single
	             {
	                SET_DATA_TO_CONSTANT(p->lse_rhs,g,zero);
	                SET_DATA_TO_CONSTANT(p->scratch1,g,zero);
	                SET_DATA_TO_CONSTANT(p->scratch2,g,zero);
	                dt_sub += dt;
	             }  
		            
		        /* Since all final outputs are stored in phi_next, IMPOSE_MASK works for all time stepping methods */
	            #pragma omp barrier
	            IMPOSE_MASK_PAR(p->phi, p->mask, p->phi_next, &(options->overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	            if(options->check_connectivity) IMPOSE_MASK_PAR(p->phi, p->mask_w, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	            
	            #pragma omp barrier
	           if(options->check_connectivity) IMPOSE_MASK_PAR(p->phi, p->mask_nw, p->phi, &(disconn_overlap),
	              GB_DIMS, &(cur_klo_gb), &(cur_khi_gb));
	              
            	#pragma omp barrier
	            	    
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
        
        if(options->check_connectivity) {
            trapComponents_mask(p, g, options);
            reinitializeSubcellFix3d(p->mask_w,g,options);
            reinitializeSubcellFix3d(p->mask_nw,g,options);
            
             /* merge trapped nw phase for ensuring reinitialization, and correct saturations */
	        MERGE_SETS(p->phi, p->mask_nw, g);
        }
        
	    
        fprintf(fp_out, "Reinitializing....");
        //reinitialize3d_subcell_fix_qss(p,g,options);
        reinitializeSubcellFix3d(p->phi,g,options);
        fprintf(fp_out, "Reinitialized\n");

            
        /* compute stopping criteria */
        /* max abs error */
        QSS3D_MAX_NORM_DIFF_LOCAL(&max_abs_err,p->phi,GB_DIMS, p->phi_prev,
		        GB_DIMS, FB_DIMS, &(options->err_check_zone));

        QSSLIB_REAL vol_phi_old = vol_phi;
        
        satn_old = satn;
        
        QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&vol_phi, p->phi, GB_DIMS, FB_DIMS, 
		        &(g->dx[0]),&(g->dx[1]),&(g->dx[2]), &eps);
		
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

        fprintf(fp_out, "t = %4.1f\t",t);
        fprintf(fp_out, "max_abs_err = %4.3f,\t",max_abs_err);
        fprintf(fp_out, "nw phase satn pct = %4.3f\n", satn*100);
        
        printf("t = %4.1f\t",t);
        printf("max_abs_err = %4.3f,\t",max_abs_err);
        printf("nw phase satn = %4.3f\n", satn*100);
        fflush(fp_out);

        if (options->use_satn_stop) {
            if (fabsf(100*(satn - satn_old)) < 0.001)
                break;
        }
        
    } /* End outer loop */
    
    //MERGE_SETS(p->phi, p->mask_nw, g);
    
    return t;
}




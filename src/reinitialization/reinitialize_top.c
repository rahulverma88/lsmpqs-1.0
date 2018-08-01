/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file reinitialize_top.c
 *
 * Description: Function definitions for reinitialization algorithms. Both 2D and 3D.
 *      reinitializeMedium* implements the reinitialization algorithm utilizing the 
 *      Godunov discretization. The SubcellFix variant is supposed to be more accurate,
 *      with less mass losses.
 */

#include <stdlib.h>
#include <omp.h>

#include "QSSLIB_config.h"
#include "qss_general_util.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_tvd_runge_kutta2d.h"
#include "qss_macros.h"
#include "qss_reinitialization3d.h"
#include "qss_reinitialization2d.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_spatial_derivatives2d.h"
#include "qss_initialization2d.h"
#include "qss_initialization3d.h"

void reinitializeMedium2d(
     QSSLIB_REAL    *data,
     Grid           *grid,
     Options        *options)
{  
    QSSLIB_REAL cfl_number = 0.5;
    QSSLIB_REAL t_r, dt_r, tmax_r;
    
    QSSLIB_REAL *data0, *lse_rhs, *data_next, *data_stage1;
    
    
    int    use_data0_for_sign = 0;    
    int    bdry_location_idx = 9; /* all boundaries */
   
    /* writing shortcuts */
    Grid             *g = grid;
    Options          *o = options;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
    tmax_r = options->tmax_r;
    
     /* allocate memory for arrays */
    data0 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    lse_rhs = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    data_next = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    data_stage1 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    
    COPY_DATA(data0,data,g)
    
    #pragma omp parallel default(none) shared(data, lse_rhs, data_next, g, cfl_number, \
        options, dt_r, data0, data_stage1, t_r, bdry_location_idx, use_data0_for_sign, tmax_r) 
    {    
        /* Set up variables for multi-threading */
        int cur_thread, cur_jlo_fb, cur_jhi_fb, num_threads, nslices, i;
                       
        cur_thread = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        nslices = g->jhi_fb - g->jlo_fb + 1;

        cur_jlo_fb = g->jlo_fb + nslices*cur_thread/num_threads;
        cur_jhi_fb = g->jlo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            
        while(t_r < tmax_r )
        {
            COMPUTE_REINITIALIZATION_EQN_RHS_2D(lse_rhs, data, data0, GB_DIMS_2D, 
                FB_DIMS_PAR_2D, &((g->dx)[0]), &((g->dx)[1]), &use_data0_for_sign);
       
	        #pragma omp barrier

            QSS2D_RK1_STEP(data_stage1,GB_DIMS_2D,data,GB_DIMS_2D,lse_rhs,
		        GB_DIMS_2D, FB_DIMS_PAR_2D, &dt_r);
		    
	        #pragma omp barrier		        
            /* boundary conditions */
            #pragma omp single
            {
                signedLinearExtrapolationBCqss(data_stage1,g,bdry_location_idx);
            }	 	 	    
      
            COMPUTE_REINITIALIZATION_EQN_RHS_2D(lse_rhs, data_stage1, data0, GB_DIMS_2D, 
                FB_DIMS_PAR_2D, &((g->dx)[0]), &((g->dx)[1]), &use_data0_for_sign);
        
	        #pragma omp barrier

            QSS2D_TVD_RK2_STAGE2(data_next,GB_DIMS_2D,data_stage1,GB_DIMS_2D,data,GB_DIMS_2D,lse_rhs,
		        GB_DIMS_2D, FB_DIMS_PAR_2D, &dt_r);

		    #pragma omp barrier

   	   
            /* boundary conditions */
            #pragma omp single
            {
                signedLinearExtrapolationBCqss(data_next,g,bdry_location_idx); 
            }
            
            #pragma omp single
            {
                COPY_DATA(data,data_next,g);
            }
       
            #pragma omp single
            {
                t_r = t_r + dt_r;   
            }
        }
    }
    
    free(lse_rhs); free(data0); free(data_stage1); free(data_next);
}	 

void reinitializeSubcellFix2d(
    QSSLIB_REAL *data,
    Grid *g,
    Options *o)
{  

    QSSLIB_REAL   *distance0 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ); //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    QSSLIB_REAL   *lse_rhs = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    QSSLIB_REAL   *data_next = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,data,g);
    
    QSS2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS_2D, FB_DIMS_2D,
		    &((g->dx)[0]), &((g->dx)[1]));
		    
    #pragma omp parallel default(none) shared(data, lse_rhs, data_next, g, cfl_number, \
        o, dt_r, distance0, copy, t_r, bdry_location_idx, n_steps) 
    {    
        /* Set up variables for multi-threading */
        int cur_thread, cur_jlo_fb, cur_jhi_fb, num_threads, nslices, i;
        int cur_jlo_gb, cur_jhi_gb;
                       
        cur_thread = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        nslices = g->jhi_fb - g->jlo_fb + 1;

        cur_jlo_fb = g->jlo_fb + nslices*cur_thread/num_threads;
        cur_jhi_fb = g->jlo_fb + nslices*(cur_thread + 1)/num_threads - 1;
		 
		#pragma omp barrier
		
        while(t_r < o->tmax_r )
        {
            #pragma omp single
            n_steps++;
 
            QSS2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(lse_rhs,
		        data, GB_DIMS_2D, copy, distance0, FB_DIMS_PAR_2D,
		        &((g->dx)[0]), &((g->dx)[1]));
		 
	        #pragma omp barrier

            QSS2D_RK1_STEP(data_next, GB_DIMS_2D, data, GB_DIMS_2D, lse_rhs,
                GB_DIMS_2D, FB_DIMS_PAR_2D, &dt_r);
                
	        #pragma omp barrier
		   	 
	        #pragma omp single
	        {
             signedLinearExtrapolationBCqss(data_next,g,bdry_location_idx);	 	 
            }
      
            #pragma omp single
            {
                COPY_DATA(data,data_next,g);	
            }
        
            #pragma omp single
            {
	            SET_DATA_TO_CONSTANT(lse_rhs,g,0); 
	        }
       
            #pragma omp single
            {
                t_r = t_r + dt_r;   
            }
        }
     }
     free(distance0); free(copy); free(lse_rhs); free(data_next);
}


void reinitializeMedium3d(
     QSSLIB_REAL    *data,
     Grid           *grid,
     Options        *options)
{  
    QSSLIB_REAL cfl_number = 0.5;
    QSSLIB_REAL t_r, dt_r, tmax_r;
    
    QSSLIB_REAL *data0, *lse_rhs, *data_next, *data_stage1;
    
    int    use_data0_for_sign = 0;    
    int    bdry_location_idx = 9; /* all boundaries */
   
      /* writing shortcuts */
    Grid             *g = grid;
    Options          *o = options;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
    tmax_r = options->tmax_r;
    
     /* allocate memory for arrays */   
    data0 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    lse_rhs = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    data_next = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    data_stage1 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    
    COPY_DATA(data0, data, g);

    #pragma omp parallel default(none) shared(data, lse_rhs, data_next, g, cfl_number, \
        options, dt_r, data0, data_stage1, t_r, bdry_location_idx, use_data0_for_sign, tmax_r) 
    {    
        int cur_thread, cur_klo_fb, cur_khi_fb, num_threads, nslices, i;
        int cur_klo_gb, cur_khi_gb;
                       
        cur_thread = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        nslices = g->khi_fb - g->klo_fb + 1;

        cur_klo_fb = g->klo_fb + nslices*cur_thread/num_threads;
        cur_khi_fb = g->klo_fb + nslices*(cur_thread + 1)/num_threads - 1;
            
            
        while(t_r < tmax_r )
        {
            COMPUTE_REINITIALIZATION_EQN_RHS_3D(lse_rhs, data, data0, GB_DIMS, 
                FB_DIMS_PAR, &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]), &use_data0_for_sign);
        
            #pragma omp barrier
            
            QSS3D_RK1_STEP(data_stage1,GB_DIMS,data,GB_DIMS,lse_rhs,
		            GB_DIMS, FB_DIMS_PAR, &dt_r);
		            
		    #pragma omp barrier
		        
            #pragma omp single
            {
             signedLinearExtrapolationBCqss(data_stage1,g,bdry_location_idx);	 	 	   
            } 
      
            COMPUTE_REINITIALIZATION_EQN_RHS_3D(lse_rhs, data_stage1, data0, GB_DIMS, 
                FB_DIMS_PAR, &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]), &use_data0_for_sign);
            
            #pragma omp barrier
            
            QSS3D_TVD_RK2_STAGE2(data_next,GB_DIMS,data_stage1,GB_DIMS,data,GB_DIMS,lse_rhs,
		        GB_DIMS, FB_DIMS_PAR, &dt_r);
   	   
   	        #pragma omp barrier
   	        
            #pragma omp single 
            {
                signedLinearExtrapolationBCqss(data_next,g,bdry_location_idx); 
            }
            
            #pragma omp single
            {
                COPY_DATA(data,data_next,g);
            }
            
            #pragma omp single
            {
                t_r = t_r + dt_r;   
            }
                
            #pragma omp barrier
        }
    }
        
    free(lse_rhs); free(data0); free(data_stage1); free(data_next);
}	 

void reinitializeSubcellFix3d(
    QSSLIB_REAL *data,
    Grid *g,
    Options *o)
{  

    QSSLIB_REAL   *distance0 = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ); //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    QSSLIB_REAL   *lse_rhs = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    QSSLIB_REAL   *data_next = (QSSLIB_REAL*) malloc(g->num_gridpts * DSZ);
    
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,data,g);
    
    QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS, FB_DIMS,
		    &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		    
    #pragma omp parallel default(none) shared(data, lse_rhs, data_next, g, \
        o, dt_r, distance0, copy, t_r, bdry_location_idx, n_steps) 
    {    
        /* Set up variables for multi-threading */
        int cur_thread, cur_klo_fb, cur_khi_fb, num_threads, nslices, i;
        int cur_klo_gb, cur_khi_gb;
        cur_thread = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        nslices = g->khi_fb - g->klo_fb + 1;
        cur_klo_fb = g->klo_fb + nslices*cur_thread/num_threads;
        cur_khi_fb = g->klo_fb + nslices*(cur_thread + 1)/num_threads - 1;
                    
		         
        #pragma omp barrier
        
        while(t_r < o->tmax_r )
        {
            #pragma omp single
            n_steps++;
 
            QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(lse_rhs,
		        data, GB_DIMS, copy, distance0, FB_DIMS_PAR,
		        &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
	        #pragma omp barrier
	        
            QSS3D_RK1_STEP(data_next, GB_DIMS, data, GB_DIMS, lse_rhs,
                GB_DIMS, FB_DIMS_PAR, &dt_r);
		   	
		   	#pragma omp barrier
		   	
	        #pragma omp single
	        {
                signedLinearExtrapolationBCqss(data_next,g,bdry_location_idx);	 
            }	 
      
            #pragma omp single
            {
                COPY_DATA(data,data_next,g);
            }
            	 
            #pragma omp single
	        {
	         SET_DATA_TO_CONSTANT(lse_rhs,g,0);
	        }
	        
            #pragma omp single
            {
                t_r = t_r + dt_r;   
            }
        }
        
        
     }
     free(distance0); free(copy); free(lse_rhs); free(data_next);
}

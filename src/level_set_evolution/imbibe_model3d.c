/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file imbibe_model3d.c

    Top level function for calling the imbibition constant curvature model in 3D. 
    
    This function does several things:
    
    1. Initializes and re-initializes all masks (solid space, disconnected components).
    2. Sets up the contact angle model - for constant theta, it sets up the variational a
        and b matrices. For spatially varying theta, it reads in a "theta.gz" file which 
        provides values of theta at different points in the domain. Based on the theta
        information, this function then determines overlap values at different points in 
        the domain between the pore space and the solid space.
    3. Calculates gradients of the solid space mask- this will be required in the level set
        computations, and needs to be only computed once since our mask doesn't change.
    4. For zero theta, it calls the less-memory intensive constCurvModel3dNoVar, while for
        all other cases, constCurvModel3d is called.
    5. It updates the "a" term for the next step.
    6. It writes out the data steps for each step.
             
*/

/* System headers */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* Local headers */
#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_general_util.h"
#include "qss_reinitialization3d.h"
#include "qss_grid.h"
#include "qss_macros.h"
#include "constCurvModel3d.h"
#include "connectivity.h"

/* Main driver for constant curvature level set method model */

void imbibe_model3d(Options *options,QSS_DataArrays *p, Grid *g, FILE *fp_out) 
{  
    char fname[256];
    int n1[3];
    QSSLIB_REAL a0 = options->a;
    QSSLIB_REAL da = options->b * options->dc;
    QSSLIB_REAL inv_dx, inv_dy, inv_dz;
    int step = options->init_step + 1;
    QSSLIB_REAL t_step = 0;
    
    /* Initialize disconnected component masks */
    initializeDisconnectedMasks(p->mask_w, g);
    initializeDisconnectedMasks(p->mask_nw, g);
    
    COPY_DATA(p->mask_disconn_init, p->mask_nw, g);
    
    /* Create reservoir inlet */
    createReservoirInlet3d(p, g);
    
    /* Set amin */
    QSSLIB_REAL amin = options->amin;

    /* Compute gradients of mask */
    QSS3D_CENTRAL_GRAD_ORDER2(p->mask_x, p->mask_y, p->mask_z,
        GB_DIMS, p->mask, GB_DIMS, FB_DIMS,
        &((g->dx)[0]),&((g->dx)[1]),&((g->dx)[2]));
    
    if (options->use_var_theta) {
        sprintf(fname, "theta.gz");
        p->theta = readDataArrayQSS(n1, fname);
    }
    
    setSpaceDerivFunc(g,options);
    
    /* Set overlap:
     For 3d, if theta is constant, overlap is constant.
     if theta is variational, then theta is read in and 
     overlap is set according to local contact angle
    */
    setThetaOverlap(p, g, options);
    
    /* Impose reservoir inlet, if in options */
    if (options->reservoir_inlet) {
        IMPOSE_MIN(p->phi, p->phi, p->phi_extra, g);
        IMPOSE_MASK(p->phi, p->mask, p->phi, g);
    }
    
    /* Get inlet index and outlet index*/
    getMainInd_new(options, p, g);
    
    while (a0 > amin) 
    { 
        printf("\nStep %d: a = %f\n", step, a0);
        fprintf(fp_out, "\nStep %d: a = %f\n", step, a0);
        
        /* Compute variational curvature term b and advective term vel */ 
        if ( (options->theta > 0) || (options->use_var_theta) )   
            QSS3D_SET_VAR_CURV_ADV(p->mask, p->mask_x, p->mask_y, p->mask_z,
     	        p->curvature_coeff, p->external_velocity_x, p->external_velocity_y,
     	        p->external_velocity_z, &(options->b_max_over_dx), &(options->max_U_over_dx),
                &(options->b), GB_DIMS, FB_DIMS, &(g->dx[0]));
        else {
            inv_dx = 1/g->dx[0];
            inv_dy = 1/g->dx[1];
            inv_dz = 1/g->dx[2];
            options->b_max_over_dx = 2 * options->b * 
                    (inv_dx*inv_dx + inv_dy*inv_dy + inv_dz*inv_dz);
            options->max_U_over_dx = 0;
        }

        options->conserve_imbibe = 1;
        
        /* compute variational a. variational b and vel are already computed. */
        if ( (options->theta > 0) || (options->use_var_theta) ) {
            if ( ~(options->use_var_theta) ) {
                QSS3D_SET_VAR_NORM(p->mask, p->mask_x, p->mask_y, p->mask_z,
     	            p->normal_velocity, &(a0), &(options->theta), GB_DIMS, FB_DIMS, &(g->dx[0]));
     	    } else {
     	        QSS3D_SET_VAR_NORM_THETA(p->mask, p->mask_x, p->mask_y, p->mask_z,
     	            p->normal_velocity, &(a0), p->theta, GB_DIMS, FB_DIMS, &(g->dx[0]));
     	    }
     	} 
     	
     	if(options->check_connectivity) trapComponents_mask(p, g, options);
     	
        if ( (options->theta > 0) || (options->use_var_theta) )   
            t_step = constCurvModel3d(options, p, g, fp_out);
        else
            t_step = constCurvModel3dNoVar(options, p, g, fp_out, a0);
        
        if (t_step == -1)
            break;
            
        a0 -= da;
        
        /* Impose reservoir inlet, if in options */
        if (options->reservoir_inlet) {
            IMPOSE_MIN(p->phi, p->phi, p->phi_extra, g);
        }
        
        /* Impose original mask for saving data. Masks imposed inside will be different for theta > 30 degrees */
        IMPOSE_MASK(p->phi, p->mask, p->phi, g);
        
        //if (t_step < options->tmax)
        sprintf(fname,"data_step_%d",step);
	    //else
		  //  sprintf(fname,"data_step_%d_no_eq",step);

	    step += 1;
	    writeDataArrayQSS(p->phi,g,fname,GZIP);
    }
      
}



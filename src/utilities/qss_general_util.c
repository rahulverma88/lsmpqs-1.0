#include <stdlib.h>

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

void reinitialize3d_subcell_fix_qss(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{  

    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,p->phi,g);
     		    		    
    QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->phi, GB_DIMS, copy, distance0, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
       QSS3D_RK1_STEP(p->phi_next, GB_DIMS, p->phi, GB_DIMS,p->lse_rhs,
         GB_DIMS, FB_DIMS, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->phi,p->phi_next,g);	 
       
      t_r = t_r + dt_r;   
     }
     
     //free(distance0);
     //free(copy);
}

void setSpaceDerivFunc(Grid *g, Options *options)
{
    if (g->num_dims == 2) {
            switch(options->order_space_accur) {
            case 2: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_2D_SECOND_CONST;
                    else
                        options->space_deriv_func = GET_RHS_2D_SECOND;
                    break;
            case 3: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_2D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_2D_THIRD;
                    break;
            case 4: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_2D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_2D_THIRD;
                    break;
            case 5: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_2D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_2D_THIRD;
                    break;
            default: options->space_deriv_func = GET_RHS_2D_THIRD;
        }
    }
    else {
    
        switch(options->order_space_accur) {
            case 2: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_3D_SECOND_CONST;
                    else
                        options->space_deriv_func = GET_RHS_3D_SECOND;
                    break;
            case 3: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_3D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_3D_THIRD;
                    break;
            case 4: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_3D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_3D_THIRD;
                    break;
            case 5: 
                    if ( (options->theta == 0) && (~options->use_var_theta) )
                        options->space_deriv_func = GET_RHS_3D_THIRD_CONST;
                    else
                        options->space_deriv_func = GET_RHS_3D_FIFTH;
                    break;
            default: options->space_deriv_func = GET_RHS_3D_THIRD;
        }
    }
}

void setThetaOverlap(QSS_DataArrays *p, Grid *g, Options *options)
{
    
    QSSLIB_REAL theta, overlap;
    if (~options->use_var_theta) theta = options->theta;
    
    int i, j, k, idx_gb, nx_gb, nxy_gb;
    
    nx_gb = g->grid_dims_ghostbox[0]; nxy_gb = nx_gb*(g->grid_dims_ghostbox[1]);
    
    QSSLIB_REAL rad_factor = 3.14159/180;
    
    if (g->num_dims == 3) {
    if (~options->use_var_theta)
    {
        if (theta <= 20*rad_factor)
            options->overlap = 0;
        else if (theta <= 40*rad_factor)
            options->overlap = 0.3 * g->dx[0];
        else if (theta <= 50*rad_factor)
            options->overlap = 0.5 * g->dx[0];
        else if (theta > 50*rad_factor)
            options->overlap = 1 * g->dx[0];
        else if (theta > 90*rad_factor)
            options->overlap = 2 * g->dx[0];
    } else {
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                for( i = g->ilo_fb; i <= g->ihi_fb; i++)
                {
                    idx_gb = i + j*nx_gb + k*nxy_gb;
                    theta = p->theta[idx_gb];
                    if (theta <= 20*rad_factor)
                        overlap = 0;
                    else if (theta <= 40*rad_factor)
                        overlap = 0.3 * g->dx[0];
                    else if (theta <= 50*rad_factor)
                        overlap = 0.5 * g->dx[0];
                    else if (theta > 50*rad_factor)
                        overlap = 1 * g->dx[0];
                    else if (theta > 90*rad_factor)
                        overlap = 2 * g->dx[0];
                        
                    p->overlap[idx_gb] = overlap;
                }
            }
        }                   
    }
    } else {
        options->overlap = 3*g->dx[0];
    }
        
}
void reinitialize2d_subcell_fix_qss(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{  

    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,p->phi,g);
     		    		    
    QSS2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS_2D, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->phi, GB_DIMS_2D, copy, distance0, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
       QSS2D_RK1_STEP(p->phi_next, GB_DIMS_2D, p->phi, GB_DIMS_2D,p->lse_rhs,
         GB_DIMS_2D, FB_DIMS_2D, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->phi,p->phi_next,g);	 
       
      t_r = t_r + dt_r;   
     }
     
     //free(distance0);
     //free(copy);
}

void qss_reinitialize_mask(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{
    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,p->mask,g)
     		    		    
    QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->mask, GB_DIMS, copy, distance0, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
       QSS3D_RK1_STEP(p->phi_next, GB_DIMS, p->mask, GB_DIMS,p->lse_rhs,
         GB_DIMS, FB_DIMS, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->mask,p->phi_next,g)	 
       
      t_r = t_r + dt_r;   
     }
     
     //free(distance0);
     //free(copy);
     
     return;

}

void qss_reinitialize_mask_no_bc(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{
    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx_1 = 7; 
    int    bdry_location_idx_2 = 8;
    
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,p->mask,g)
     		    		    
    QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->mask, GB_DIMS, copy, distance0, FB_DIMS,
		 &((g->dx)[0]), &((g->dx)[1]),&((g->dx)[2]));
		 
       QSS3D_RK1_STEP(p->phi_next, GB_DIMS, p->mask, GB_DIMS,p->lse_rhs,
         GB_DIMS, FB_DIMS, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx_1);	 	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx_2);
      
      COPY_DATA(p->mask,p->phi_next,g)	 
       
      t_r = t_r + dt_r;   
     }
     
     //free(distance0);
     //free(copy);
     
     return;

}

void qss_reinitialize_mask2d(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{
    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
      
    COPY_DATA(copy,p->mask,g)
     		    		    
    QSS2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS_2D, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->mask, GB_DIMS_2D, copy, distance0, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
       QSS2D_RK1_STEP(p->phi_next, GB_DIMS_2D, p->mask, GB_DIMS_2D,p->lse_rhs,
         GB_DIMS_2D, FB_DIMS_2D, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->mask,p->phi_next,g)	 
       
      t_r = t_r + dt_r;   
     }
     
     //free(distance0);
     //free(copy);
     
     return;

}



void signedLinearExtrapolationBCqss(
  QSSLIB_REAL *phi,
  Grid *grid,
  int bdry_location_idx)
{
  int num_dims = grid->num_dims;
  if (num_dims == 2) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: {
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        QSS2D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } else if (num_dims == 3) {

    switch (bdry_location_idx) { 
      case 0: 
      case 1: 
      case 2: 
      case 3: 
      case 4: 
      case 5: {
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &bdry_location_idx);
        break;
      }
      case 6: {
        int tmp_bdry_location_idx = 0;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 7: {
        int tmp_bdry_location_idx = 2;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 8: {
        int tmp_bdry_location_idx = 4;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      case 9: {
        int tmp_bdry_location_idx = 0;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 1;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 2;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 3;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 4;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        tmp_bdry_location_idx = 5;
        QSS3D_SIGNED_LINEAR_EXTRAPOLATION(
          phi,
          &(grid->ilo_gb), &(grid->ihi_gb), 
          &(grid->jlo_gb), &(grid->jhi_gb), 
          &(grid->klo_gb), &(grid->khi_gb), 
          &(grid->ilo_fb), &(grid->ihi_fb), 
          &(grid->jlo_fb), &(grid->jhi_fb), 
          &(grid->klo_fb), &(grid->khi_fb), 
          &tmp_bdry_location_idx);

        break;
      };

      default: { /* DO NOTHING */ }

    }

  } /* end switch on num_dims */

}

void initializeDisconnectedMasks(QSSLIB_REAL *data, Grid *g)
{
    QSSLIB_REAL    center_x, center_y, center_z, radius;

    if (g->num_dims == 2){
        center_y = 0.5*(g->x_lo[1] + g->x_hi[1]);
        center_x = 0.5*(g->x_lo[0] + g->x_hi[0]);
        radius = 2*(g->x_hi[0] - g->x_lo[0]);
        createCircle(data, center_x, center_y, radius, -1, g);
    } else if(g->num_dims == 3){
        center_z = 0.5*(g->x_lo[2] + g->x_hi[2]);
        center_y = 0.5*(g->x_lo[1] + g->x_hi[1]);
        center_x = 0.5*(g->x_lo[0] + g->x_hi[0]);
        radius = 2*(g->x_hi[0] - g->x_lo[0]);
        createSphere(data, center_x, center_y, center_z, radius, -1, g);
    }
        
}
 
void qss_reinitializeDisconnectedMask2d(
    QSS_DataArrays *p,
    Grid *g,
    Options *o)
{
    QSSLIB_REAL   *distance0 = p->scratch1; //if doing 2nd order subcell fix, will need a different array
    QSSLIB_REAL   *copy      = p->scratch2;
     
    QSSLIB_REAL cfl_number = 0.5; // watch out, 0.9 might not do well in some cases
    QSSLIB_REAL t_r, dt_r;
        
    int    bdry_location_idx = 9; /* all boundaries */
    int    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
    
    /* Do wetting mask */
    COPY_DATA(copy,p->mask_w,g)
     		    		    
    QSS2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS_2D, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->mask_w, GB_DIMS_2D, copy, distance0, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
       QSS2D_RK1_STEP(p->phi_next, GB_DIMS_2D, p->mask_w, GB_DIMS_2D,p->lse_rhs,
         GB_DIMS_2D, FB_DIMS_2D, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->mask_w,p->phi_next,g)	 
       
      t_r = t_r + dt_r;   
     }
     
     /* Do non-wetting mask */
    n_steps = 0;
   
    t_r = 0;
    dt_r = cfl_number * (g->dx)[0];
    COPY_DATA(copy,p->mask_nw,g)
     		    		    
    QSS2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(distance0, copy, GB_DIMS_2D, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
    while(t_r < o->tmax_r )
    {
        n_steps++;
 
       QSS2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(p->lse_rhs,
		 p->mask_nw, GB_DIMS_2D, copy, distance0, FB_DIMS_2D,
		 &((g->dx)[0]), &((g->dx)[1]));
		 
       QSS2D_RK1_STEP(p->phi_next, GB_DIMS_2D, p->mask_nw, GB_DIMS_2D,p->lse_rhs,
         GB_DIMS_2D, FB_DIMS_2D, &dt_r);
		   	 
      signedLinearExtrapolationBCqss(p->phi_next,g,bdry_location_idx);	 	 
      
      COPY_DATA(p->mask_nw,p->phi_next,g)	 
       
      t_r = t_r + dt_r;   
     }
     //free(distance0);
     //free(copy);
     
     return;

}

void createReservoirInlet3d(
    QSS_DataArrays *p,
    Grid *g) 
{

  QSSLIB_REAL corner_x, corner_y, corner_z;
  QSSLIB_REAL side_length_x, side_length_y, side_length_z;
  
  corner_x = g->x_lo_ghostbox[0];
  corner_y = g->x_lo_ghostbox[1];
  corner_z = g->x_lo_ghostbox[2];
  
  side_length_x = 3*g->dx[0];
  side_length_y = g->x_hi_ghostbox[1] - g->x_lo_ghostbox[1];
  side_length_z = g->x_hi_ghostbox[2] - g->x_lo_ghostbox[2];
  
  createBox(p->phi_extra, corner_x, corner_y, corner_z,
    side_length_x, side_length_y, side_length_z,
    -1, g);

}
    

int return_1(){ return 1;}
int return_0(){return 0;}
double return_0_double(){return 0;};

void readSpheresBinaryFile(
    char *file_name,
    int  *pn_sphere,
    QSSLIB_REAL *x_lo,
    QSSLIB_REAL *x_hi,
    QSSLIB_REAL **pcenterx,
    QSSLIB_REAL **pcentery,
    QSSLIB_REAL **pcenterz,
    QSSLIB_REAL **pradius)
{    
    int   zip_status;
    char   *file_base;
    FILE *fp;
    
    int    n_sphere;
    QSSLIB_REAL   *centerx, *centery, *centerz, *radius;
    
    checkUnzipFile(file_name,&zip_status,&file_base);    
    fp = fopen(file_base,"r");
    
    if( fp )
    {
       fread(&n_sphere,sizeof(int),1,fp);
       fread(x_lo,DSZ,3,fp);
       fread(x_hi,DSZ,3,fp);

       centerx = (QSSLIB_REAL *)malloc(n_sphere*DSZ);
       centery = (QSSLIB_REAL *)malloc(n_sphere*DSZ);
       centerz = (QSSLIB_REAL *)malloc(n_sphere*DSZ);
       radius  = (QSSLIB_REAL *)malloc(n_sphere*DSZ);

       fread(centerx,DSZ,n_sphere,fp);
       fread(centery,DSZ,n_sphere,fp);
       fread(centerz,DSZ,n_sphere,fp);
       fread(radius,DSZ,n_sphere,fp);

       fclose(fp);
       zipFile(file_base,zip_status);
    }
    else
    {
        n_sphere = 0;
        printf("\nCould not open file %s\n",file_name);
    }
   
    free(file_base); 
    *pn_sphere = n_sphere;
    *pcenterx = centerx;
    *pcentery = centery;
    *pcenterz = centerz;
    *pradius = radius;
}

void writeSpheresBinaryFile(
    char *file_base, 
    int zip_status,
    int  n_sphere,
    QSSLIB_REAL *x_lo,
    QSSLIB_REAL *x_hi,
    QSSLIB_REAL *centerx,
    QSSLIB_REAL *centery,
    QSSLIB_REAL *centerz,
    QSSLIB_REAL *radius)
{    
    FILE *fp;   
    
    fp = fopen(file_base,"w");
    
    if( fp )
    {
       fwrite(&n_sphere,sizeof(int),1,fp);
       fwrite(x_lo,DSZ,3,fp);
       fwrite(x_hi,DSZ,3,fp);

       fwrite(centerx,DSZ,n_sphere,fp);
       fwrite(centery,DSZ,n_sphere,fp);
       fwrite(centerz,DSZ,n_sphere,fp);
       fwrite(radius,DSZ,n_sphere,fp);

       fclose(fp);
       zipFile(file_base,zip_status);
    }
    else
    {
        printf("\nCould not open file %s\n",file_base);
    }
}


/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file qss_initialization2d.c

    Function definitions for 2D initialization.
             
*/
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "qss_initialization3d.h"
#include "qss_macros.h"


void createPlane(
  QSSLIB_REAL *phi, 
  QSSLIB_REAL normal_x, QSSLIB_REAL normal_y, QSSLIB_REAL normal_z,
  QSSLIB_REAL point_x, QSSLIB_REAL point_y, QSSLIB_REAL point_z,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces3d() */
  createIntersectionOfHalfSpaces3d(
    phi, 1, 
    &normal_x, &normal_y, &normal_z,
    &point_x, &point_y, &point_z, 
    grid);
}


void createIntersectionOfHalfSpaces3d(
  QSSLIB_REAL *phi, 
  int num_half_spaces,
  QSSLIB_REAL *normal_x, QSSLIB_REAL *normal_y, QSSLIB_REAL *normal_z,
  QSSLIB_REAL *point_x, QSSLIB_REAL *point_y, QSSLIB_REAL *point_z,
  Grid *grid)
{
  int    i, j, k, l;
  int    idx, nxy;
  QSSLIB_REAL x, y, z;
  QSSLIB_REAL dot_prod, norm;
  QSSLIB_REAL max;
  QSSLIB_REAL signed_dist_to_plane;
  
  nxy = grid->grid_dims_ghostbox[0]*grid->grid_dims_ghostbox[1];
  for (k = 0; k < grid->grid_dims_ghostbox[2]; k++)
  {
    for (j = 0; j < grid->grid_dims_ghostbox[1]; j++)
    {
      for (i = 0; i < grid->grid_dims_ghostbox[0]; i++) 
      {
        idx = i+j*grid->grid_dims_ghostbox[0] + k*nxy;
        x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
        y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
        z = (grid->x_lo_ghostbox)[2] + (grid->dx)[2]*k;
        max = -FLT_MAX;

        for (l = 0; l < num_half_spaces; l++)
        { 
          dot_prod = (x - point_x[l])*normal_x[l] + 
                     (y - point_y[l])*normal_y[l] +
                     (z - point_z[l])*normal_z[l];
          norm = normal_x[l]*normal_x[l] + normal_y[l]*normal_y[l]+ 
                 normal_z[l]*normal_z[l];
          norm = sqrt(norm);

          signed_dist_to_plane = dot_prod/norm;
          if (signed_dist_to_plane > max) max = signed_dist_to_plane;
        }

        phi[idx] = max;  

      }
    }
  } /* end loop over grid */ 
}


void createPolyhedron3d(
  QSSLIB_REAL *phi, 
  int num_planes,
  QSSLIB_REAL *normal_x, QSSLIB_REAL *normal_y, QSSLIB_REAL *normal_z,
  QSSLIB_REAL *point_x, QSSLIB_REAL *point_y, QSSLIB_REAL *point_z,
  Grid *grid)
{
  /* use createIntersectionOfHalfSpaces3d() */
  createIntersectionOfHalfSpaces3d(
    phi, num_planes, 
    normal_x, normal_y, normal_z,
    point_x, point_y, point_z, 
    grid);
}


void createSphere(
  QSSLIB_REAL   *phi,
  QSSLIB_REAL   center_x,
  QSSLIB_REAL   center_y,
  QSSLIB_REAL   center_z,
  QSSLIB_REAL   radius,
  int      inside_flag,
  Grid     *grid)
{
  /* use createIntersectionOfSpheres() */
  createIntersectionOfSpheres(
    phi, 1,
    &center_x, &center_y, &center_z,
    &radius, &inside_flag,
    grid);
}


void createIntersectionOfSpheres(
  QSSLIB_REAL   *phi,
  int       num_spheres,
  QSSLIB_REAL   *center_x,
  QSSLIB_REAL   *center_y,
  QSSLIB_REAL   *center_z,
  QSSLIB_REAL   *radius,
  int      *inside_flag,
  Grid     *grid)
{
  int    i, j, k, l;
  int    idx, nx, nxy;
  QSSLIB_REAL x, y, z;
  QSSLIB_REAL max;
  QSSLIB_REAL signed_dist_to_sphere;
  
  nx = (grid->grid_dims_ghostbox)[0];
  nxy = (grid->grid_dims_ghostbox)[0]*(grid->grid_dims_ghostbox)[1];
  
  for (k = 0; k < (grid->grid_dims_ghostbox)[2]; k++)
  {
    for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
    {
      for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
      {
        idx = i + j*nx + k*nxy;
        x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
        y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
        z = (grid->x_lo_ghostbox)[2] + (grid->dx)[2]*k;
    
        max = -FLT_MAX;

        for (l = 0; l < num_spheres; l++)
        {
          signed_dist_to_sphere = sqrt( (x-center_x[l])*(x-center_x[l])
                                       +(y-center_y[l])*(y-center_y[l])
                                       +(z-center_z[l])*(z-center_z[l]) ) 
                                - radius[l];

          if (inside_flag[l] >= 0) {
            signed_dist_to_sphere = -signed_dist_to_sphere;
          }

          if (signed_dist_to_sphere > max) max = signed_dist_to_sphere;
        }
        phi[idx] = max;
      }
    }
  } /* end of loop over grid */
}


void createCylinder(
  QSSLIB_REAL *phi,
  QSSLIB_REAL  tangent_x,
  QSSLIB_REAL  tangent_y,
  QSSLIB_REAL  tangent_z,
  QSSLIB_REAL  point_x,
  QSSLIB_REAL  point_y,
  QSSLIB_REAL  point_z,
  QSSLIB_REAL  radius,
  int     inside_flag,
  Grid    *grid)
{
  /* use createIntersectionOfCylinders() */
  createIntersectionOfCylinders(
    phi, 1,
    &tangent_x, &tangent_y, &tangent_z,
    &point_x, &point_y, &point_z,
    &radius,
    &inside_flag,
    grid);
}


void createIntersectionOfCylinders(
  QSSLIB_REAL   *phi,
  int       num_cylinders,
  QSSLIB_REAL   *tangent_x,
  QSSLIB_REAL   *tangent_y,
  QSSLIB_REAL   *tangent_z,
  QSSLIB_REAL   *point_x,
  QSSLIB_REAL   *point_y,
  QSSLIB_REAL   *point_z,
  QSSLIB_REAL   *radius,
  int      *inside_flag,
  Grid     *grid)
{
  int     i, j, k, l, nx, nxy;  
  
  nx = (grid->grid_dims_ghostbox)[0];
  nxy = (grid->grid_dims_ghostbox)[0]*(grid->grid_dims_ghostbox)[1];
  
  for (k = 0; k < (grid->grid_dims_ghostbox)[2]; k++)
  {
    for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
    {
      for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
      {
        QSSLIB_REAL  x, y, z, max;
        int idx = i+j*nx + k*nxy;

        x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
        y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
        z = (grid->x_lo_ghostbox)[2] + (grid->dx)[2]*k;
    
        max = -FLT_MAX;
        for (l = 0; l < num_cylinders; l++)
        {
          QSSLIB_REAL signed_dist_to_cylinder;
          QSSLIB_REAL norm_sq_x_minus_p;
          QSSLIB_REAL x_minus_p_dot_tangent;

          norm_sq_x_minus_p = (x - point_x[l])*(x - point_x[l])
                            + (y - point_y[l])*(y - point_y[l])
                            + (z - point_z[l])*(z - point_z[l]);

          x_minus_p_dot_tangent = (x - point_x[l])*tangent_x[l]
                                + (y - point_y[l])*tangent_y[l]
                                + (z - point_z[l])*tangent_z[l];
          x_minus_p_dot_tangent /= sqrt( tangent_x[l]*tangent_x[l]
                                       + tangent_y[l]*tangent_y[l]
                                       + tangent_z[l]*tangent_z[l]);

          signed_dist_to_cylinder = 
            sqrt( norm_sq_x_minus_p 
                - x_minus_p_dot_tangent*x_minus_p_dot_tangent) - radius[l];

          if (inside_flag[l] >= 0) {
            signed_dist_to_cylinder = -signed_dist_to_cylinder;
          }

          if (signed_dist_to_cylinder > max) max = signed_dist_to_cylinder; 
        }
        phi[idx] = max;
      }
    }
  } /* end of loop over grid */

}


void createHyperboloid(
  QSSLIB_REAL *phi,
  QSSLIB_REAL  tangent_x,
  QSSLIB_REAL  tangent_y,
  QSSLIB_REAL  tangent_z,
  QSSLIB_REAL  center_x,
  QSSLIB_REAL  center_y,
  QSSLIB_REAL  center_z,
  QSSLIB_REAL  alpha,
  QSSLIB_REAL  beta,
  int     inside_flag,
  Grid    *grid)
{
  /* use createIntersectionOfHyperboloids() */
  createIntersectionOfHyperboloids(
    phi, 1,
    &tangent_x, &tangent_y, &tangent_z,
    &center_x, &center_y, &center_z,
    &alpha, &beta,
    &inside_flag,
    grid);
}


void createIntersectionOfHyperboloids(
  QSSLIB_REAL *phi, int num_hyperboloids,
  QSSLIB_REAL *tangent_x, QSSLIB_REAL *tangent_y, QSSLIB_REAL *tangent_z,
  QSSLIB_REAL *center_x, QSSLIB_REAL *center_y, QSSLIB_REAL *center_z,
  QSSLIB_REAL *alpha, QSSLIB_REAL *beta,
  int *inside_flag,
  Grid *grid)
{
  int     i, j, k, l, nx, nxy;
  
  nx = (grid->grid_dims_ghostbox)[0];
  nxy = (grid->grid_dims_ghostbox)[0]*(grid->grid_dims_ghostbox)[1];
  
  for (k = 0; k < (grid->grid_dims_ghostbox)[2]; k++)
  {
    for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
    {
      for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
      {
        QSSLIB_REAL  x, y, z, max;
        int idx = i+j*nx + k*nxy;

        x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
        y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
        z = (grid->x_lo_ghostbox)[2] + (grid->dx)[2]*k;
      
        max = -FLT_MAX;
        for (l = 0; l < num_hyperboloids; l++)
        {
          QSSLIB_REAL level_set_function_value;
          QSSLIB_REAL norm_sq_x_minus_center;
          QSSLIB_REAL sq_dist_to_axis;
          QSSLIB_REAL dist_along_axis;

          norm_sq_x_minus_center = (x - center_x[l])*(x - center_x[l])
                                 + (y - center_y[l])*(y - center_y[l])
                                 + (z - center_z[l])*(z - center_z[l]);

          dist_along_axis = (x - center_x[l])*tangent_x[l]
                          + (y - center_y[l])*tangent_y[l]
                          + (z - center_z[l])*tangent_z[l];
          dist_along_axis /= sqrt( tangent_x[l]*tangent_x[l]
                                 + tangent_y[l]*tangent_y[l]
                                 + tangent_z[l]*tangent_z[l]);

          sq_dist_to_axis = 
              norm_sq_x_minus_center - dist_along_axis*dist_along_axis;

          level_set_function_value = 
            - dist_along_axis*dist_along_axis/alpha[l]/alpha[l]
            + sq_dist_to_axis/beta[l]/beta[l] - 1;

          if (inside_flag[l] >= 0) {
            level_set_function_value = -level_set_function_value;
          }

          if (level_set_function_value > max) max = level_set_function_value; 

        }
        phi[idx] = max;
      }
    }
  } /* end loop over grid */

}


void createCone(
  QSSLIB_REAL *phi,
  QSSLIB_REAL  tangent_x,
  QSSLIB_REAL  tangent_y,
  QSSLIB_REAL  tangent_z,
  QSSLIB_REAL  center_x,
  QSSLIB_REAL  center_y,
  QSSLIB_REAL  center_z,
  QSSLIB_REAL  alpha,
  QSSLIB_REAL  beta,
  int     inside_flag,
  Grid    *grid)
{
  /* use createIntersectionOfCones() */
  createIntersectionOfCones(
    phi, 1,
    &tangent_x, &tangent_y, &tangent_z,
    &center_x, &center_y, &center_z,
    &alpha, &beta,
    &inside_flag,
    grid);
}


void createIntersectionOfCones(
  QSSLIB_REAL *phi, int num_cones,
  QSSLIB_REAL *tangent_x, QSSLIB_REAL *tangent_y, QSSLIB_REAL *tangent_z,
  QSSLIB_REAL *center_x, QSSLIB_REAL *center_y, QSSLIB_REAL *center_z,
  QSSLIB_REAL *alpha, QSSLIB_REAL *beta,
  int *inside_flag,
  Grid *grid)
{
  int     i, j, k, l, nx, nxy;
  
  nx = (grid->grid_dims_ghostbox)[0];
  nxy = (grid->grid_dims_ghostbox)[0]*(grid->grid_dims_ghostbox)[1];
  
  for (k = 0; k < (grid->grid_dims_ghostbox)[2]; k++)
  {
    for (j = 0; j < (grid->grid_dims_ghostbox)[1]; j++)
    {
      for (i = 0; i < (grid->grid_dims_ghostbox)[0]; i++) 
      {
        QSSLIB_REAL  x, y, z, max;
        int idx = i+j*nx + k*nxy;

        x = (grid->x_lo_ghostbox)[0] + (grid->dx)[0]*i;
        y = (grid->x_lo_ghostbox)[1] + (grid->dx)[1]*j;
        z = (grid->x_lo_ghostbox)[2] + (grid->dx)[2]*k;
      
        max = -FLT_MAX;
        for (l = 0; l < num_cones; l++)
        {
          QSSLIB_REAL level_set_function_value;
          QSSLIB_REAL norm_sq_x_minus_center;
          QSSLIB_REAL dist_to_axis;
          QSSLIB_REAL dist_along_axis;

          norm_sq_x_minus_center = (x - center_x[l])*(x - center_x[l])
                                 + (y - center_y[l])*(y - center_y[l])
                                 + (z - center_z[l])*(z - center_z[l]);

          dist_along_axis = (x - center_x[l])*tangent_x[l]
                          + (y - center_y[l])*tangent_y[l]
                          + (z - center_z[l])*tangent_z[l];
          dist_along_axis /= sqrt( tangent_x[l]*tangent_x[l]
                                 + tangent_y[l]*tangent_y[l]
                                 + tangent_z[l]*tangent_z[l]);

          dist_to_axis = 
              sqrt(norm_sq_x_minus_center - dist_along_axis*dist_along_axis);

          level_set_function_value = 
            - dist_along_axis*dist_along_axis/alpha[l]/alpha[l]
            + dist_to_axis*dist_to_axis/beta[l]/beta[l];

          if (inside_flag[l] >= 0) {
            level_set_function_value = -level_set_function_value;
          }

          if (level_set_function_value > max) max = level_set_function_value; 

        }
        phi[idx] = max;
      }
    }
  } /* end loop over grid */

}


void createBox(
  QSSLIB_REAL *phi,
  QSSLIB_REAL corner_x, QSSLIB_REAL corner_y, QSSLIB_REAL corner_z,
  QSSLIB_REAL side_length_x, QSSLIB_REAL side_length_y, QSSLIB_REAL side_length_z,
  int inside_flag,
  Grid *grid)
{
  /* use createIntersectionOfBoxes() */
  createIntersectionOfBoxes(
   phi, 1,
   &corner_x, &corner_y, &corner_z,
   &side_length_x, &side_length_y, &side_length_z,
   &inside_flag,grid);
}     


void createIntersectionOfBoxes(
  QSSLIB_REAL *phi,
  int num_cuboids,   
  QSSLIB_REAL *corner_x, QSSLIB_REAL *corner_y, QSSLIB_REAL *corner_z,
  QSSLIB_REAL *side_length_x, QSSLIB_REAL *side_length_y, QSSLIB_REAL *side_length_z,
  int *inside_flag,
  Grid *grid)
{   
  QSSLIB_REAL  point_x[6], point_y[6], point_z[6];
  QSSLIB_REAL  normal_x[6], normal_y[6], normal_z[6];
  int     i, l, num_planes;
     
  QSSLIB_REAL  *phi1 = 
    (QSSLIB_REAL *)malloc(grid->num_gridpts*sizeof(QSSLIB_REAL)); 
        
  for(l = 0; l < num_cuboids; l++)
  {

    /* Each cuboid is intersection of 6 half spaces */
    for(i = 0; i < 3; i++)
    {
      point_x[i] = corner_x[l];
      point_y[i] = corner_y[l];
      point_z[i] = corner_z[l];
    
      if (i == 0) normal_x[i] = -1;
      else        normal_x[i] = 0;
    
      if (i == 1) normal_y[i] = -1;
      else        normal_y[i] = 0;
    
      if (i == 2) normal_z[i] = -1;
      else        normal_z[i] = 0;
    }
       
    for(i = 3; i < 6; i++)
    {
      point_x[i] = corner_x[l] + side_length_x[l];
      point_y[i] = corner_y[l] + side_length_y[l];
      point_z[i] = corner_z[l] + side_length_z[l];
    
      if (i == 3) normal_x[i] = 1;
      else        normal_x[i] = 0;
    
      if (i == 4) normal_y[i] = 1;
      else        normal_y[i] = 0;
    
      if (i == 5) normal_z[i] = 1;
      else        normal_z[i] = 0;
    }
      
    num_planes = 6;
    createIntersectionOfHalfSpaces3d(phi1,num_planes,
				     normal_x,normal_y,normal_z,
				     point_x,point_y,point_z,
				     grid);
    if(inside_flag[l] >= 1) 
    { 
	NEGATE_DATA(phi1,grid)
    }

    if(l == 0 ) COPY_DATA(phi,phi1,grid)
    else        IMPOSE_MASK(phi,phi1,phi,grid)				      
  }
    
  free(phi1);
}

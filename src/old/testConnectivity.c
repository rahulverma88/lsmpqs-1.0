/* System headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "QSSLIB_config.h"
#include "qss_util3d.h"
#include "qss_grid.h"
#include "qss_data_arrays.h"
#include "qss_options.h"
#include "qss_initialization3d.h"
#include "qss_macros.h"
#include "connectivity.h"

#define NUM_SPH 50

int  main()
{
    /* structure containing all arrays */
  QSS_DataArrays *data_arrays; 
  /* grid structure */
  Grid *grid;   
  
  TrappedPhase  *tp_nw = NULL, *tp_w = NULL;
  
  int     i;
  char    fname[256], fname_data_in[256], fname_mask_in[256], fname_grid_in[256];
  time_t t;
  
  sprintf(fname_data_in,"data_init.gz");
  sprintf(fname_mask_in,"mask.gz");
  sprintf(fname_grid_in,"grid.gz");
  
  grid = readGridFromBinaryFile(fname_grid_in);

  QSSLIB_REAL   mask_sign = -1;

  FILE    *fp_out; 
  
  fp_out = fopen("out_file","w");
  
  
  QSSLIB_REAL    center_x[NUM_SPH], center_y[NUM_SPH], center_z[NUM_SPH], radius[NUM_SPH];
  int inside_flags[NUM_SPH];
  
  data_arrays = allocateQSSDataArrays();
  allocateMemoryForQSSDataArrays(data_arrays,grid); 
  Options *options = createOptionsDefault();
  
  /* Intializes random number generator */
   srand((unsigned) time(&t));
   
  radius[0] = 0.8*(grid->x_hi[0] - grid->x_lo[0]);
  center_x[0] = (grid->x_lo[0] - radius[0]*0.7); //-0.1
  center_y[0] = 0.5*(grid->x_lo[1] + grid->x_hi[1]);
  center_z[0] = 0; //0.5*(grid->x_lo[2] + grid->x_hi[2]);
  
  //createSphere(data_arrays->phi, center_x[0], center_y[0], center_z[0], radius[0], -1, grid);

  
  for (i = 0; i < NUM_SPH; i++)
  {
        radius[i] = (100*grid->dx[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_x[i] = grid->x_lo[0] + (grid->x_hi[0] - grid->x_lo[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_y[i] = grid->x_lo[1] + (grid->x_hi[1] - grid->x_lo[1])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_z[i] = grid->x_lo[2] + (grid->x_hi[2] - grid->x_lo[2])*(QSSLIB_REAL)(rand() % 100)/100.0;
        inside_flags[i] = 1;
        
        //printf("center for sphere %d, of radius %.3f = (%.3f, %.3f, %.3f)\n",i, radius[i], center_x[i], center_y[i], center_z[i]);
        //createSphere(data_arrays->phi, center_x[i], center_y[i], center_z[i], radius[i], inside_flags[i], grid);
  }
  


  createIntersectionOfSpheres(data_arrays->phi, NUM_SPH, center_x, center_y, center_z, radius, inside_flags, grid);
  data_arrays->mask = readDataArrayQSS(grid->grid_dims_ghostbox,fname_mask_in);
  
  NEGATE_DATA(data_arrays->phi, grid);
  
  IMPOSE_MASK(data_arrays->phi, data_arrays->mask, data_arrays->phi, grid);
  
  sprintf(fname, "connectivityTestImage");
  writeDataArrayQSS(data_arrays->phi, grid, fname, GZIP);
    
  tp_nw = allocateAndInitializeTrappedPhase(options,grid,fp_out,data_arrays->phi,-1.0,
                                                              data_arrays->mask,-1.0);
  tp_w = allocateAndInitializeTrappedPhase(options,grid,fp_out,data_arrays->phi, 1.0,
                                                              data_arrays->mask,-1.0);
  unsigned char *connectivity;
  
  connectivity = (unsigned char *)malloc((grid->num_gridpts)*sizeof(unsigned char));
  
  IMPOSE_UCHAR(connectivity, data_arrays->phi, grid);
  sprintf(fname, "ucharImage");
  writeDataArrayUchar(connectivity, grid, fname, GZIP);
   
  IMPOSE_UCHAR(connectivity, data_arrays->mask, grid);
  sprintf(fname, "ucharMask3d");
  writeDataArrayUchar(connectivity, grid, fname, GZIP);
  
  //printf("checking for trapped nw blobs\n");
  //checkForTrappedBlobs(tp_nw,data_arrays,grid,mask_sign,fp_out);
  //printf("checking for trapped w blobs\n");
  //checkForTrappedBlobs(tp_w,data_arrays,grid,mask_sign,fp_out);
  //printf("trapping checking done\n");
  
  destroyTrappedPhase(tp_nw);
  destroyTrappedPhase(tp_w);
  
  free(connectivity);

  fclose(fp_out);
  
  return 0;
}     
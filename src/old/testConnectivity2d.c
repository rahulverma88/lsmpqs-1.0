/* System headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "QSSLIB_config.h"
#include "qss_grid.h"
#include "qss_data_arrays.h"
#include "qss_options.h"
#include "qss_initialization2d.h"
#include "qss_macros.h"
#include "connectivity.h"

#define NUM_SPH_MASK 100
#define NUM_SPH 50

int  main()
{
    /* structure containing all arrays */
  QSS_DataArrays *data_arrays; 
  /* grid structure */
  Grid *grid;   
  time_t t;

  /* Initializes random number generator */
  srand((unsigned) time(&t));
  
  //TrappedPhase  *tp_nw = NULL, *tp_w = NULL;
  
  int     i;
  char    fname[256], fname_data_in[256], fname_grid_in[256];
  
  QSSLIB_REAL x_lo[3] = {0, 0, 0}, x_hi[3] = {10, 10, 1};
  
  Options *options = createOptionsDefault();

  grid = createGridSetDx(2,options->dx,x_lo,x_hi,
                     (QSSLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE)options->accuracy_id);
   
  QSSLIB_REAL   mask_sign = -1;

  FILE    *fp_out; 
  
  fp_out = fopen("out_file","w");
  
  QSSLIB_REAL    center_x_mask[NUM_SPH_MASK], center_y_mask[NUM_SPH_MASK], radius_mask[NUM_SPH_MASK];
  int inside_flags_mask[NUM_SPH_MASK];
  
  QSSLIB_REAL    center_x[NUM_SPH], center_y[NUM_SPH], radius[NUM_SPH];
  int inside_flags[NUM_SPH];
  
  data_arrays = allocateQSSDataArrays();
  allocateMemoryForQSSDataArrays(data_arrays,grid); 
  
  for (i = 0; i < NUM_SPH_MASK; i++)
  {
        radius_mask[i] = (10*grid->dx[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_x_mask[i] = grid->x_lo[0] + (grid->x_hi[0] - grid->x_lo[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_y_mask[i] = grid->x_lo[1] + (grid->x_hi[1] - grid->x_lo[1])*(QSSLIB_REAL)(rand() % 100)/100.0;
        inside_flags_mask[i] = 1;
        
  } 

  createIntersectionOfCircles(data_arrays->mask, NUM_SPH_MASK, center_x_mask, center_y_mask, radius_mask, inside_flags_mask, grid);
  
  for (i = 0; i < NUM_SPH; i++)
  {
        radius[i] = (100*grid->dx[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_x[i] = grid->x_lo[0] + (grid->x_hi[0] - grid->x_lo[0])*(QSSLIB_REAL)(rand() % 100)/100.0;
        center_y[i] = grid->x_lo[1] + (grid->x_hi[1] - grid->x_lo[1])*(QSSLIB_REAL)(rand() % 100)/100.0;
        inside_flags[i] = 1;
        
  }
  


  createIntersectionOfCircles(data_arrays->phi, NUM_SPH, center_x, center_y, radius, inside_flags, grid);
  
  NEGATE_DATA(data_arrays->phi, grid);
  
  //IMPOSE_MASK(data_arrays->phi, data_arrays->mask, data_arrays->phi, grid);
  
  sprintf(fname, "connectivityTestImage2d");
  writeDataArrayQSS(data_arrays->phi, grid, fname, GZIP);
  
  sprintf(fname, "maskTestImage2d");
  writeDataArrayQSS(data_arrays->mask, grid, fname, GZIP);
    
  /*
  tp_nw = allocateAndInitializeTrappedPhase(options,grid,fp_out,data_arrays->phi,-1.0,
                                                              data_arrays->mask,-1.0);
  tp_w = allocateAndInitializeTrappedPhase(options,grid,fp_out,data_arrays->phi, 1.0,
                                                              data_arrays->mask,-1.0);
   */
                                                        
  unsigned char *connectivity;
  
  connectivity = (unsigned char *)malloc((grid->num_gridpts)*sizeof(unsigned char));
  
  IMPOSE_UCHAR(connectivity, data_arrays->phi, grid);
  sprintf(fname, "ucharImage2d");
  writeDataArrayUchar(connectivity, grid, fname, GZIP);
   
  IMPOSE_UCHAR(connectivity, data_arrays->mask, grid);
  sprintf(fname, "ucharMask2d");
  writeDataArrayUchar(connectivity, grid, fname, GZIP);
  
  //printf("checking for trapped nw blobs\n");
  //checkForTrappedBlobs(tp_nw,data_arrays,grid,mask_sign,fp_out);
  //printf("checking for trapped w blobs\n");
  //checkForTrappedBlobs(tp_w,data_arrays,grid,mask_sign,fp_out);
  //printf("trapping checking done\n");
  
  //destroyTrappedPhase(tp_nw);
  //destroyTrappedPhase(tp_w);
  
  free(connectivity);
  fclose(fp_out);
  
  return 0;
}     
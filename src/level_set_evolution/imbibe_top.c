/* Make this so that it works for both 3d and 2d. Should also act as a check for bad
   input data. At the end of the simulation, this would also write output files */
   
/* sSystem headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "QSSLIB_config.h"
#include "qss_util3d.h"
#include "qss_grid.h"
#include "drain_top.h"
#include "qss_data_arrays.h"
#include "imbibe_model3d.h"
#include "imbibe_model2d.h"

int  imbibeTop(
     Options     *options,
     char        *fname_data_in,
     char        *fname_grid_in,
     char        *fname_mask_in)                   
{ 
  /* structure containing all arrays */
  QSS_DataArrays *data_arrays; 
  /* grid structure */
  Grid *grid;   
   /* time parameters */
  time_t   time0, time1;
 
  int     n1[3], n2[3], i;
  char    fname[256];
  FILE    *fp_out; 
  
  time(&time0);

  data_arrays = allocateQSSDataArrays();


  if( (fname_data_in == NULL) || (fname_grid_in == NULL) || 
      (fname_mask_in == NULL) )
  {
     /* Allocate/initialize mask and Grid structure.
        Mask is assumed to be a level set function that is negative in the
	pore space and positive elsewhere.
        
     grid = createMaskThroatFromSpheres3d(&(data_arrays->mask),options); 
     
      Allocate/initialize (planar) initial interface position 
     normalx = 1.0; normaly = 0.0; normalz = 0.0;
     pointx  = (grid->x_lo[0]) + 20 * (grid->dx)[0];
     pointy = pointz = 0.0; 
     
     data_arrays->phi = (QSSLIB_REAL *)calloc(grid->num_gridpts,sizeof(QSSLIB_REAL));
     createPlane(data_arrays->phi,normalx,normaly,normalz,pointx,pointy,
                                                                  pointz,grid);    					
     if(options->do_mask)
     {
        IMPOSE_MASK(data_arrays->phi,data_arrays->mask,data_arrays->phi,grid)
     }	
     */				
  }
  else
  {
     /* Read input data arrays and grid */
     data_arrays->phi = readDataArrayQSS(n1,fname_data_in);
     data_arrays->mask = readDataArrayQSS(n2,fname_mask_in);
     grid = readGridFromBinaryFile(fname_grid_in);
     allocateMemoryForQSSDataArrays(data_arrays,grid); 

     /* Verify that input data files describe the same geometry */
     for(i = 0; i < 3; i++)
     {
        if( (n1[i] != (grid->grid_dims_ghostbox)[i] ) || 
	    (n2[i] != (grid->grid_dims_ghostbox)[i] ))
	{
	   printf("\nInput data dimensions n1[%d]=%d, n2[%d]=%d,",
	                                          i,n1[i],i,n2[i]);
	   printf(" (grid->grid_dims_ghostbox)[%d]=%d",
	                             i,(grid->grid_dims_ghostbox)[i]);
	   printf(" don't match\n");
	   printf("\nTerminating...");
	   
	   return 1;			  
	}
     }
  }   
  
  /* Open output file */
  fp_out = fopen(options->outfile,"a+");
  if( options->print_details)
  {
    fprintf(fp_out,"Time Start %s\n",ctime(&time0));
    fprintf(fp_out,"Constant Curvature Model Level Set Method simulation\n");
      
    printOptions(options,fp_out);
    printGrid(grid,fp_out);
  }
  
  fflush(fp_out);
  
  /* Run the curvature model, only 3d supported so far */
  if( grid->num_dims == 3 )
    imbibe_model3d(options,data_arrays,grid,fp_out);	
  else if ( grid->num_dims == 2)
    imbibe_model2d(options,data_arrays,grid,fp_out);
  					
         
  /* If desired, output initial data */
  if( options->save_data )
  {
    sprintf(fname,"%s/data_final",options->path);
    writeDataArrayQSS(data_arrays->phi,grid,fname,GZIP);    
  }  
  
  
  /* clean up memory */
  destroyQSSDataArrays(data_arrays);
  destroyGrid(grid); 
  fclose(fp_out);

  
  return 0;
}


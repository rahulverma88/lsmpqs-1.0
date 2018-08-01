/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file drain.c

     Progresive Quasistatic (PQS) Algorithm for Drainage

     Usage:
       'drain in_file data_in grid_in mask_in' - 
            data_in, grid_in and mask_in are input binary files for
            fluid level set function, grid structure and masking level set
            function. If the files are zipped, the names should have .gz 
            extension.
     Output:
      data_init.gz   - data array storing initial level set function 
                      (NW phase == data_init < 0)
      data_stepID.gz - data array storing the NW phase level set function
                      for each simulation step ID
		      For imbibition, ID starts from 1.
      grid.gz       - binary data array storing grid information
      mask.gz       - masking level set function (solid phase == mask > 0)		      
      out_file      - options, Grid data and any other text output
      vol_frac.gz   - 1D data array storing volume fraction occupied by NW phase
                      for each step             
*/

/* System headers */
#include <stdio.h>

/* PQS package headers */
#include "qss_data_arrays.h"
#include "qss_options.h"
#include "drain_top.h"


/* Main driver for constant curvature level set method model */

int main(int argc, char **argv) 
{
   /* input filename storage */
   char  *in_fname, fname[256];
   char  *fname_data_in = (char *)NULL;
   char  *fname_grid_in = (char *)NULL;
   char  *fname_mask_in = (char *)NULL;
   
   int return_status;
   
   /* Structure containing input options and parameters.
      See lsm_options.h for details on the structure elements.
   */
   Options *options; 
   
   /* Initialize the problem */   
   if( argc == 1 )
   {   /* input file not provided, set all options to default */
       options =  createOptionsDefault();
   }
   else if (argc == 2)
   {   /* set options according to input file */
       in_fname = argv[1];
       sprintf(fname,"%s",argv[1]);
       options = createOptionsFromInputFile(fname);
   }
   else if( argc >= 4 )
   { /* read data from provided input files */
     in_fname = argv[1];
     sprintf(fname,"%s",argv[1]);
     options = createOptionsFromInputFile(fname);
//    printf("b = %f\n",options->b);
     
     fname_data_in = argv[2];
     fname_grid_in = argv[3];
     fname_mask_in = argv[4];
   }
   else
   {
     printf("\nRunning options:");
     printf("\n\t./drain");
     printf("\n\t./drain input_file");
     printf("\n\t./drain input_file data_init grid mask");
     printf("\n"); 
   }
    
   return_status = drainTop(options,fname_data_in,fname_grid_in,
                                                              fname_mask_in);   
   /* clean up memory */
   destroyOptions(options);
   return return_status;
}

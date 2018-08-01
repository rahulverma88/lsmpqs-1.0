/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file convert2Raw.c

    Converts all level set data steps between init_step and final_step to raw file format.
    Essentially, it strips the headers from the level set files, so they can be read by Paraview.
    Additionally, it takes as input the phase of interest (wetting/non-wetting: w/n).

    It assumes "grid.gz" and "mask.gz" are present in the current directory, and refer to
    the grid binary file and mask level set file, respectively.
    
    The main difference between makeRaw and convert2Raw is that makeRaw produces ubc files (0/1),
    whereas convert2Raw produces floating point output - essentially just stripping headers
    from the level set output files.
    
    Usage:
        'convert2Raw init_step final_step phase'
    
    Input:
        init_step: initial data step
        final_step: final data step
        phase: can take values 'w' and 'n', referring to the wetting and non-wetting phases, 
            respectively. phi < 0 would refer to the non-wetting phase, and phi > 0 the wetting.
    
    Output:
        data_step_<step>_<phase>.raw  : Raw file which contains the converted level set to raw file
            including the ghost cells.
        mask.raw: Raw file containing converted mask to ubc        		                
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_data_arrays.h"
#include "omp.h"
#include "qss_macros.h"
#include "qss_reinitialization3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_general_util.h"
#include "qss_util3d.h"
#include "qss_util2d.h"


int main(int argc, char **argv)
{
    /* input filename storage */
   char  fname[256];
   int     n1[3], n2[3], i;
   
   /* structure containing all arrays */
    QSS_DataArrays *p; 
    /* grid structure */
    Grid *g; 
  
    int init_step = atoi(argv[1]);
    int final_step = atoi(argv[2]);
    char phase = argv[3][0];

    p = allocateQSSDataArrays();
    g = readGridFromBinaryFile("grid.gz");
    allocateMemoryForQSSDataArrays(p,g); 
    p->mask = readDataArrayQSS(n2,"mask.gz");

    NEGATE_DATA(p->mask, g);
    sprintf(fname,"mask.raw");
    writeDataArrayRaw(p->mask, g, fname, 0);

    
    for (i = init_step; i <= final_step; i++) {
    	sprintf(fname,"data_step_%d.gz", i);
	p->phi = readDataArrayQSS(n2,fname);

	if(phase == 'n') {
		sprintf(fname,"data_step_%d_nw.raw", i);
    		writeDataArrayRaw(p->phi, g, fname, 0);
    	
	} else if (phase == 'w') {
		COPY_DATA(p->scratch1, p->phi, g);
		NEGATE_DATA(p->scratch1, g);
		IMPOSE_MASK(p->scratch1,p->mask,p->scratch1,g);
		
		sprintf(fname,"data_step_%d_w.raw", i);
    		writeDataArrayRaw(p->phi, g, fname, 0);
	}
    }

    destroyQSSDataArrays(p);

}
    


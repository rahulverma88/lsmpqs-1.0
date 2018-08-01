/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file makeRaw.c

    Converts all level set data steps between init_step and final_step to ubc raw file format.
    Additionally, it takes as input the phase of interest (wetting/non-wetting: w/n).

    It assumes "grid.gz" and "mask.gz" are present in the current directory, and refer to
    the grid binary file and mask level set file, respectively.
    
    The main difference between makeRaw and convert2Raw is that makeRaw produces ubc files (0/1),
    whereas convert2Raw produces floating point output - essentially just stripping headers
    from the level set output files.
    
    Usage:
        'makeRaw init_step final_step phase'
    
    Input:
        init_step: initial data step
        final_step: final data step
        phase: can take values 'w' and 'n', referring to the wetting and non-wetting phases, 
            respectively. phi < 0 would refer to the non-wetting phase, and phi > 0 the wetting.
    
    Output:
        data_step_<step>_<phase>_uchar.raw  : Raw file which contains the converted level set to ubc file
        data_step_<step>_<phase>_uchar_fb.raw  : Raw file which contains the converted level set to ubc file,
            but only the "fill box", that is it removes the ghost cells.
        mask_uchar.raw: Raw file containing converted mask to ubc
        mask_uchar_fb.raw: Raw file containing mask to ubc, but only the "fill box".
        		                
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

void copyDataFB(unsigned char *data_fb, unsigned char *data_full, Grid *g)
{
    int i, j, k, idx, idx_gb;
    int num_gc = 3;
    int nx_vox = g->ihi_fb - g->ilo_fb + 1;
    int ny_vox = g->jhi_fb - g->jlo_fb + 1;
    int nx_gb = g->ihi_gb - g->ilo_gb + 1;
    int nxy_gb = nx_gb *( g->jhi_gb - g->jlo_gb + 1);
    
    for( k = g->klo_fb; k <= g->khi_fb; k++)
    {
        for( j = g->jlo_fb; j <= g->jhi_fb; j++)
        {
            for( i = g->ilo_fb; i <= g->ihi_fb; i++)
            {
                idx = (i - num_gc) + (j - num_gc)*nx_vox + (k - num_gc)*nx_vox*ny_vox;
                idx_gb = i + j*nx_gb + k*nxy_gb;
                data_fb[idx] = data_full[idx_gb];
            }
        }
    }
            
}
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
    int num_fb_pts = g->grid_dims[0]*g->grid_dims[1]*g->grid_dims[2];
    
    unsigned char *bin_data = (unsigned char *)malloc(num_fb_pts*UCSZ);
    
    IMPOSE_UCHAR(p->phi_bin, p->mask, g);
    copyDataFB(bin_data, p->phi_bin, g);
    
    sprintf(fname,"mask.raw");
    writeDataArrayUchar(p->phi_bin, g, fname, 0);

    sprintf(fname, "mask_fb.raw");
    writeDataArrayUcharFB(bin_data, g, fname, 0);
    
    for (i = init_step; i <= final_step; i++) {
    	sprintf(fname,"data_step_%d.gz", i);
	p->phi = readDataArrayQSS(n2,fname);

	if(phase == 'n') {
		IMPOSE_UCHAR(p->phi_bin, p->phi, g);
        copyDataFB(bin_data, p->phi_bin, g);
    
		sprintf(fname,"data_step_%d_nw_uchar.raw", i);
    	writeDataArrayUchar(p->phi_bin, g, fname, 0);
    	
    	sprintf(fname,"data_step_%d_nw_uchar_fb.raw", i);
    	writeDataArrayUchar(bin_data, g, fname, 0);
    	
    	
	} else if (phase == 'w') {
		COPY_DATA(p->scratch1, p->phi, g);
		NEGATE_DATA(p->scratch1, g);
		IMPOSE_MASK(p->scratch1,p->mask,p->scratch1,g);

		IMPOSE_UCHAR(p->phi_bin, p->scratch1, g);
		copyDataFB(bin_data, p->phi_bin, g);
		
		sprintf(fname,"data_step_%d_w_uchar_fb.raw", i);
    	writeDataArrayUchar(bin_data, g, fname, 0);

		sprintf(fname,"data_step_%d_w_uchar.raw", i);
    	writeDataArrayUchar(p->phi_bin, g, fname, 0);
	}
    }

    destroyQSSDataArrays(p);

}
    


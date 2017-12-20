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
    


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
    int second_step = atoi(argv[2]);

    p = allocateQSSDataArrays();
    g = readGridFromBinaryFile("grid.gz");
    allocateMemoryForQSSDataArrays(p,g); 
   


    sprintf(fname,"data_step_%d.gz", init_step);
    p->mask = readDataArrayQSS(n2,fname);
    NEGATE_DATA(p->mask, g);


    sprintf(fname,"data_step_%d.gz", second_step);
    p->phi = readDataArrayQSS(n2,fname);

    for (i = 0; i < g->num_gridpts; i++)
	p->phi[i] -= 0.5*g->dx[0];

    IMPOSE_MASK(p->phi,p->mask,p->phi,g);
    	
    sprintf(fname,"diff_steps_%d_%d.raw", init_step, second_step);
    writeDataArrayRaw(p->phi, g, fname, 0);

    destroyQSSDataArrays(p);

}
    

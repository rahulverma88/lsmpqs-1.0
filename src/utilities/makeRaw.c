/* Probably obsolete */
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
  
    p = allocateQSSDataArrays();
    g = readGridFromBinaryFile("grid.gz");
    allocateMemoryForQSSDataArrays(p,g); 
    p->mask = readDataArrayQSS(n2,"mask.gz");
    p->phi = readDataArrayQSS(n2,"data_step_29.gz");
    
    sprintf(fname,"maskRaw.raw");
    writeDataArrayRaw(p->mask, g, fname, 0);
    
    sprintf(fname,"dataStep29Raw.raw");
    writeDataArrayRaw(p->phi, g, fname, 0);
    
}
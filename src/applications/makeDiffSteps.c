#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "QSSLIB_config.h"
#include "qss_options.h"
#include "qss_data_arrays.h"
#include "omp.h"
#include "qss_macros.h"
#include "connectivity.h"
#include "qss_reinitialization3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_general_util.h"
#include "qss_util3d.h"
#include "qss_util2d.h"

/*
    Create lenses trapped between two level set steps at different curvatures.
    Needs two arguments: argument 1 is the "mask", and argument 2 is the level set
    which is masked. So if the first level set occupies 70 % of the pore space, and the second
    occupies 50 %, then the trapped lens should occupy 20% of the pore space.
*/
int main(int argc, char **argv)
{
    /* input filename storage */
   char  fname[256];
   int     n1[3], n2[3];
   int idx_gb, nx_gb, nxy_gb, i, j, k, ind;
   int main_comp_val, inlet_ind, outlet_ind;
   
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

    QSSLIB_REAL eps = 0;
    IMPOSE_INT_EPS(p->phi_bin, p->phi, g, eps);
    nx_gb = g->grid_dims_ghostbox[0];
    nxy_gb = g->grid_dims_ghostbox[1]*nx_gb;
    
    /* add a layer of 1's in slice ilo_fb - 1, to simulate connectivity to reservoir */
    /* That way, any nw blob connected to inlet should become part of  nw-connected phase */
    if (g->num_dims == 3) {
        i = g->ilo_fb - 1;
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                idx_gb = i + j*nx_gb + k*nxy_gb;
                p->phi_bin[idx_gb] = 1; /* make fb-1 slice 1s to simulate inlet */
            }
        }
        
        i = g->ihi_fb + 1;
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                idx_gb = i + j*nx_gb + k*nxy_gb;
                p->phi_bin[idx_gb] = 1; /* make fb+1 slice 1s to simulate outlet */
            }
        }
        
        findConnectivity3d(p, g); 
        /* The main component is the one connected to the inlet and outlet - the one which percolates through */
        i = g->ilo_fb - 1;
        j = 0.5*(g->jlo_fb + g->jhi_fb);
        k = 0.5*(g->klo_fb + g->khi_fb); 
        inlet_ind = i + j*nx_gb + k*nxy_gb;
        
        i = g->ihi_fb + 1;
        outlet_ind = i + j*nx_gb + k*nxy_gb;
        
        if( p->connectivity[inlet_ind] != p->connectivity[outlet_ind])
        {
            printf("No percolation!\n");
        } else
        {
            main_comp_val = p->connectivity[inlet_ind];
            for (ind = 0; ind < g->num_gridpts; ind++)
            {
                if ( (p->connectivity[ind] != main_comp_val) & (p->phi[ind] <= 0))
                    p->phi[ind] *= -1;
            }
            
            sprintf(fname,"diff_steps_%d_%d_connected_only.raw", init_step, second_step);
            writeDataArrayRaw(p->phi, g, fname, 0);
        }    
         
    }
    
    destroyQSSDataArrays(p);

}
    

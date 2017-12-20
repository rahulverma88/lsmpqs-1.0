
/* System headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "QSSLIB_config.h"
#include "qss_initialization3d.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_data_arrays.h"
#include "qss_macros.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_reinitialization3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_general_util.h"

int main(int argc, char *argv[])
{
    
    LSMLIB_REAL  x_lo[3], x_hi[3];
    int  n[3], n_local[3], i, j, k, idx, idx1, dim, idx_gb;
    char    fname[256];	

    sprintf(fname,"%s",argv[1]);

    int nxyz;
    QSS_DataArrays *p; 
    p = allocateQSSDataArrays();
    Options *options;
    options =  createOptionsDefault();
    int  nx, nx_pad, nxy, nxy_pad, nx_gb, nxy_gb, nxyz_gb;

    int nx_vox = atoi(argv[2]), ny_vox = atoi(argv[3]), nz_vox = atoi(argv[4]);
    
    int     num_gc = 3;
    
    /* n[3] refers to dimensions of fill box */
    n[0] = nx_vox; n[1] = ny_vox; n[2] = nz_vox;
    
    nx = n[0]; nxy = nx*n[1];
    nxyz = nxy*n[2];
    
    nx_gb = n[0] + 2*num_gc; nxy_gb = nx_gb*(n[1] + 2*num_gc);
    nxyz_gb = nxy_gb*(n[2] + 2*num_gc);
        
    if( n[2] > 1) dim = 3;
    else          dim = 2;

    unsigned char *data;
    FILE *fp = fopen(fname,"r");
    
    /* data gets the exact geometry file, which should just be the fill box */
    data = (unsigned char *)malloc((nx_vox*ny_vox*nz_vox)*sizeof(unsigned char));

    fread(data, sizeof(unsigned char), (nx_vox*ny_vox*nz_vox), fp); 

    fclose(fp);
    
    Grid *g;

    /* x_lo refers to fill box. For ghostbox, it is x_lo_ghostbox */
    x_lo[0] = 0 + num_gc*options->dx; x_hi[0] = x_lo[0] + n[0]*options->dx;
    x_lo[1] = 0 + num_gc*options->dx; x_hi[1] = x_lo[1] + n[1]*options->dx;
    x_lo[2] = 0 + num_gc*options->dx; x_hi[2] = x_lo[2] + n[2]*options->dx;
            
    g = createGridSetGridDims(dim,n,x_lo,x_hi,(QSSLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE)options->accuracy_id);
    sprintf(fname,"grid");
    writeGridToBinaryFile(g,fname,GZIP);
        
    allocateMemoryForQSSDataArrays(p,g);       	
        
    for( k = g->klo_fb; k <= g->khi_fb; k++)
    {
        for( j = g->jlo_fb; j <= g->jhi_fb; j++)
        {
            for( i = g->ilo_fb; i <= g->ihi_fb; i++)
            {
                idx = (i - num_gc) + (j - num_gc)*nx_vox + (k - num_gc)*nx_vox*ny_vox;
                idx_gb = i + j*nx_gb + k*nxy_gb;
                if( data[idx] == 1)  p->mask[idx_gb] =  options->dx;  //grain space
                else                  p->mask[idx_gb] = -options->dx; //pore space
                
            }
        }
    }
    
    
    qss_reinitialize_mask(p, g, options);

    sprintf(fname,"mask");
    writeDataArrayQSS(p->mask,g,fname,GZIP);
    
    destroyQSSDataArrays(p);
    destroyOptions(options);
    destroyGrid(g);     
    free(data);

    return 0;
}


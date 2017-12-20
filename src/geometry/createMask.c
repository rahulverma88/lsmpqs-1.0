
/* System headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "QSSLIB_config.h"

#include "qss_macros.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "qss_data_arrays.h"
#include "qss_util3d.h"
#include "qss_reinitialization3d.h"
#include "qss_tvd_runge_kutta3d.h"
#include "qss_spatial_derivatives3d.h"
#include "qss_initialization3d.h"
#include "qss_general_util.h"

#define DSZ  sizeof(QSSLIB_REAL)


int main(int argc, char *argv[])
{
    
    QSSLIB_REAL  x_lo[3], x_hi[3];
    int  n[3], n_local[3], i, j, k, idx, idx1, dim;
    char    fname[256];	
    char  *fname_ubc = (char *)NULL;

    sprintf(fname,"%s",argv[1]);

    size_t nxyz;
    QSS_DataArrays *p; 
    p = allocateQSSDataArrays();
    Options *options;
    options =  createOptionsDefault();
    int  nx, nx_pad, nxy, nxy_pad;

    int nx_vox = atoi(argv[2]), ny_vox = atoi(argv[3]), nz_vox = atoi(argv[4]);
    
    int     num_gc = 3;
    n[0] = nx_vox - num_gc*2; n[1] = ny_vox - num_gc*2; n[2] = nz_vox - num_gc*2;

    nx = n[0] + 2*num_gc; nxy = nx*(n[1]+2*num_gc);
    nxyz = nxy*(n[2]+2*num_gc);
    
    if( n[2] > 1) dim = 3;
    else          dim = 2;

    unsigned char *data;
    FILE *fp = fopen(fname,"r");
    data = (unsigned char *)malloc(nxyz*sizeof(unsigned char));

    fread(data, sizeof(unsigned char), nxyz, fp); 

    fclose(fp);
    
    Grid *g;

    
    x_lo[0] = 0 + num_gc*options->dx; x_hi[0] = x_lo[0] + n[0]*options->dx;
    x_lo[1] = 0 + num_gc*options->dx; x_hi[1] = x_lo[1] + n[1]*options->dx;
    x_lo[2] = 0 + num_gc*options->dx; x_hi[2] = x_lo[2] + n[2]*options->dx;
        
//Accuracy_settings_menu[options->accuracy_id].num_ghostcells;
    
    g = createGridSetGridDims(dim,n,x_lo,x_hi,(QSSLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE)options->accuracy_id);
    sprintf(fname,"grid");
    writeGridToBinaryFile(g,fname,GZIP);
        
    allocateMemoryForQSSDataArrays(p,g);     
    
    QSSLIB_REAL    center_x, center_y, center_z, radius, normalx, normaly, normalz;
    QSSLIB_REAL    pointx, pointy, pointz, x_hi_first;
    
    center_y = 0.5*(x_lo[1] + x_hi[1]);
    
    radius = 0.8*(x_hi[0] - x_lo[0]);
    center_x = (x_lo[0] - radius*0.7); //-0.1

    x_hi_first = n[2]*options->dx;
    center_z = 0.5*(x_lo[2] + x_hi[2]);
//    printf("center_x = %lf, center_y = %lf, center_z = %lf\n",center_x,center_y,center_z);

    createSphere(p->phi,center_x,center_y,center_z,radius,-1,g);    	
        
    for( k = g->klo_gb; k <= g->khi_gb; k++)
    {
        for( j = g->jlo_gb; j <= g->jhi_gb; j++)
        {
            for( i = g->ilo_gb; i <= g->ihi_gb; i++)
            {
                //idx = i + j*nx + k*nxy;
                idx1 = i + j*nx + k*nxy;
                 	//printf("%d\n",data[idx1]);
                if( data[idx1] == 1)  p->mask[idx1] =  options->dx;  //grain space
                else                  p->mask[idx1] = -options->dx; //pore space
                
            }
        }
    }

    //qss_reinitialize_mask1(p, g, options);
    IMPOSE_MASK(p->phi,p->mask,p->phi,g)  

    sprintf(fname,"mask");
    writeDataArrayQSS(p->mask,g,fname,GZIP);

    sprintf(fname,"data_init");
    writeDataArrayQSS(p->phi,g,fname,GZIP);
    
    destroyQSSDataArrays(p);
    destroyOptions(options);
    destroyGrid(g);     
    free(data);

    return 0;
}

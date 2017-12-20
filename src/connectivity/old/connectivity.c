#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "connectivity.h"
#include "qss_macros.h"

void doUnion(int a, int b, int *component)
{
    // get the root component of a and b, and set the one's parent to the other
    while (component[a] != a)
        a = component[a];
    while (component[b] != b)
        b = component[b];
    component[b] = a;
}

void unionCoords2d(int x, int y, int x2, int y2, int *component, int *input, int h, int w)
{
    int ind1 = x*h + y;
    int ind2 = x2*h + y2;
    if (y2 < h && x2 < w && input[ind1] && input[ind2] && y2 >= 0 && x2 >= 0)
        doUnion(ind1, ind2, component);
}

void findConnectivity2d(QSS_DataArrays *p, Grid *g)
{
    int i, j, x, y;
    int w, h, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
               
    int *component = (int *)malloc((w*h)*sizeof(int));

    for (i = 0; i < w*h; i++)
        component[i] = i;
        
    for (x = 0; x < w; x++)
    for (y = 0; y < h; y++)
    {
        unionCoords2d(x, y, x+1, y, component, p->phi_bin, h, w);
        unionCoords2d(x, y, x, y+1, component, p->phi_bin, h, w);
        unionCoords2d(x, y, x+1, y+1, component, p->phi_bin, h, w);
        unionCoords2d(x, y, x-1, y+1, component, p->phi_bin, h, w);
        unionCoords2d(x, y, x+1, y-1, component, p->phi_bin, h, w);
        unionCoords2d(x, y, x-1, y-1, component, p->phi_bin, h, w);
    }

    /* update connectivity array */
    for (x = 0; x < w; x++)
    {
        for (y = 0; y < h; y++)
        {
            c = x*h + y;
            if (p->phi_bin[c] == 0)
            {
                p->connectivity[c] = 0;
                continue;
            }
            while (component[c] != c) c = component[c];
            
            c1 = x*h + y;
            p->connectivity[c1] = component[c];

        }
    }
    
    free(component);
}

void unionCoords3d(int x, int y, int z, int x2, int y2, int z2, int *component, \
        int *input, int h, int w, int d)
{
    int ind1, ind2, nxy;
    
    nxy = h*w;
    ind1 = y + x*h + z*nxy;
    ind2 = y2 + x2*h + z2*nxy;
    
    if (y2 < h && x2 < w && z2 < d && input[ind1] && input[ind2] \
            && y2 >= 0 && x2 >= 0 && z2 >= 0)
        doUnion(ind1, ind2, component);
}

void findConnectivity3d(QSS_DataArrays *p, Grid *g)
{
    int i, j, k, x, y, z;
    int w, h, d, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
    d = g->grid_dims_ghostbox[2];

    int *component = (int *)malloc((w*h*d)*sizeof(int));
    
    for (i = 0; i < w*h*d; i++)
        component[i] = i;
    
    for (z = 0; z < d; z++){
        for (x = 0; x < w; x++){
            for (y = 0; y < h; y++){
                /* 6-connectivity */
                unionCoords3d(x, y, z, x+1, y, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y, z+1, component, p->phi_bin, h, w, d);
                
                unionCoords3d(x, y, z, x-1, y, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y, z-1, component, p->phi_bin, h, w, d);
                
                /* 18-connectivity */
                unionCoords3d(x, y, z, x-1, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z, component, p->phi_bin, h, w, d);

                unionCoords3d(x, y, z, x-1, y, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y-1, z+1, component, p->phi_bin, h, w, d);
                        
                unionCoords3d(x, y, z, x+1, y, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x, y+1, z-1, component, p->phi_bin, h, w, d);
                
                /* 26-connectivity */
                unionCoords3d(x, y, z, x-1, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z-1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z-1, component, p->phi_bin, h, w, d);

                unionCoords3d(x, y, z, x-1, y-1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y-1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x+1, y+1, z+1, component, p->phi_bin, h, w, d);
                unionCoords3d(x, y, z, x-1, y+1, z+1, component, p->phi_bin, h, w, d);
            }
        }
    }
    
    for (z = 0; z < d; z++){
        for (x = 0; x < w; x++){
            for (y = 0; y < h; y++){
                c = y + x*h + z*w*h;
                if (p->phi_bin[c] == 0)
                {
                    p->connectivity[c] = 0;
                    continue;
                }
                while (component[c] != c) c = component[c];
            
                c1 = y + x*h + z*w*h;
                p->connectivity[c1] = component[c];
            }
        }
    }
}

void getMainInd(Options *o, QSS_DataArrays *p, Grid *g)
{
    int i;
    for (i = 0; i < g->num_gridpts; i++)
        if (p->phi[i] < 0) {
            o->phi_nw_ind = i;
            break;
        }
        
        
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);
    
    for (i = 0; i < g->num_gridpts; i++)
        if (p->scratch1[i] < 0) {
            o->phi_w_ind = i;
            break;
        }

}
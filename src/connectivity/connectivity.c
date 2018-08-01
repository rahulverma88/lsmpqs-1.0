/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file connectivity.c

    Contains base routines to find connected components given a binary 0/1 2D/3D array,
    using union-find algorithms.
    
    The code was adapted from the following forum:
    https://stackoverflow.com/questions/44467899/connected-component-labeling-with-diagonal-connections-using-union-find
    
    Functions/Structures:
        (Applicable to both 2D and 3D):
        union-find struct: contains size of the 2D/3D array, and a pointer to a parent.
        make_sets: initializes the union-find structure
        find: returns index of parent
        do_union: combines sub-trees, making sure one branch does not become too long
        getMainInd_new: gets the index of the non-wetting/wetting inlet/outlet.
        getMainInd: is obsolete. There for historical reasons.
        
        (2D):
        unionCoords2d: creates union-find tree structure for 2D array, taking care of the boundaries.
        findConnectivity2d: finds connected components for 2D array, using unionCoords2d and find.
        
        (3D):
        unionCoords3d: creates union-find tree structure for 3D array, taking care of the boundaries.
        findConnectivity3d: finds connected components for 3D array, using unionCoords3d and find.
    
*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "connectivity.h"
#include "qss_macros.h"

typedef struct _union_find{
    int* parent;
    int size;
} union_find;

union_find make_sets(int size) {
    union_find result;
    result.parent = malloc(sizeof(int) * size);
    result.size = size;
    int i;
    for (i = 0; i < size; ++i) {
        result.parent[i] = size;
    }

    return result;
}

int find(union_find uf, int i) {
    if (uf.parent[i] < uf.size)
        return uf.parent[i] = find(uf, uf.parent[i]);
    return i;
}

void do_union(union_find uf, int i, int j) {
    int pi = find(uf, i);
    int pj = find(uf, j);
    if (pi == pj) {
        return;
    }
    if (pi < pj) {
        // link the smaller group to the larger one
        uf.parent[pi] = pj;
    } else if (pi > pj) {
        // link the smaller group to the larger one
        uf.parent[pj] = pi;
    } else {
        // equal rank: link arbitrarily and increase rank
        uf.parent[pj] = pi;
        ++uf.parent[pi];  
    }
}

void unionCoords2d(int x, int y, int x2, int y2, union_find component, unsigned char *input, int h, int w)
{
    int ind1 = x*h + y;
    int ind2 = x2*h + y2;
    if (y2 < h && x2 < w && input[ind1] && input[ind2] && y2 >= 0 && x2 >= 0)
        do_union(component, ind1, ind2);
}

void findConnectivity2d(QSS_DataArrays *p, Grid *g)
{
    int i, j, x, y;
    int w, h, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
               
    union_find component = make_sets(w*h);
        
    for (x = 0; x < w; x++)
        for (y = 0; y < h; y++){
            unionCoords2d(x, y, x+1, y, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x, y+1, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x-1, y, component, p->phi_bin, h, w);
            unionCoords2d(x, y, x, y-1, component, p->phi_bin, h, w);
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
            p->connectivity[c] = find(component, c);

        }
    }
    
    free(component.parent);
}

void unionCoords3d(int x, int y, int z, int x2, int y2, int z2, union_find component, \
        unsigned char *input, int h, int w, int d)
{
    int ind1, ind2, nxy;
    
    nxy = h*w;
    ind1 = y + x*h + z*nxy;
    ind2 = y2 + x2*h + z2*nxy;
    
    if (y2 < h && x2 < w && z2 < d && input[ind1] && input[ind2] \
            && y2 >= 0 && x2 >= 0 && z2 >= 0)
        do_union(component, ind1, ind2);
}

void findConnectivity3d(QSS_DataArrays *p, Grid *g)
{
    int i, j, k, x, y, z;
    int w, h, d, c, c1;
    
    w = g->grid_dims_ghostbox[1];
    h = g->grid_dims_ghostbox[0];
    d = g->grid_dims_ghostbox[2];

    union_find component = make_sets(w*h*d);
    
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
            
                p->connectivity[c] = find(component,c);
            }
        }
    }
    
    free(component.parent);
}

void getMainInd_new(Options *o, QSS_DataArrays *p, Grid *g)
{
    int i, j, k, idx, nw_x_ind, nw_y_ind, nw_z_ind, w_x_ind, w_y_ind, w_z_ind, ind;
    
    if (g->num_dims == 2) {
    
        if (o->center_inlet) {
            /* Juanes2d case: inlet is in the center of domain, 
              wetting outlet is on all four sides */
    
            nw_x_ind = 0.5*(g->ilo_fb + g->ihi_fb);
            nw_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
    
            o->phi_nw_ind = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]);
    
            w_x_ind = g->ilo_fb - 1;
            w_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            o->phi_w_ind = w_x_ind + w_y_ind * (g->grid_dims_ghostbox[0]);
            
        } else {
            nw_x_ind = g->ilo_fb - 1;
            nw_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
    
            o->phi_nw_ind = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]);
    
            w_x_ind = g->ilo_fb + 1;
            w_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            o->phi_w_ind = w_x_ind + w_y_ind * (g->grid_dims_ghostbox[0]);
        
        
        }
    } else if (g->num_dims == 3) {
    
        if (o->center_inlet) {
            /* inlet (or non-wetting index point) is at center of domain */
            nw_x_ind = 0.5*(g->ilo_fb + g->ihi_fb);
            nw_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            nw_z_ind = 0.5*(g->klo_fb + g->khi_fb);
    
            ind = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]) 
                + nw_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
    
            o->phi_nw_ind = ind;
    
            /* outlet (or wetting index point) is in last + 1 x-fill box boundary */
            w_x_ind = g->ihi_fb + 1;
            w_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            w_z_ind = 0.5*(g->klo_fb + g->khi_fb);
    
            ind = w_x_ind + w_y_ind * (g->grid_dims_ghostbox[0]) 
                + w_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
        
            o->phi_w_ind = ind;
                
        } else {
            /* inlet (or non-wetting index point) is in initial - 1 x- fill box boundary */
            nw_x_ind = g->ilo_fb - 1;
            nw_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            nw_z_ind = 0.5*(g->klo_fb + g->khi_fb);
    
            ind = nw_x_ind + nw_y_ind * (g->grid_dims_ghostbox[0]) 
                + nw_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
    
            o->phi_nw_ind = ind;
    
            /* outlet (or wetting index point) is in last + 1 x-fill box boundary */
            w_x_ind = g->ihi_fb + 1;
            w_y_ind = 0.5*(g->jlo_fb + g->jhi_fb);
            w_z_ind = 0.5*(g->klo_fb + g->khi_fb);
    
            ind = w_x_ind + w_y_ind * (g->grid_dims_ghostbox[0]) 
                + w_z_ind * (g->grid_dims_ghostbox[0]) * (g->grid_dims_ghostbox[1]);
        
            o->phi_w_ind = ind;
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
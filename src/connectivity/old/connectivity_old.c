#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "connectivity.h"

void doUnion(int a, int b, int *component)
{
    // get the root component of a and b, and set the one's parent to the other
    while (component[a] != a)
        a = component[a];
    while (component[b] != b)
        b = component[b];
    component[b] = a;
}

void unionCoords2d(int i, int j, int i2, int j2, QSS_DataArrays *p, Grid *g, int *component)
{
    int ind1, ind2;
    ind1 = i + j*(g->grid_dims_ghostbox)[0];
    ind2 = i2 + j2*(g->grid_dims_ghostbox)[0];
    
    if (i2 < (g->grid_dims_ghostbox)[0] && j2 < (g->grid_dims_ghostbox)[1] && \
    p->phi_bin[ind1] && p->phi_bin[ind2]){// && i2 > 0 && j2 > 0) {
        //printf("(i,j) = (%d,%d); (i2,j2) = (%d, %d); p->phi_bin[ind1] = %d; component[ind1] = %d\n",i,j,i2,j2,p->phi_bin[ind1],component[ind1]);
        doUnion(ind1, ind2, component);
    }
}

void findConnectivity2d(QSS_DataArrays *p, Grid *g)
{
    int i, j, ind;
    int *component = (int *)malloc((g->num_gridpts)*sizeof(int));
    
    for (i = 0; i < g->num_gridpts; i++)
        component[i] = i;
    
    
    for (j = 0; j < (g->grid_dims_ghostbox)[1]; j++){
        for (i = 0; i < (g->grid_dims_ghostbox)[0]; i++){
            unionCoords2d(i, j, i+1, j, p, g, component);
            unionCoords2d(i, j, i, j+1, p, g, component);
            unionCoords2d(i, j, i+1, j+1, p, g, component);
            unionCoords2d(i, j, i-1, j+1, p, g, component);
            unionCoords2d(i, j, i+1, j-1, p, g, component);
            unionCoords2d(i, j, i-1, j-1, p, g, component);
        }
    }
    
    for (j = 0; j < (g->grid_dims_ghostbox)[1]; j++){
        for (i = 0; i < (g->grid_dims_ghostbox)[0]; i++){
            ind = i + j*(g->grid_dims_ghostbox)[0];
            if (p->phi_bin[ind] == 0)
            {
                p->connectivity[ind] = 0;
                continue;
            }
            while (component[ind] != ind) ind = component[ind];
            p->connectivity[ind] = component[ind];
        }
    }

}


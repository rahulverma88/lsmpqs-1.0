#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "QSSLIB_config.h"
#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "qss_macros.h"
#include "connectivity.h"
#include "qss_general_util.h"
#include <omp.h>

void trapComponents(QSS_DataArrays *p, Grid *g, Options *o, QSSLIB_REAL val){
    char fname[256];
    int main_comp_val, main_comp_ind;
    int idx_gb, nx_gb, nxy_gb, i, j, k;
    
    nx_gb = g->grid_dims_ghostbox[0];
    nxy_gb = g->grid_dims_ghostbox[1]*nx_gb;
    
    unsigned char *bin_vel_w = (unsigned char *)malloc((g->num_gridpts)*sizeof(unsigned char));
    unsigned char *bin_vel_nw = (unsigned char *)malloc((g->num_gridpts)*sizeof(unsigned char));
    unsigned char *bin_velocity = (unsigned char *)malloc((g->num_gridpts)*sizeof(unsigned char));

    QSSLIB_REAL eps = 0;//4*g->dx[0];
    
    /* Find connectivity of non-wetting phase first */
    IMPOSE_INT_EPS(p->phi_bin, p->phi, g, eps);

    /* add a layer of 1's in slice ilo_fb - 1, to simulate connectivity to reservoir */
    /* That way, any nw blob connected to inlet should become part of  nw-connected phase */
    if (g->num_dims == 3) {
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                 i = g->ilo_fb - 1;
                 idx_gb = i + j*nx_gb + k*nxy_gb;
                 p->phi_bin[idx_gb] = 1; // make slice 1 
            }
        } 
    
    }
    //printf("Finding non-wetting phase connectivity...");
    if (g->num_dims == 2)
        findConnectivity2d(p, g);
    else if (g->num_dims == 3)
        findConnectivity3d(p, g);
    //printf("Done\n");

    /* The main component index is the index of the inlet/outlet */
    /* The main component is the one connected to the inlet/outlet */
    main_comp_ind = o->phi_nw_ind;
    main_comp_val = p->connectivity[main_comp_ind];

    /* Get a binary matrix, where cells which belong to a disconnected component are zero,
       and everything else is 1 */    
    GET_BIN_VELOCITY(bin_vel_nw, p->connectivity, g, main_comp_val);
    
    
    //sprintf(fname, "phi_nw_bin");
    //writeDataArrayUchar(p->phi_bin, g, fname, GZIP);       
    //sprintf(fname, "phi_nw_conn");
    //writeDataArrayInt(p->connectivity, g, fname, GZIP);
    
    SET_DATA_TO_CONSTANT(p->connectivity, g, 0)
    
    /*--------------------------------------------------------------*/
    /* Now do it for wetting phase */
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);

    /* Generate binary wetting phase */
    IMPOSE_INT_EPS(p->phi_bin, p->scratch1, g, eps);            
            
    /* add a layer of 1's in slice ihi_fb + 1, to simulate connectivity to outlet */
    /* That way, any w blob connected to outlet should become part of  w-connected phase */
    if (g->num_dims == 3) {
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                 i = g->ihi_fb + 1;
                 idx_gb = i + j*nx_gb + k*nxy_gb;
                 p->phi_bin[idx_gb] = 1; // make slice 1 
            }
        } 
    
    }
    
    /* Find connectivity for wetting phase */  
    //printf("Finding wetting phase connectivity...");          
    if (g->num_dims == 2)
        findConnectivity2d(p, g);
    else if (g->num_dims == 3)
        findConnectivity3d(p, g);
    //printf("Done\n");
    
    main_comp_ind = o->phi_w_ind;
    main_comp_val = p->connectivity[main_comp_ind];
    
    GET_BIN_VELOCITY(bin_vel_w, p->connectivity, g, main_comp_val);
    
    
    //sprintf(fname, "phi_w_bin");
    //writeDataArrayUchar(p->phi_bin, g, fname, GZIP);
    //sprintf(fname, "phi_w_conn");
    //writeDataArrayInt(p->connectivity, g, fname, GZIP);
    
    /* Perform AND operation between binary w and nw to get a general matrix */
    IMPOSE_AND(bin_velocity, bin_vel_w, bin_vel_nw, g);
    //sprintf(fname, "bin_velocity");
    //writeDataArrayInt(bin_velocity, g, fname, GZIP);
    
    /* Setting velocities around the disconnected components to zero */    
    IMPOSE_CONN(p->normal_velocity, bin_velocity, g, val);
    //IMPOSE_CONN(p->curvature_coeff, bin_velocity, g, 0);
    //IMPOSE_CONN(p->external_velocity_x, bin_velocity, g, 0);
    //IMPOSE_CONN(p->external_velocity_y, bin_velocity, g, 0);
    
    if(g->num_dims == 3)
        IMPOSE_CONN(p->external_velocity_z, bin_velocity, g, 0);
    
    
    free(bin_velocity);
    free(bin_vel_w);
    free(bin_vel_nw);

}

void trapComponents_mask(QSS_DataArrays *p, Grid *g, Options *o){
    char fname[256];
    int main_comp_val, main_comp_ind;
    int idx_gb, nx_gb, nxy_gb, i, j, k;
    
    nx_gb = g->grid_dims_ghostbox[0];
    nxy_gb = g->grid_dims_ghostbox[1]*nx_gb;
    unsigned char *bin_disconn = (unsigned char *)malloc((g->num_gridpts)*sizeof(unsigned char));

    QSSLIB_REAL eps = 0;

    MERGE_SETS(p->phi, p->mask_nw, g);
    
    /* Find connectivity of non-wetting phase first */
    IMPOSE_INT_EPS(p->phi_bin, p->phi, g, eps);

    /* add a layer of 1's in slice ilo_fb - 1, to simulate connectivity to reservoir */
    /* That way, any nw blob connected to inlet should become part of  nw-connected phase */
    if (g->num_dims == 3) {
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                 i = g->ilo_fb - 1;
                 idx_gb = i + j*nx_gb + k*nxy_gb;
                 p->phi_bin[idx_gb] = 1; // make slice 1 
            }
        } 
    } 
    
    
    if (g->num_dims == 2)
        findConnectivity2d(p, g);
    else if (g->num_dims == 3)
        findConnectivity3d(p, g);

    /* The main component index is the index of the inlet/outlet
     The main component is the one connected to the inlet/outlet */
    main_comp_ind = o->phi_nw_ind;
    main_comp_val = p->connectivity[main_comp_ind];
   
    /* Get a binary matrix, where cells which belong to a disconnected component are 1,
       and everything else is 0 */    
    GET_BIN_DISCONN(bin_disconn, p->connectivity, g, main_comp_val);  
    
    /* Reset mask_nw */
   // COPY_DATA(p->mask_nw, p->mask_disconn_init, g); 
   /* It is set to -20*dx as that should be approximately the size of a pore,
    and should not affect things near the interface much */
    SET_DATA_TO_CONSTANT(p->mask_nw, g, -20*g->dx[0]);
    /* Rebuild mask_nw */
    MAKE_DISCONN_MASK(p->mask_nw, bin_disconn, p->phi, g);
        
    
    SET_DATA_TO_CONSTANT(p->connectivity, g, 0)
    
    /*--------------------------------------------------------------*/
    /* Now do it for wetting phase */
    /* Copy wetting phase into scratch1, to generate a masked copy of the wetting phase */
    COPY_DATA(p->scratch1, p->phi, g);
    NEGATE_DATA(p->scratch1, g);
    IMPOSE_MASK(p->scratch1, p->mask, p->scratch1, g);
    
    /* This may be unnecessary */
    IMPOSE_MASK(p->scratch1, p->mask_nw, p->scratch1, g);

    /* Since all of the wetting phase is in scratch1, and is masked, merging is unnecessary */
    /* Generate binary wetting phase */
    IMPOSE_INT_EPS(p->phi_bin, p->scratch1, g, eps);            
    
    /* add a layer of 1's in slice ihi_fb + 1, to simulate connectivity to outlet */
    /* That way, any w blob connected to outlet should become part of  w-connected phase */
    if (g->num_dims == 3) {
        for( k = g->klo_fb; k <= g->khi_fb; k++)
        {
            for( j = g->jlo_fb; j <= g->jhi_fb; j++)
            {
                 i = g->ihi_fb + 1;
                 idx_gb = i + j*nx_gb + k*nxy_gb;
                 p->phi_bin[idx_gb] = 1; // make slice 1 
            }
        } 
    
    } else {
    
    /* for Juanes 2d case, add wetting phase on all four sides of domain */
        for( j = g->jlo_gb; j <= g->jhi_gb; j++)
        {
            i = g->ilo_fb - 1;
            idx_gb = i + j*nx_gb;
            p->phi_bin[idx_gb] = 1; 
        }
        for( j = g->jlo_gb; j <= g->jhi_gb; j++)
        {
            i = g->ihi_fb + 1;
            idx_gb = i + j*nx_gb;
            p->phi_bin[idx_gb] = 1; 
        }
        
        for( i = g->ilo_gb; i <= g->ihi_gb; i++)
        {
            j = g->jlo_fb - 1;
            idx_gb = i + j*nx_gb;
            p->phi_bin[idx_gb] = 1; 
        }
        
        for( i = g->ilo_gb; i <= g->ihi_gb; i++)
        {
            j = g->jhi_fb + 1;
            idx_gb = i + j*nx_gb;
            p->phi_bin[idx_gb] = 1; 
        }
    }
    
    /* Find connectivity for wetting phase */  
    if (g->num_dims == 2)
        findConnectivity2d(p, g);
    else if (g->num_dims == 3)
        findConnectivity3d(p, g);
    
    main_comp_ind = o->phi_w_ind;
    main_comp_val = p->connectivity[main_comp_ind];
    
    GET_BIN_DISCONN(bin_disconn, p->connectivity, g, main_comp_val);
    
    //sprintf(fname, "disconn_wet");
    //writeDataArrayInt(p->connectivity, g, fname, GZIP);
    /* Reset mask_w */
    //initializeDisconnectedMasks(p->mask_w, g);
    //COPY_DATA(p->mask_w, p->mask_disconn_init, g); 
    SET_DATA_TO_CONSTANT(p->mask_w, g, -20*g->dx[0]);
    
        /* Rebuild mask_w */
    MAKE_DISCONN_MASK(p->mask_w, bin_disconn, p->scratch1, g);
    
//    qss_reinitializeDisconnectedMask2d(p, g, o);
    
    /*
    sprintf(fname, "phi_w_bin");
    writeDataArrayUchar(p->phi_bin, g, fname, GZIP);
    sprintf(fname, "phi_w_conn");
    writeDataArrayInt(p->connectivity, g, fname, GZIP);
    */
    
    free(bin_disconn);

}

/******************************************************************************
 *
 *   Author:   Masa Prodanovic
 *   Copyright (c) 2009, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file connectivity.c
    File contains routines for examining connectivity (number of disconnected 
    components etc.) of desired phase of a level set function (phi>0 or phi<0).
    
    The functions here operate on the digitized version of the phase so 
    connectivity is based on the way voxels are connected in digitized space:
       in 2D - 4 and 8 connectivity is implemented
       in 3D - 6 and 26 connectivity is implemented
*/
//#include "lsmlib_headers.h"
#include <stdlib.h>
#include <glib.h>
#include <math.h>
#include "connectivity.h"

/* Various utilities for finding voxels on volume boundary */

void find_max_conn_comp2d(
     char    dir,
     int     pos,
     QSSLIB_REAL   *phi,
     Grid   *g,
     int    *max_conn_sz,
     int    *min_conn_sz)
{
    int i,j, sz, max, min, idx;
    
    max = 0;
    min = g->num_gridpts;
    sz = 0;
    
    if(dir == 'x' )
    {
       
       for(j= 0; j < (g->grid_dims_ghostbox)[1];j++ )
       {
           idx = pos + j*(g->grid_dims_ghostbox)[0];
	   if( phi[idx] < 0 )
	   {
	      sz++;
	   }
	   else 
	   {
	       if(sz > max) max = sz;
	       if(sz > 0 && sz < min) min = sz;
	       //printf("\nsz %d",sz);
	       sz = 0; //reset
	   }    
       }
    }
    else
    {
       for(i= 0; i < (g->grid_dims_ghostbox)[0]; i++ )
       {
           idx = i + pos*(g->grid_dims_ghostbox)[0];
	   if( phi[idx] < 0 )
	   {
	      sz++;
	   }
	   else 
	   {
	       if(sz > max) max = sz;
	       if(sz > 0 && sz < min) min = sz; 
	       sz = 0; //reset
	   }    
       }
    }
    
    
    *max_conn_sz = max;
    *min_conn_sz = min;
}

void find_max_conn_comp3d(
     char    dim,
     int     pos,
     QSSLIB_REAL *phi,
     Grid   *g,
     int    *max_conn_sz,
     int    *min_conn_sz)
{
    int  i, j, k, sz, max, min, num, m, n, mn;
    int  nxy = (g->grid_dims_ghostbox)[0]*(g->grid_dims_ghostbox)[1], nx = (g->grid_dims_ghostbox)[0];
    int stencil[4], is_in[4];
    unsigned char *slice;
    int   idx, idx1, idx2, l;
    int   *new_list_item, *new_list_item1, *new_list_item2;
    GList *burn_list = NULL, *burn_list_ptr;
      
    num = 0;
    if(dim == 'x' )
    {
       m = (g->grid_dims_ghostbox)[1];
       n = (g->grid_dims_ghostbox)[2];
       mn = m*n;
       
       slice = (unsigned char *)calloc(mn,sizeof(unsigned char));
       
       idx1 = 0;
       for(k = 0; k < n; k++)
       {      
	 for(j= 0; j < m;j++,idx1++)
	 {    
	     idx = pos + j*nx + k*nxy;    
	     if( phi[idx] < 0 )
	     {
		slice[idx1] = UNBURN;
		num++;
	     }
	 }
       }
    }
    else if(dim == 'y' )
    {
       m = (g->grid_dims_ghostbox)[0];
       n = (g->grid_dims_ghostbox)[2];
       mn = m*n;
       
       slice = (unsigned char *)calloc(mn,sizeof(unsigned char));
       
       idx1 = 0;
       for(k = 0; k < n; k++)
       {      
	 for(i= 0; i < m;j++,idx1++)
	 {    
	     idx = i + pos*nx + k*nxy;    
	     if( phi[idx] < 0 )
	     {
		slice[idx1] = UNBURN;
		num++;
	     }
	 }
       }
    }
    else //(dim == 'z' )
    {
       m = (g->grid_dims_ghostbox)[0];
       n = (g->grid_dims_ghostbox)[1];
       mn = m*n;
       
       slice = (unsigned char *)calloc(mn,sizeof(unsigned char));
       
       idx1 = 0;
       for(j = 0; j < n; j++)
       {      
	 for(i= 0; i < m;i++,idx1++)
	 {    
	     idx = i + j*nx + pos*nxy;    
	     if( phi[idx] < 0 )
	     {
		slice[idx1] = UNBURN;
		num++;
	     }
	 }
       }
    }
    //printf("\ntotal unburn %d",num);
    set_4_stencil(&m,stencil);   
  
    max = 0;   min = g->num_gridpts;
    for(idx = 0; idx < mn; idx ++)
    {
        if( slice[idx] == UNBURN )
	{
          //initiate list with idx
	   sz = 1;
	   slice[idx] = TMPBURN;
	   new_list_item = (int *)malloc(sizeof(int));
	   *new_list_item = idx;
	   burn_list = g_list_prepend(burn_list,new_list_item);

	   while(sz)
	   {
	      burn_list_ptr = g_list_first(burn_list);
	      new_list_item1 = (burn_list_ptr->data);
	      idx1 = *new_list_item1;
	      //printf("\nidx1 %d",idx1);
	      check_if_4_neighbor_in_volume(idx1,m,n,1,mn,is_in);
	      
	      for(l = 0; l < 4; l++)
	      {
	          if(is_in[l])
		  {
	            idx2 = idx1 + stencil[l];
	            if( slice[idx2] == UNBURN )
		    {
		       slice[idx2] = TMPBURN;
		       new_list_item2 = (int *)malloc(sizeof(int));
	               *new_list_item2 = idx2;
		       burn_list = g_list_append(burn_list,new_list_item2);
		    }
		  } 
	      }
	      //remove examined item from the list
	      burn_list = g_list_remove(burn_list,new_list_item1);
	      free(new_list_item1);
	      
	      sz = g_list_length(burn_list);
	   }
	   
	   sz = 0;
	   for(i = 0; i < mn; i++)
	   {
	       if( slice[i] == TMPBURN ) 
	       {
	           sz++; slice[i] = BURN;
	       }
	   }
	   //printf("\ncomponent size %d",sz);
	   if( sz > max ) max = sz;
	   if(sz < min ) min = sz;
	}
    }
    
    *max_conn_sz = max;
    *min_conn_sz = min;
    
    g_list_free(burn_list);
        //Funky: sometimes I get abort errors if i free the memory?!
    free(slice);
}


int find_opp_side2d(
    QSSLIB_REAL    *phi,
    QSSLIB_REAL    *mask,
    Grid           *g,
    int            **pvox_ind)
{

     int i, j, nvox, *vox_ind, idx;
     int nx = (g->grid_dims_ghostbox)[0];
     int ny = (g->grid_dims_ghostbox)[1];
     int n_alloc = 2*( nx + ny);
     QSSLIB_REAL zero=0.0;
     
     nvox = 0;
     
     vox_ind = (int *)malloc(n_alloc*sizeof(int));
     
     for(j = 0; j < ny; j++ )
     {
           i = 0;
           idx = i + j*nx;
	   if(( phi[idx] > zero ) && ( mask[idx] < zero ))
	   {
	      vox_ind[nvox++] = idx;
	     // printf("\n(%d,%d)",i,j);	      
	   }
	  
	   i = nx-1;
	   idx = i + j*nx;
	   if(( phi[idx] > zero) && ( mask[idx] < zero))
	   {
	      vox_ind[nvox++] = idx;
	      // printf("\n(%d,%d)",i,j);
	      
	   }
     }
     
     for(i= 1; i < nx-1;i++ )
     {
           j = 0;
           idx = i; 
	   if(( phi[idx] > zero ) && ( mask[idx] < zero ))
	   {
	      vox_ind[nvox++] = idx;
	       //printf("\n(%d,%d)",i,j);
	   }
	  
	   j = ny-1;
	   idx = i + j*nx;
	   if(( phi[idx] > 0 ) && ( mask[idx] < zero ))
	   {
	      vox_ind[nvox++] = idx;
	      //printf("\n(%d,%d)",i,j);
	   }
     }
     
     vox_ind = (int *)realloc(vox_ind,nvox*sizeof(int));
     
     *pvox_ind = vox_ind;
     return nvox;
}


int find_opp_side3d(
    QSSLIB_REAL *phi,
    QSSLIB_REAL *mask,
    Grid *g,
    int **pvox_ind)
{

     int i, j, k,  nvox, *vox_ind, idx;
     int nxy = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[1], nxz = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[2];
     int nyz = (g->grid_dims_ghostbox)[1] * (g->grid_dims_ghostbox)[2];
     int  n_alloc = 2*(nxy + nyz + nxz);
     nvox = 0;
     
     vox_ind = (int *)malloc(n_alloc*sizeof(int));
     
     for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     {
       for(j= 0; j < (g->grid_dims_ghostbox)[1];j++ )
       {
             idx = 0 + j*(g->grid_dims_ghostbox)[0] + k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }

	     idx = (g->grid_dims_ghostbox)[0]-1 + j*(g->grid_dims_ghostbox)[0]+ k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
       }
     }
     
     for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     {
       for(i= 1; i < (g->grid_dims_ghostbox)[0]-1;i++ )
       {
             idx = i + 0 + k*nxy ;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }

	     idx = i + ((g->grid_dims_ghostbox)[1]-1)*(g->grid_dims_ghostbox)[0] + k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
       }
     }
     
     for(j = 1; j < (g->grid_dims_ghostbox)[1]-1; j++)
     {
       for(i= 1; i < (g->grid_dims_ghostbox)[0]-1;i++ )
       {
             idx = i + j*(g->grid_dims_ghostbox)[0] + 0*nxy ;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }

	     idx = i + j*(g->grid_dims_ghostbox)[0] + ( (g->grid_dims_ghostbox)[2] - 1 )*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
       }
     }
     
     vox_ind = (int *)realloc(vox_ind,nvox*sizeof(int));
     
     *pvox_ind = vox_ind;
     return nvox;
}


int  find_opp_side3d_x(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int **pvox_ind)
{

     int  j, k,  nvox, *vox_ind, idx;
     int nxy = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[1];
     int nxz = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[2];
     int nyz = (g->grid_dims_ghostbox)[1] * (g->grid_dims_ghostbox)[2];
     int  n_alloc = 2*(nxy + nyz + nxz);
     nvox = 0;
     
     vox_ind = (int *)malloc(n_alloc*sizeof(int));
     
     for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     {
       for(j= 0; j < (g->grid_dims_ghostbox)[1]; j++ )
       {

             //TEMP - skip looking at x=0 side!
/*
             idx = 0 + j*(g->grid_dims_ghostbox)[0] + k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
*/
	     idx = (g->grid_dims_ghostbox)[0]-1 + j*(g->grid_dims_ghostbox)[0]+ k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
       }
     }
     
     vox_ind = (int *)realloc(vox_ind,nvox*sizeof(int));
     
     *pvox_ind = vox_ind;
     return nvox;
}



int  find_yz_plane_voxels(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int position,
     int **pvox_ind)
{

     int  j, k,  nvox, *vox_ind, idx;
     int nxy = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[1];
     int nxz = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[2];
     int nyz = (g->grid_dims_ghostbox)[1] * (g->grid_dims_ghostbox)[2];
     int  n_alloc = 2*(nxy + nyz + nxz);
     nvox = 0;
     
     vox_ind = (int *)malloc(n_alloc*sizeof(int));
     
     for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     {
       for(j= 0; j < (g->grid_dims_ghostbox)[1]; j++ )
       {
             idx = position + j*(g->grid_dims_ghostbox)[0] + k*nxy;
	     if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	     {
		vox_ind[nvox++] = idx;
	     }
       }
     }
     
     vox_ind = (int *)realloc(vox_ind,nvox*sizeof(int));
     
     *pvox_ind = vox_ind;
     return nvox;
}


int  find_multiple_yz_plane_voxels(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int position_start,
     int position_end,
     int **pvox_ind)
{

     int  i,j, k,  nvox, *vox_ind, idx;
     int nxy = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[1];
     int nxz = (g->grid_dims_ghostbox)[0] * (g->grid_dims_ghostbox)[2];
     int nyz = (g->grid_dims_ghostbox)[1] * (g->grid_dims_ghostbox)[2];
     int  n_alloc = 2*(nxy + nyz + nxz);
     nvox = 0;
     
     vox_ind = (int *)malloc(n_alloc*sizeof(int));
     
     for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     {
       for(j= 0; j < (g->grid_dims_ghostbox)[1]; j++ )
       {
             for(i = position_start; i <=position_end; i++)
	     {
               idx = i + j*(g->grid_dims_ghostbox)[0] + k*nxy;
	       if(( phi[idx] > 0 ) && ( mask[idx] < 0 ))
	       {
		  vox_ind[nvox++] = idx;
	       }
	     }  
       }
     }
     
     vox_ind = (int *)realloc(vox_ind,nvox*sizeof(int));
     
     *pvox_ind = vox_ind;
     return nvox;
}   


void	set_4_stencil(
	int	*n,
	int	*sten)
{	
        int    nx = n[0];
	
			sten[3] =   nx;
	sten[0] = -1;			sten[1] = 1;
			sten[2] =  -nx;
}


void	set_8_stencil(
	int	*n,
	int	*sten)
{
        int    nx = n[0];
	
	sten[0] = -nx-1;      sten[1] = -nx;      sten[2] = -nx+1;
	sten[3] = -1;			            sten[4] =  1;
	sten[5] =  nx-1;      sten[6] =  nx;      sten[7] =  nx+1;
}


void	set_6_stencil(
	int	*n,
	int	*sten)
{
	int    nx = n[0], nxy = n[0]*n[1];
	
			sten[5] =  nxy;

			sten[3] =   nx;
	sten[0] = -1;			sten[1] = 1;
			sten[2] =  -nx;

			sten[4] = -nxy;
}


void	set_26_stencil(
	int	*n,
	int	*sten)
{
        int     nx = n[0], nxy = n[0]*n[1];
	int	nxm1, nxp1;

	nxm1 = nx-1;		nxp1 = nx+1;

	sten[ 0] = -nxy-nxp1;  sten[ 1] = -nxy-nx;  sten[ 2] = -nxy-nxm1;
	sten[ 3] = -nxy-1;     sten[ 4] = -nxy;     sten[ 5] = -nxy+1;
	sten[ 6] = -nxy+nxm1;  sten[ 7] = -nxy+nx;  sten[ 8] = -nxy+nxp1;
	sten[ 9] = -nxp1;      sten[10] = -nx;      sten[11] = -nxm1;
	sten[12] = -1;			            sten[13] =  1;
	sten[14] =  nxm1;      sten[15] =  nx;      sten[16] =  nxp1;
	sten[17] =  nxy-nxp1;  sten[18] =  nxy-nx;  sten[19] =  nxy-nxm1;
	sten[20] =  nxy-1;     sten[21] =  nxy;     sten[22] =  nxy+1;
	sten[23] =  nxy+nxm1;  sten[24] =  nxy+nx;  sten[25] =  nxy+nxp1;
}


void	check_if_4_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume)
{
	int	i;
	int	ind, ix, iy;

	for( i = 0;  i < 4;  i++ ) is_in_volume[i] = 1;


	ind = pos;	iy = ind/nx;    ix = ind%nx;

	/* Note: if() else if() NOT desired in case ny == 1, nx == 1 */

	if( iy == ny-1 ) { is_in_volume[3] = 0; }
	if( iy ==    0 ) { is_in_volume[2] = 0; }

	if( ix == nx-1 ) { is_in_volume[1] = 0; }
	if( ix ==    0 ) { is_in_volume[0] = 0; }
}


void	check_if_8_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume)
{
	int	i;
	int	ind, ix, iy;

	for( i = 0;  i < 8;  i++ ) is_in_volume[i] = 1;


	ind = pos;	iy = ind/nx;    ix = ind%nx;

	/* Note: if() else if() NOT desired in case ny == 1, nx == 1 */

	if( iy == ny-1 ) 
	{ 
	    is_in_volume[5] = is_in_volume[6] = is_in_volume[7] = 0; 
	}
	if( iy ==    0 ) 
	{ 
	    is_in_volume[0] = is_in_volume[1] = is_in_volume[2] = 0; 
	}

	if( ix == nx-1 ) 
	{ 
	    is_in_volume[2] = is_in_volume[4] = is_in_volume[7] = 0; 
	}
	
	if( ix ==    0 )
	{ 
	    is_in_volume[0] = is_in_volume[3] = is_in_volume[5] = 0; 
	}
}


void	check_if_6_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume)
{
	int	i;
	int	ind, ix, iy, iz;

	for( i = 0;  i < 6;  i++ ) is_in_volume[i] = 1;


	ind = pos;		iz = ind/nxy;
	ind = ind%nxy;		iy = ind/nx;
				ix = ind%nx;

	/* Note: if() else if() not desired in case nz == 1, ny == 1, nx == 1 */

	if( iz == nz-1 ) { is_in_volume[5] = 0; }
	if( iz ==    0 ) { is_in_volume[4] = 0; }

	if( iy == ny-1 ) { is_in_volume[3] = 0; }
	if( iy ==    0 ) { is_in_volume[2] = 0; }

	if( ix == nx-1 ) { is_in_volume[1] = 0; }
	if( ix ==    0 ) { is_in_volume[0] = 0; }
}


void	check_if_26_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume)
{
	int	i;
	int	ind, ix, iy, iz;

	for( i = 0;  i < 26;  i++ ) is_in_volume[i] = 1;

	ind = pos;		iz = ind/nxy;
	ind = ind%nxy;		iy = ind/nx;
				ix = ind%nx;

	/* Note: if() else if() not desired in case nz = 1, ny = 1, nx = 1 */

	if( iz == nz-1 )
	{
		for( i = 17;  i < 26;  i++ ) is_in_volume[i] = 0;
	}
	if( iz == 0  )
	{
		for( i =  0;  i <  9;  i++ ) is_in_volume[i] = 0;
	}

	if( iy == ny-1 )
	{
		for( i =  6;  i <  9;  i++ ) is_in_volume[i] = 0;
		for( i = 14;  i < 17;  i++ ) is_in_volume[i] = 0;
		for( i = 23;  i < 26;  i++ ) is_in_volume[i] = 0;
	}
	if( iy == 0 )
	{
		for( i =  0;  i <  3;  i++ ) is_in_volume[i] = 0;
		for( i =  9;  i < 12;  i++ ) is_in_volume[i] = 0;
		for( i = 17;  i < 20;  i++ ) is_in_volume[i] = 0;
	}

	if( ix == nx-1 )
	{
		for( i =  2;  i < 12;  i += 3 ) is_in_volume[i] = 0;
		for( i = 13;  i < 26;  i += 3 ) is_in_volume[i] = 0;
	}
	if( ix == 0 )
	{
		for( i =  0;  i < 13;  i += 3 ) is_in_volume[i] = 0;
		for( i = 14;  i < 26;  i += 3 ) is_in_volume[i] = 0;
	}
}


int findPhaseConnectedComponents(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size)
{
    int   i, j, k, set_size, num_comp;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int           stencil[stencil_size], is_in[stencil_size];
    unsigned char *data;
    int           idx, idx1, idx2, idx3, l, *new_list_item, *new_list_item2, list_size;
    int           *component_size, component_size_alloc = 20, component_size_realloc = 20;
    GList         *burn_list = NULL, *burn_list_ptr;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
    
    set_size = 0;
    
    data = (unsigned char *)calloc(g->num_gridpts,sizeof(unsigned char));
    
    /* Initialize appropriate phase voxels to UNBURN value */
    for(i = 0; i < g->num_gridpts; i++)
    {
       if( phi[i]*phase > 0 ) 
       {  
          data[i] = UNBURN;
	  set_size++;
       }	  
    }
    
    //printf("\nTotal %d voxels in the set",set_size); fflush(stdout);
    
    SET_STENCIL_FUNCTIONS()
    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
            
    num_comp = 0; /* number of components */
    component_size = (int *)calloc(component_size_alloc,sizeof(int));
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( data[idx] == UNBURN )
	  {  
             /* initiate list with voxel 'idx' */
	     list_size = 1;
	     data[idx] = TMPBURN;
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);

             /* loop until there are no more voxels in the list */
	     while(list_size)
	     {
	        /* get the first voxel from the list */
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);                
			
		/* loop over all neighbors of the first voxel*/		
		for(l = 0; l < stencil_size; l++)
		{
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		     *  to the list and mark with temporary mark TMPBURN
		     */  
	           if( ( is_in[l] ) && ( data[idx2] == UNBURN ) )
		   {
			data[idx2] = TMPBURN;
			
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
		   }
		}
		
		/* remove the examined voxel from the list */
		burn_list = g_list_remove(burn_list,new_list_item);
		free(new_list_item);
		list_size = g_list_length(burn_list);
	     }

             /* replace temporary mark TMPBURN by the permanent one (BURN) */
	     component_size[num_comp] = 0;
	     for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( data[idx3] == TMPBURN ) 
		 {
	             component_size[num_comp]++;
		     data[idx3] = BURN;
		 }
	     }
	     
	     /* increase count the connected components */
	     num_comp++;
	     
	     /* re-allocate memory if needed */ 
	     if( num_comp > component_size_alloc)
	     {
	       component_size_alloc += component_size_realloc;
	       component_size =(int*)realloc(component_size,component_size_alloc*sizeof(int));	     
	     }
	  }
       }
        
    component_size =(int*)realloc(component_size,num_comp*sizeof(int));
     	
    free(data);
    g_list_free(burn_list);
    
    *pcomponent_size = component_size;
    return num_comp;
}



int findPhaseConnectedComponentsControlVolume(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     QSSLIB_REAL *control_vol,
     QSSLIB_REAL control_vol_phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size)
{
    int   i, j, k, set_size, num_comp;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int           stencil[stencil_size], is_in[stencil_size];
    unsigned char *data;
    int           idx, idx1, idx2, idx3, l, *new_list_item, *new_list_item2, list_size;
    int           *component_size, component_size_alloc = 20, component_size_realloc = 20;
    GList         *burn_list = NULL, *burn_list_ptr;
    QSSLIB_REAL        zero = 0.0;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
    
    set_size = 0;
    
    data = (unsigned char *)calloc(g->num_gridpts,sizeof(unsigned char));
    
    /* Initialize appropriate phase voxels to UNBURN value 
      TO_DO - this multiplication business is not the fastest way to do it, speed up*/
    for(i = 0; i < g->num_gridpts; i++)
    {
       if( phi[i]*phase > zero && control_vol[i]*control_vol_phase > zero) 
       {  
          data[i] = UNBURN;
	  set_size++;
       }	  
    }
    
    //printf("\nTotal %d voxels in the set",set_size); fflush(stdout);
    
    SET_STENCIL_FUNCTIONS()

    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
            
    num_comp = 0; /* number of components */
    component_size = (int *)calloc(component_size_alloc,sizeof(int));
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( data[idx] == UNBURN )
	  {  
             /* initiate list with voxel 'idx' */
	     list_size = 1;
	     data[idx] = TMPBURN;
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);

             /* loop until there are no more voxels in the list */
	     while(list_size)
	     {
	        /* get the first voxel from the list */
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);                
			
		/* loop over all neighbors of the first voxel*/		
		for(l = 0; l < stencil_size; l++)
		{
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		     *  to the list and mark with temporary mark TMPBURN
		     */  
	           if( ( is_in[l] ) && ( data[idx2] == UNBURN ) )
		   {
			data[idx2] = TMPBURN;
			
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
		   }
		}
		
		/* remove the examined voxel from the list */
		burn_list = g_list_remove(burn_list,new_list_item);
		free(new_list_item);
		list_size = g_list_length(burn_list);
	     }

             /* replace temporary mark TMPBURN by the permanent one (BURN) */
	     component_size[num_comp] = 0;
	     for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( data[idx3] == TMPBURN ) 
		 {
	             component_size[num_comp]++;
		     data[idx3] = BURN;
		 }
	     }
	     
	     /* increase count the connected components */
	     num_comp++;
	     
	     /* re-allocate memory if needed */ 
	     if( num_comp > component_size_alloc)
	     {
	       component_size_alloc += component_size_realloc;
	       component_size =(int*)realloc(component_size,component_size_alloc*sizeof(int));	     
	     }
	  }
       }
        
    component_size =(int*)realloc(component_size,num_comp*sizeof(int));
     	
    free(data);
    g_list_free(burn_list);
    
    *pcomponent_size = component_size;
    return num_comp;
}


int processPhaseConnectedComponents(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     Grid   *g,
     int     stencil_size,
     int     min_size)
{
    int   i, j, k, set_size, num_comp;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    unsigned char *data;
    int   idx, idx1, idx2, idx3, l, *new_list_item, *new_list_item2, sz;
    GList *burn_list = NULL, *burn_list_ptr;
    
    int   *component_size, component_size_alloc = 20, component_size_realloc = 20;
    int   *comp_start_vox, max_component_size, max_comp_id;
    
    QSSLIB_REAL opposite_phase_val;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
  
    set_size = 0;
    
    data = (unsigned char *)calloc(g->num_gridpts,sizeof(unsigned char));
    
    /* Initialize appropriate phase voxels to UNBURN value */
    for(i = 0; i < g->num_gridpts; i++)
    {
       if( phi[i]*phase > 0 ) 
       {  
          data[i] = UNBURN;
	  set_size++;
       }	  
    }
    
    opposite_phase_val = -phase*g->dx[0];
    //printf("\nTotal %d voxels in the set",set_size); fflush(stdout);
    
    SET_STENCIL_FUNCTIONS()


    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  
            
    num_comp = 0; /* number of components */
    component_size = (int *)calloc(component_size_alloc,sizeof(int));
    comp_start_vox = (int *)calloc(component_size_alloc,sizeof(int));
    
    max_component_size = 0;
    max_comp_id = 0;
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( data[idx] == UNBURN )
	  {  
             /* initiate list with voxel 'idx' */
	     sz = 1;
	     data[idx] = TMPBURN;
	     comp_start_vox[num_comp] = idx;
	     
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);

             /* loop until there are no more voxels in the list */
	     while(sz)
	     {
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in); 
				
		for(l = 0; l < stencil_size; l++)
		{
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		       to the list */  
	           if( ( is_in[l] ) && ( data[idx2] == UNBURN ) )
		   {
			data[idx2] = TMPBURN;
			
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
		   }
		}
		/* remove voxel 'idx' from the list */
		burn_list = g_list_remove(burn_list,new_list_item);
		free(new_list_item);
		sz = g_list_length(burn_list);
	     }

	     component_size[num_comp] = 0;
	     for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( data[idx3] == TMPBURN ) 
		 {
	             component_size[num_comp]++;
		 }
	     }
	     
	     printf("\nsize %d",component_size[num_comp]);
	     if(component_size[num_comp] > max_component_size)
	     {
	       max_component_size = component_size[num_comp];
	       max_comp_id = num_comp;
	     }
	     
	     if( component_size[num_comp] < min_size )
	     {
	       //set voxels of this component to the opposite phase
	         for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	         {
		   if( data[idx3] == TMPBURN ) 
		   {
	             phi[idx3] = opposite_phase_val;
		   }
	         }	     
	     } 
	     
	      for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( data[idx3] == TMPBURN ) 
		 {
	            data[idx3] = BURN;
		 }
	     }
	     
	     num_comp++; /* increase count the connected components */
	     if( num_comp > component_size_alloc)
	     {
	       component_size_alloc += component_size_realloc;
	       component_size =(int*)realloc(component_size,component_size_alloc*sizeof(int));
	       comp_start_vox =(int*)realloc(component_size,component_size_alloc*sizeof(int));	     
	     }
	  }
       }
        
    component_size =(int*)realloc(component_size,num_comp*sizeof(int));
     	
    printf("\nmax_size %d total comp %d",max_component_size,num_comp);	
    free(data);
    g_list_free(burn_list);
    free(component_size);
    return num_comp;
}


int labelPhaseConnectedComponents(
     QSSLIB_REAL  *phi,
     QSSLIB_REAL  phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size,
     int     **pcomp_start_vox,
     int     **plabel_phi,
     int     min_label)
{
    int   i, j, k, set_size, num_comp;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   *label_phi;
    int   idx, idx1, idx2, idx3, l, *new_list_item, *new_list_item2, list_size;
    GList *burn_list = NULL, *burn_list_ptr;
    
    int   *component_size, component_size_alloc = 20, component_size_realloc = 20;
    int   *comp_start_vox, label;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
  
    set_size = 0;
    
    label_phi = (int *)calloc(g->num_gridpts,sizeof(int));
    
    /* Initialize appropriate phase voxels to UNBURN value */
    for(i = 0; i < g->num_gridpts; i++)
    {
       if( phi[i]*phase > 0 ) 
       {  
          label_phi[i] = INT_UNBURN;
	  set_size++;
       }	  
    }
    
    //printf("\nTotal %d voxels in the set",set_size); fflush(stdout);
        
    SET_STENCIL_FUNCTIONS()



    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  
            
    num_comp = 0; /* number of components */
    component_size = (int *)calloc(component_size_alloc,sizeof(int));
    comp_start_vox = (int *)calloc(component_size_alloc,sizeof(int));
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( label_phi[idx] == INT_UNBURN )
	  {  
	     label = num_comp + min_label;
	     
             /* initiate list with voxel 'idx' */
	     list_size = 1;	     	     
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);
	     
	     /* label the starting voxel */
	     label_phi[idx] = label;

             /* record the starting voxel for this component */
             comp_start_vox[num_comp] = idx;
	     

             /* loop until there are no more voxels in the list */
	     while(list_size)
	     {
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);
					
		/* loop over all neighborhood_voxels */		
		for(l = 0; l < stencil_size; l++)
		{		   
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		       to the list */  
	           if( ( is_in[l] ) && ( label_phi[idx2] == INT_UNBURN ) )
		   {			
			/* add newly acquired voxel to the list */
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
			
			/* label newly acquired voxel */
			label_phi[idx2] = label;

		   }
		}
		
		/* remove examined voxel from the list */
		burn_list = g_list_remove(burn_list,new_list_item);
		free(new_list_item);
		list_size = g_list_length(burn_list);
	     }

	     component_size[num_comp] = 0;
	     for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( label_phi[idx3] == label ) 
		 {
	             component_size[num_comp]++;
		 }
	     }
	     
	     /* increase count of the connected components */
	     num_comp++;
	     
	     if( num_comp > component_size_alloc)
	     {
	       component_size_alloc += component_size_realloc;
	       component_size =(int*)realloc(component_size,component_size_alloc*sizeof(int));
	       comp_start_vox =(int*)realloc(component_size,component_size_alloc*sizeof(int));	     
	     }
	  }
       }
        
    component_size =(int*)realloc(component_size,num_comp*sizeof(int));
    comp_start_vox =(int*)realloc(comp_start_vox,num_comp*sizeof(int)); 	
    g_list_free(burn_list);
    
    *plabel_phi = label_phi;
    *pcomponent_size = component_size;
    *pcomp_start_vox = comp_start_vox;
    
    return num_comp;
}


int labelPhaseConnectedComponentsControlVolume(
     QSSLIB_REAL  *phi,
     QSSLIB_REAL  phase,
     QSSLIB_REAL *control_vol,
     QSSLIB_REAL  control_vol_sign,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size,
     int     **pcomp_start_vox,
     int     **plabel_phi,
     int     min_label)
{
    int   i, j, k, set_size, num_comp;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   *label_phi;
    int   idx, idx1, idx2, idx3, l, *new_list_item, *new_list_item2, list_size;
    GList *burn_list = NULL, *burn_list_ptr;
    
    int   *component_size, component_size_alloc = 20, component_size_realloc = 20;
    int   *comp_start_vox, label;
    
    QSSLIB_REAL  zero = 0.0;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
  
    set_size = 0;
    
    label_phi = (int *)calloc(g->num_gridpts,sizeof(int));
    
    /* Initialize appropriate phase voxels to UNBURN value */
    for(i = 0; i < g->num_gridpts; i++)
    {
       if( (phi[i]*phase >= zero) && (control_vol[i]*control_vol_sign >= zero) ) 
       {  
          label_phi[i] = INT_UNBURN;
	  set_size++;
       }	  
    }
    
    //printf("\nTotal %d voxels in the set",set_size); fflush(stdout);
        
    SET_STENCIL_FUNCTIONS()




    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  
            
    num_comp = 0; /* number of components */
    component_size = (int *)calloc(component_size_alloc,sizeof(int));
    comp_start_vox = (int *)calloc(component_size_alloc,sizeof(int));
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( label_phi[idx] == INT_UNBURN )
	  {  
	     label = num_comp + min_label;
	     
             /* initiate list with voxel 'idx' */
	     list_size = 1;	     	     
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);
	     
	     /* label the starting voxel */
	     label_phi[idx] = label;

             /* record the starting voxel for this component */
             comp_start_vox[num_comp] = idx;
	     

             /* loop until there are no more voxels in the list */
	     while(list_size)
	     {
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);
					
		/* loop over all neighborhood_voxels */		
		for(l = 0; l < stencil_size; l++)
		{		   
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		       to the list */  
	           if( ( is_in[l] ) && ( label_phi[idx2] == INT_UNBURN ) )
		   {			
			/* add newly acquired voxel to the list */
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
			
			/* label newly acquired voxel */
			label_phi[idx2] = label;

		   }
		}
		
		/* remove examined voxel from the list */
		burn_list = g_list_remove(burn_list,new_list_item);
		free(new_list_item);
		list_size = g_list_length(burn_list);
	     }

	     component_size[num_comp] = 0;
	     for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	     {
		 if( label_phi[idx3] == label ) 
		 {
	             component_size[num_comp]++;
		 }
	     }
	     
	     /* increase count of the connected components */
	     num_comp++;
	     
	     if( num_comp > component_size_alloc)
	     {
	       component_size_alloc += component_size_realloc;
	       component_size =(int*)realloc(component_size,component_size_alloc*sizeof(int));
	       comp_start_vox =(int*)realloc(comp_start_vox,component_size_alloc*sizeof(int));	     
	     }
	  }
       }
        
    component_size =(int*)realloc(component_size,num_comp*sizeof(int));
    comp_start_vox =(int*)realloc(comp_start_vox,num_comp*sizeof(int)); 	
    g_list_free(burn_list);
    
    *plabel_phi = label_phi;
    *pcomponent_size = component_size;
    *pcomp_start_vox = comp_start_vox;
    
    return num_comp;
}


/*! \fn examineComponentNeighborhood() examines the neighborhood of a component labeled by
*   labelPhaseConnectedComponents() for presence of grain voxels (mask > 0) as well as
*   the exterior (i.e. whether the component touches volume boundary.
*   Components that neither neighbor grain nor exterior are assumed trapped.
*
*   Arguments:


*   output_QSSLIB_REAL_array(in) - 0 or 1 depending on whether a QSSLIB_REAL array with neg. values where
*                            the blob is should be output
*   basename(in)            - file basename (component_label will be added to the basename)
*/
  
void examineComponentNeighborhood(
    int     component_label,
    int     component_start_vox,
    int     *label_phi,
    QSSLIB_REAL   *phi,
    QSSLIB_REAL   *mask,
    Grid    *g,
    int     stencil_size, 
    int     *touching_exterior,
    int     *touching_grain,
    int     output_QSSLIB_REAL_array,
    char    *basename)
{    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   idx, idx1, idx2, l, *new_list_item, *new_list_item2, list_size;
    GList *burn_list = NULL, *burn_list_ptr;
    
    char  filename[256];
    QSSLIB_REAL *phi_comp;
    
     /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();
    
    SET_STENCIL_FUNCTIONS()




   
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  
    if(output_QSSLIB_REAL_array)
    {      
       phi_comp = (QSSLIB_REAL *)malloc(g->num_gridpts*sizeof(QSSLIB_REAL));
      
       for(idx=0; idx < g->num_gridpts; idx++) 
       {
          phi_comp[idx] = fabs( phi[idx] );
       }	  
    }

    *touching_grain = 0;
    *touching_exterior = 0;
    
    /* initiate list with voxel 'idx' */
    idx = component_start_vox;
    list_size = 1;	     	     
    new_list_item = (int *)malloc(sizeof(int));
    *new_list_item = idx;
    burn_list = g_list_prepend(burn_list,new_list_item);
    
    label_phi[idx] = INT_BURN;
    if(output_QSSLIB_REAL_array)
    {
       phi_comp[idx] = -fabs( phi[idx] );
    }
    
    /* loop until there are no more voxels in the list */
    while(list_size)
    {
       burn_list_ptr = g_list_first(burn_list);
       new_list_item = (burn_list_ptr->data);
       idx1 = *new_list_item;

       /* check if neighbors of the first voxel are in the volume */
       check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);
       
       

       /* loop over all neighborhood_voxels */		
       for(l = 0; l < stencil_size; l++)
       {		   
	  idx2 = idx1 + stencil[l];		     

	   /* if the neighbor is in the volume and UNBURN, add it
	      to the list */  
	  if( is_in[l]  )
	  {	
	       if( label_phi[idx2] == component_label )
	       {	
	         /* add newly acquired voxel to the list */
	         new_list_item2 = (int *)malloc(sizeof(int));
	         *new_list_item2 = idx2;
	         burn_list = g_list_append(burn_list,new_list_item2);
		 
		 label_phi[idx2] = INT_BURN;
		 if(output_QSSLIB_REAL_array)
		 {
		   phi_comp[idx2] = - fabs( phi[idx2] );
	         }
	       }
	       else if( mask[idx2] > 0)
	       {
	          *touching_grain = 1;
	       }
	  }
	  else *touching_exterior = 1;	  
       }

       /* remove examined voxel from the list */
       burn_list = g_list_remove(burn_list,new_list_item);
       free(new_list_item);
       list_size = g_list_length(burn_list);
    }
	   
    if(output_QSSLIB_REAL_array)
    {
        //QSSLIB_REAL volume_comp, eps = 1.5*g->dx[0];
	
        sprintf(filename,"%s%d",basename,component_label);
        writeDataArrayQSS(phi_comp,g,filename,GZIP);
	
	/*	   
         Volume computation probably not very accurate, need to reinitialize the array
	 (better off with counting voxels, if not reinitialized)
	    LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(&volume_comp,
	    phi_comp,
            &(g->ilo_gb), &(g->ihi_gb), &(g->jlo_gb), &(g->jhi_gb),
	    &(g->klo_gb), &(g->khi_gb),
            &(g->ilo_fb), &(g->ihi_fb), &(g->jlo_fb), &(g->jhi_fb),
            &(g->klo_fb), &(g->khi_fb),
            &(g->dx[0]),&(g->dx[1]),&(g->dx[2]),
	    &eps);
	printf("volume %g\n",volume_comp);*/
	free(phi_comp);
    }
}    



int numberLabeledComponentsInPhaseNeighborhood(
    QSSLIB_REAL  *phi,
    QSSLIB_REAL   phase,
    int    *labels_data,
    int    max_num_labels,
    unsigned char **plabel_found,
    Grid   *g,
    int    stencil_size)
{
      int   i, j, k;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   idx, idx1, l;
    int   num_labels_found;
   
    unsigned char   *label_found;
    
     /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();
    
     SET_STENCIL_FUNCTIONS()





    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity  sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  

    label_found = (unsigned char *)calloc(max_num_labels,sizeof(unsigned char));

    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( phi[idx]*phase > 0 )
	  {		
		/* check if neighbors of the voxel 'idx' are in the volume */
		check_if_neighbor_in_volume_func(idx,nx,ny,nz,nxy,is_in); 
				
		for(l = 0; l < stencil_size; l++)
		{
	           idx1 = idx + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		       to the list */  
	           if( ( is_in[l] ) && ( labels_data[idx1] > 0) )
		   {
			label_found[ labels_data[idx1] - 1] = 1;
		   }
		}
	   }
	 }
	 
    num_labels_found = 0;
    for(i = 0; i < max_num_labels; i++) num_labels_found += label_found[i];
     	 
    *plabel_found = label_found;	 
    return   num_labels_found;
}    

/*! \fn existVoxelsInPhaseNeighborhood() checks  neighborhood of the level set function
*  phase (phase1 of phi1) and returns 1 if a (different) level set function 
*  phase (phase2 of phi2) voxels exist in its neighborhood (and otherwise 0).
*
*  Arguments
*    phi1 - level set function array
*    phase1 - phase of phi1 to examine (positive or negative)
*    phi2 - level set function array
*    phase2 - phase of phi2 to consider (positive or negative)
*    labels_data - array with integer labels in [1,max_num_labels]
*    max_num_labels - max number of integer labels
*    plabel_found - pointer to newly allocated array of max_num_labels elements
*                   *plabel_found[i] = 1 if label(component) i+1 is identified 
*                   in the neighborhood of the phase
*    g - Grid structure
*    stencil_size - stencil size, i.e. pixel/voxel connectivity (4,8,6 or 26) 
*    
*/
  
int existVoxelsInPhaseNeighborhood(
    QSSLIB_REAL  *phi1,
    QSSLIB_REAL   phase1,
    QSSLIB_REAL  *phi2,
    QSSLIB_REAL   phase2,
    Grid   *g,
    int    stencil_size)
{
      int   i, j, k;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   idx, idx1, l;
    int   voxels_found = 0;
   
     /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();
    
     SET_STENCIL_FUNCTIONS()






    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity  sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
  

    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( phi1[idx]*phase1 > 0 )
	  {		
		/* check if neighbors of the voxel 'idx' are in the volume */
		check_if_neighbor_in_volume_func(idx,nx,ny,nz,nxy,is_in); 
				
		for(l = 0; l < stencil_size; l++)
		{
	           idx1 = idx + stencil[l];		     
		     
		    /* if the neighbor is in the desired phase
		       set the voxels_found to 1 and return*/  
	           if( ( is_in[l] ) && ( phi2[idx1]*phase2 > 0) )
		   {
			voxels_found = 1;
			return voxels_found;
		   }
		}
	   }
	 }
	 
   return   voxels_found;
}    

/* NOT USED AT THE MOMENT */
QSSLIB_REAL  *expandPositiveSet(QSSLIB_REAL *phi,Grid *g,int stencil_size)
{
    int   i, j, k;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int   stencil[stencil_size], is_in[stencil_size];
    int   idx, idx1, l;
    
    QSSLIB_REAL *phi1;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();
    
    SET_STENCIL_FUNCTIONS()







  
    phi1 = (QSSLIB_REAL *)calloc(g->num_gridpts,sizeof(QSSLIB_REAL));
    COPY_DATA(phi1,phi,g);
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( phi[idx] > 0 )
	  {		
		/* check if neighbors of the current voxel are in the volume */
		if( stencil_size == 6)
                   check_if_6_neighbor_in_volume(idx,nx,ny,nz,nxy,is_in);
                else
                   check_if_26_neighbor_in_volume(idx,nx,ny,nz,nxy,is_in);
				
		for(l = 0; l < stencil_size; l++)
		{
	           idx1 = idx + stencil[l];		     
		     
		    /* if the neighbor is in the volume and UNBURN, add it
		       to the list */  
	           if( ( is_in[l] ) && ( phi[idx1] < 0) )
		   {
			phi1[idx1 ] = -phi[idx1];
		   }
		}
	   }
	 }
	 
   return phi1;	 
}



int  findPhaseCenterOfMass(
   QSSLIB_REAL   *phi,
   QSSLIB_REAL   phase,
   Grid    *g,
   QSSLIB_REAL *comx,
   QSSLIB_REAL *comy,
   QSSLIB_REAL *comz)
{
   int nvox, idx, i, j, k;
   
   *comx = *comy = *comz = 0.0;
   nvox = 0;
   idx = 0;
   for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
     for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
       for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
       {
          if( (phase*phi[idx]) > 0.0 )
	  {
	   nvox++;
           *comx += i;	
           *comy += j;	
           *comz += k;	
          }
       }
       
   *comx = *comx / (QSSLIB_REAL)nvox;
   *comx  = g->x_lo_ghostbox[0] + (*comx)*g->dx[0];
   
   *comy = *comy / (QSSLIB_REAL)nvox;
   *comy  = g->x_lo_ghostbox[1] + (*comy)*g->dx[1];
   
   *comz = *comz / (QSSLIB_REAL)nvox;
   *comz  = g->x_lo_ghostbox[2] + (*comz)*g->dx[2];
   
   return nvox; 
}

/* Functions for connectivity and phase trapping: specific to motion in confined pore spaces */

TrappedPhase *allocateAndInitializeTrappedPhase(
         Options *options,
	 Grid *g,
	 FILE *fp_out,
	 QSSLIB_REAL *phi,
	 QSSLIB_REAL phase,
	 QSSLIB_REAL *mask,
	 QSSLIB_REAL  mask_phase)
{
    TrappedPhase *tp = NULL;
      
    QSSLIB_REAL   zero = 0.0, opposite_phase = -phase;
    int    idx, n1[3], *comp_size;
    char   phase_name[5], fname[256];
    
    int init_from_file = 0;

    tp =  (TrappedPhase *)malloc(sizeof(TrappedPhase));

    /*allocate arrays in the structure */    
    tp->label_trapped =(unsigned char *)calloc(g->num_gridpts,sizeof(unsigned char));
    tp->num_vox = 0;
    tp->num_comp = 0;
    tp->phase = phase;
    tp->num_merge_check = 0;
    tp->num_merge_check_alloc = 256;
    tp->merge_check = (int *)malloc(tp->num_merge_check_alloc*sizeof(int));

    if (phase < 0)
    { 
       tp->connectivity = 6;   //NW fluid
       tp->min_comp_size = 1;
    }
    else
    {
       tp->connectivity = 6;  // W fluid
       tp->min_comp_size = 1; //wetting phase susceptible to more noise, thus 50
    }

    /*initialize arrays*/
    if( (init_from_file) && ((phase < 0) || (phase > 0) ) )
    {  /* initialize trapped phase from file */

	  tp->phi_trapped = readDataArrayQSS(n1,fname);
       
       for(idx=0; idx < g->num_gridpts; idx++)
       {
          if( tp->phi_trapped[idx]*tp->phase > zero ) 
	  {  
	     tp->label_trapped[idx] = TRAP; 
	     tp->num_vox++;
	     
	     //IS THIS THE CORRECT WAY??
	     if(mask_phase*tp->phase < 0)
	     {
	        mask[idx] = tp->phi_trapped[idx];
	     }
	  }
       }
       fprintf(fp_out,"\nTrapped %s phase initialized from %s, %d voxels.",phase_name,
                                                                   fname,tp->num_vox);

       tp->num_comp = findPhaseConnectedComponentsControlVolume(phi,phase,mask,
	                      mask_phase,g,tp->connectivity,&comp_size);
       /* TO-DO: could record component seed voxels here if needed? */			

    }
    else
    {/* initialize with some 'opposite_phase' values, that will be set to 'phase'
	if trapped volumes appear during simulation */
      tp->phi_trapped =  (QSSLIB_REAL *)malloc(g->num_gridpts*sizeof(QSSLIB_REAL));
      for(idx = 0; idx < g->num_gridpts; idx ++)
      {
	tp->phi_trapped[idx] = opposite_phase*fabs(g->dx[0]);
      }

      fprintf(fp_out,"\nTrapped %s phase initialized with 0 voxels.",phase_name);
    }

    return tp;
}

/* SET_DATA_TRAPPED_MARK
   if the examined voxel (that belongs to main phase) already has TRAP mark
   then it recorded for later merge check 
   otherwise it is set to NONTRAP as a candidate for trapping
*/   
#define SET_DATA_TRAPPED_MARK()\
{\
   if( data_trapped[idx] == TRAP )\
   {\
      if( tp->num_merge_check >= tp->num_merge_check_alloc )\
      {\
	tp->num_merge_check_alloc += 256;\
	tp->merge_check = (int*)realloc(tp->merge_check,tp->num_merge_check_alloc*sizeof(int));\
      }\
      \
      tp->merge_check[ tp->num_merge_check ] = idx;\
      (tp->num_merge_check)++;\
   }\
   else\
   {\
      data_trapped[idx] = NONTRAP;\
   }\
}     
     		
int findPhaseConnectedComponentsAndRecordTrappedControlVolume(
     QSSLIB_REAL   *phi,
     QSSLIB_REAL   *control_vol,
     QSSLIB_REAL    control_vol_phase,
     TrappedPhase  *tp,
     Grid          *g)
{
    int   i, j, k, num_comp_trap;
    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
    
    int       idx, idx1, idx2, idx3, l, *new_list_item;
    int       *new_list_item2, list_size, comp_size;
    GList     *burn_list = NULL, *burn_list_ptr;
    QSSLIB_REAL      zero = 0.0, abs_val;
    
    int   touching_exterior;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();
    
    unsigned char *data_trapped = tp->label_trapped;
    QSSLIB_REAL   *phi_trapped =  tp->phi_trapped, phase = tp->phase;  
    QSSLIB_REAL   opposite_phase = -phase;
    int  stencil_size = tp->connectivity;
    int       stencil[stencil_size], is_in[stencil_size];

    tp->num_merge_check = 0;
        
    /* Initialize all appropriate phase voxels to NONTRAP value 
       This is not very elegantly done (too many if's), but
       should be fast
    */                          
    if (control_vol_phase > zero && phase > zero )
    { 
      for(idx = 0; idx < g->num_gridpts; idx++)
      {	 
	 if( (control_vol[idx] > zero) && (phi[idx] > zero) ) 
	 {  
	    SET_DATA_TRAPPED_MARK()
	 }	  
      }
    }  
    else if (control_vol_phase < zero && phase > zero )  
    {
      for(idx = 0; idx < g->num_gridpts; idx++)
      {
	 if( (control_vol[idx] < zero) && (phi[idx] > zero) ) 
	 {  
           SET_DATA_TRAPPED_MARK()
	 }	  
      }
    }  
    else if (control_vol_phase > zero && phase < zero )  
    {
      for(idx = 0; idx < g->num_gridpts; idx++)
      {
	 if( (control_vol[idx] > zero) && (phi[idx] < zero) ) 
	 {  
           SET_DATA_TRAPPED_MARK()
	 }
     }
    }	 
    else //(control_vol_phase < zero && phase < zero )  
    { 
      for(idx = 0; idx < g->num_gridpts; idx++)
      {
	 if( (control_vol[idx] < zero) && (phi[idx] < zero) ) 
	 {  
           SET_DATA_TRAPPED_MARK() 
	 }	  
      }
    }
    
    SET_STENCIL_FUNCTIONS()
    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to 
    *  neighbors of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
         
    /* initialize number of components */	    
    num_comp_trap = 0; 
    
    idx = 0;
    for(k = 0; k < (g->grid_dims_ghostbox)[2]; k++)
      for(j = 0; j < (g->grid_dims_ghostbox)[1]; j++)
        for(i = 0; i < (g->grid_dims_ghostbox)[0]; i++, idx++)
        {
          if( data_trapped[idx] == NONTRAP )
	  {  
             /* initiate list with voxel 'idx' */
	     list_size = 1;
	     data_trapped[idx] = TMPTRAP;
	     new_list_item = (int *)malloc(sizeof(int));
	     *new_list_item = idx;
	     burn_list = g_list_prepend(burn_list,new_list_item);

             touching_exterior = 0;
	     comp_size = 1;
             /* loop over neighboring voxels until there are no more candidates - 
	        this is grass-fire type algorithm */
	     while(list_size)
	     {
	        comp_size++;
	        /* get the first voxel from the list */
		burn_list_ptr = g_list_first(burn_list);
		new_list_item = (burn_list_ptr->data);
		idx1 = *new_list_item;
		
		/* check if neighbors of the first voxel are in the volume */
		check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);                
			
		/* loop over all neighbors of the first voxel*/		
		for(l = 0; l < stencil_size; l++)
		{
	           idx2 = idx1 + stencil[l];		     
		     
		    /* if the neighbor is in the volume and NONTRAP, add it
		     *  to the list and mark with temporary mark TMPTRAP
		     */  
	           if( ( is_in[l] ) && ( data_trapped[idx2] == NONTRAP ) )
		   {
			data_trapped[idx2] = TMPTRAP;
			
			new_list_item2 = (int *)malloc(sizeof(int));
	        	*new_list_item2 = idx2;
			burn_list = g_list_append(burn_list,new_list_item2);
		   }
		   else if( is_in[l] == 0 )
		   {
		        touching_exterior = 1;
		   }   
		}
		
		/* remove the examined voxel from the list */
		burn_list = g_list_remove(burn_list_ptr,new_list_item);
		free(new_list_item);
		list_size = g_list_length(burn_list);
	     }

             /* Replace temporary mark TMPTRAP with the permanent one (TRAP). 
	       
	       If the component is trapped (touching exterior is 0) and
	       its size is big enough to consider it
	       then also record the voxels of trapped phase in
	       phi_trapped and remove them from the original phase phi
	    */
                 	     
	     if((touching_exterior == 0) && (comp_size >= tp->min_comp_size))	
	     {	     
	       /* increase count of trapped connected components */
	       num_comp_trap++;

	       for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	       {
		   if( data_trapped[idx3] == TMPTRAP ) 
		   {
		       (tp->num_vox)++;
		       
		       /*fresh mark for trapped voxels*/
		       data_trapped[idx3] = TRAP;
		       
		       abs_val = fabs(phi[idx3]);
		       /* add this voxel to trapped phase */
		       phi_trapped[idx3] =  phase*abs_val;
		       
		       /* remove the voxel from the invading phase */
		       phi[idx3] = opposite_phase*abs_val;
		       
		       /* this fixes problems for W phase trapping */
		       if (phase*control_vol_phase < 0 )
		       {
		           control_vol[idx3] = phase*abs_val;
		       }
		       		       
		   }
	       }
             }
	     else
	     {
	       for(idx3 = 0; idx3 < g->num_gridpts; idx3++)
	       {
		   if( data_trapped[idx3] == TMPTRAP ) 
		   {
		       data_trapped[idx3] = 0; //remove from consideration
		   }
	       }
	     }  
	  }
       }
     	
    g_list_free(burn_list);
    
    return num_comp_trap;
}


/* \fn checkForTrappedBlobs */
int checkForTrappedBlobs(
     TrappedPhase     *tp,
     QSS_DataArrays   *p,
     Grid             *g,
     QSSLIB_REAL      mask_sign,
     FILE             *fp_out)
 {
     int num_comp_trap, idx;
     char   phase_name[5];
     QSSLIB_REAL zero=0.0;
    
     if(tp->phase < 0) sprintf(phase_name,"NW");
     else              sprintf(phase_name,"W");

      /* Check if any group of tp->phase voxels in p->phi became trapped
       (i.e. isolated/disconnected from the rest) */
     num_comp_trap = findPhaseConnectedComponentsAndRecordTrappedControlVolume(
	             p->phi,p->mask,mask_sign,tp,g);
     if(num_comp_trap)
     {
        fprintf(fp_out,"\nNumber of newly trapped %s phase components %d,",
	                                             phase_name,num_comp_trap);
	fprintf(fp_out," total trapped %d voxels\n",tp->num_vox);
	//writeDataArray(p->phi,g,"phi",GZIP);
	//writeDataArray(tp->phi_trapped,g,"phi_trap",GZIP);					        
     }
     /* Any already trapped voxel touched by the moving fluid interface have been tagged
       in findPhaseConnectedComponentsAndRecordTrappedControlVolume()
	Merge relevant trapped blobs (connected component) into the main phase. 
     */
     tp->num_merged_vox = 0;
     
     if(tp->num_merge_check)
     { 
       mergeTrappedComponentIntoMainPhase(p->phi,tp,g);
     }  	
     if (tp->num_merged_vox)
     {
        fprintf(fp_out,"\nNumber of trapped %s phase voxels merged into main phase %d\n",
	                                               phase_name,tp->num_merged_vox);
     
	//writeDataArray(p->phi,g,"phi_merge",GZIP);
	//writeDataArray(tp->phi_trapped,g,"phi_merge_trap",GZIP);	
     }						       
     return 	num_comp_trap;					  
 } 
	

void mergeTrappedComponentIntoMainPhase(
     QSSLIB_REAL   *phi,
     TrappedPhase  *tp,
     Grid          *g)
{    
    int   nx = (g->grid_dims_ghostbox)[0];
    int   ny = (g->grid_dims_ghostbox)[1];
    int   nz = (g->grid_dims_ghostbox)[2];
    int   nxy = nx*ny;
 
    int   idx1, idx2, i, l, *new_list_item;
    int   *new_list_item2, list_size;
    GList *burn_list = NULL, *burn_list_ptr;
    QSSLIB_REAL  abs_val;
    
    unsigned char *data_trapped = tp->label_trapped;
    QSSLIB_REAL   *phi_trapped =  tp->phi_trapped, phase = tp->phase;
    QSSLIB_REAL   opposite_phase = -phase;
    int  stencil_size = tp->connectivity;
    int   stencil[stencil_size], is_in[stencil_size];
    int  search_seed_idx;
    
    /* function pointers */
    void         (*set_stencil_func)();
    void         (*check_if_neighbor_in_volume_func)();   
    
  
    SET_STENCIL_FUNCTIONS()
    
    /* Stencil stores shifts (in a packed 2d or 3d array) corresponding to neighbors 
    *  of a voxel in 4/8 (for 2d) or 6/26 (for 3d) connectivity sense 
    */
    set_stencil_func((g->grid_dims_ghostbox),stencil);
    
    tp->num_merged_vox = 0;
    
    for(i=0; i < tp->num_merge_check; i++)
    {
      /* tp->merge_check stores voxels that have been tagged as (previously) trapped phase
         that has been touched by the moving interface
	 Start a grass-fire algorithm from such voxels, merge all connected
	 ones to the main phase while removing them from the trapped phase.
      */
      search_seed_idx = tp->merge_check[i];
       
      if( data_trapped[ search_seed_idx ] == TRAP )
      {     
       /* initiate the list with voxel 'search_seed_idx' */
	list_size = 1;
	new_list_item = (int *)malloc(sizeof(int));
	*new_list_item = search_seed_idx;
	burn_list = g_list_prepend(burn_list,new_list_item);
	data_trapped[search_seed_idx] = NONTRAP;
        tp->num_merged_vox++;

	/* loop until you don't find all the voxels connected to the 
	   initial voxel on the list */
	while(list_size)
	{
	   /* get the first voxel from the list */
	   burn_list_ptr = g_list_first(burn_list);
	   new_list_item = (burn_list_ptr->data);
	   idx1 = *new_list_item;
  
	   /* check if neighbors of the first voxel are in the volume */
	   check_if_neighbor_in_volume_func(idx1,nx,ny,nz,nxy,is_in);                

	   /* loop over all neighbors of the first voxel */		
	   for(l = 0; l < stencil_size; l++)
	   {
	      idx2 = idx1 + stencil[l];		     

	       /* if the neighbor is in the volume and TRAP, add it
		*  to the list
		* also remove it from the trapped phase and add to the main phase
		*/  
	      if( ( is_in[l] ) && ( data_trapped[idx2] & TRAP ) )
	      {
		   tp->num_merged_vox++;
		   new_list_item2 = (int *)malloc(sizeof(int));
		   *new_list_item2 = idx2;
		   burn_list = g_list_append(burn_list,new_list_item2);

		   /* remove this voxel from trapped phase */
		   data_trapped[idx2] = NONTRAP;
		   abs_val = fabs(phi_trapped[idx2]);
		   phi_trapped[idx2] = opposite_phase*abs_val;

		   /* add this voxel to main phase */
		   /* TO-DO - not sure this works for wetting phase.... */
        	   //phi[idx2] = (phi_trapped[idx2] < phi[idx2] ) ? phi_trapped[idx2] : phi[idx2];
		   phi[idx2] = phase*abs_val;
	      }
	   }

	   /* remove the examined voxel from the list */
	   burn_list = g_list_remove(burn_list,new_list_item);
	   free(new_list_item);
	   list_size = g_list_length(burn_list);
	}
      }
    }	
    
    tp->num_vox -= tp->num_merged_vox;
 
    tp->num_merge_check = 0;
    
    g_list_free(burn_list);
}

/*! \fn destroyTrappedPhase() frees memory for TrappedPhase structure */
void destroyTrappedPhase(TrappedPhase *tp)
{

    free(tp->phi_trapped);
    free(tp->label_trapped);
    free(tp->merge_check);  
    free(tp);
}

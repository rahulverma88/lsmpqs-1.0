/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file qss_data_arrays.c

    Function definitions for creating, reading, manipulating, and writing data arrays.
             
*/
#include <stdio.h>
#include <stdlib.h>

#include "qss_data_arrays.h"

#define DSZ  sizeof(QSSLIB_REAL)
#define ISZ  sizeof(int)
#define UCSZ sizeof(unsigned char)
#define TCSZ sizeof(size_t)

#define QSSLIB_SERIAL_dummy_pointer        ((QSSLIB_REAL*)(-1))
#define QSSLIB_SERIAL_dummy_pointer_int    ((int*)(-1))
#define QSSLIB_SERIAL_dummy_pointer_uchar  ((unsigned char*)(-1))
#define QSSLIB_SERIAL_dummy_pointer_sizet  ((size_t*)(-1))

/*
    Allocates dummy pointers for all members of the QSS Data arrays structure.
*/
QSS_DataArrays *allocateQSSDataArrays(void)
{
  QSS_DataArrays *qss_data_arrays;
  int i;
  
  qss_data_arrays  = (QSS_DataArrays *)malloc(sizeof(QSS_DataArrays));
  
  qss_data_arrays->phi = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->phi_next = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->phi_extra = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->phi_prev = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->mask = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->lse_rhs = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->normal_velocity = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->curvature_coeff = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->external_velocity_x = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->external_velocity_y = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->external_velocity_z = QSSLIB_SERIAL_dummy_pointer;
  
  qss_data_arrays->mask_x = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->mask_y = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->mask_z = QSSLIB_SERIAL_dummy_pointer;
  
  qss_data_arrays->scratch1 = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->scratch2 = QSSLIB_SERIAL_dummy_pointer;

  qss_data_arrays->connectivity = QSSLIB_SERIAL_dummy_pointer_int;
  qss_data_arrays->phi_bin = QSSLIB_SERIAL_dummy_pointer_uchar;

  qss_data_arrays->mask_w = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->mask_nw = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->mask_disconn_init = QSSLIB_SERIAL_dummy_pointer;

  qss_data_arrays->theta = QSSLIB_SERIAL_dummy_pointer;
  qss_data_arrays->overlap = QSSLIB_SERIAL_dummy_pointer;

  return  qss_data_arrays;
}

/*
    Destroys memory for QSS Data arrays.
*/
void destroyQSSDataArrays(QSS_DataArrays *qss_data_arrays)
{
  if (qss_data_arrays) {
  
  free(qss_data_arrays->phi);
  
  free(qss_data_arrays->phi_next);
  
  free(qss_data_arrays->phi_extra);      
  free(qss_data_arrays->phi_prev);  
  free(qss_data_arrays->mask); 
  
  free(qss_data_arrays->lse_rhs);
  
  free(qss_data_arrays->normal_velocity);
  free(qss_data_arrays->curvature_coeff);
  free(qss_data_arrays->external_velocity_x);
  free(qss_data_arrays->external_velocity_y);
  free(qss_data_arrays->external_velocity_z);
 
  free(qss_data_arrays->mask_x);
  free(qss_data_arrays->mask_y);
  free(qss_data_arrays->mask_z);

  free(qss_data_arrays->scratch1);
  free(qss_data_arrays->scratch2);
  
  free(qss_data_arrays->connectivity);
  free(qss_data_arrays->phi_bin);
  free(qss_data_arrays->mask_w);
  free(qss_data_arrays->mask_nw);
  free(qss_data_arrays->mask_disconn_init);

  free(qss_data_arrays->theta);
  free(qss_data_arrays->overlap);


  free(qss_data_arrays);
  }
}

/*
    Allocates memory for QSS Data arrays (as opposed to simply pointing to a dummy pointer
*/
void  allocateMemoryForQSSDataArrays(
  QSS_DataArrays *qss_data_arrays,
  Grid *grid)
{
 /* Only arrays that are equal to QSSLIB_SERIAL_dummy_pointer will get memory allocated.
  *   If memory allocation is to be avoided, set the pointer to NULL,
  *   Non-NULL pointers different from QSSLIB_SERIAL_dummy_pointer are assumed allocated
  *   elsewhere and that will not be overridden.
  */
        
  if( qss_data_arrays->phi == QSSLIB_SERIAL_dummy_pointer )
    qss_data_arrays->phi = (QSSLIB_REAL*) calloc(grid->num_gridpts,DSZ);
    
  if( qss_data_arrays->phi_next  == QSSLIB_SERIAL_dummy_pointer ) 
    qss_data_arrays->phi_next = (QSSLIB_REAL*) calloc(grid->num_gridpts,DSZ);
  
  if( qss_data_arrays->phi_extra  == QSSLIB_SERIAL_dummy_pointer )
    qss_data_arrays->phi_extra = (QSSLIB_REAL*) calloc(grid->num_gridpts,DSZ);  
  
  if( qss_data_arrays->phi_prev  == QSSLIB_SERIAL_dummy_pointer ) 
    qss_data_arrays->phi_prev = (QSSLIB_REAL*) calloc(grid->num_gridpts,DSZ);
   
  if( qss_data_arrays->mask  == QSSLIB_SERIAL_dummy_pointer )
    qss_data_arrays->mask = (QSSLIB_REAL*) calloc(grid->num_gridpts,DSZ);  
    
  if( qss_data_arrays->lse_rhs  == QSSLIB_SERIAL_dummy_pointer )
    qss_data_arrays->lse_rhs = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if( qss_data_arrays->normal_velocity  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->normal_velocity = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
 
   if( qss_data_arrays->curvature_coeff  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->curvature_coeff = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
  
  if( qss_data_arrays->external_velocity_x  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->external_velocity_x = 
      (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if( qss_data_arrays->external_velocity_y  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->external_velocity_y = 
      (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
    
  if(grid->num_dims == 3)
  {
    if( qss_data_arrays->external_velocity_z  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->external_velocity_z = 
        (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);	
  }
  else qss_data_arrays->external_velocity_z = (QSSLIB_REAL *)NULL;
  
  if( qss_data_arrays->mask_x  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->mask_x = 
      (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
   
  if( qss_data_arrays->mask_y  == QSSLIB_SERIAL_dummy_pointer )  
    qss_data_arrays->mask_y = 
      (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
    
  if(grid->num_dims == 3)
  {
    if( qss_data_arrays->mask_z  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->mask_z = 
        (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);	
  }
  else qss_data_arrays->mask_z = (QSSLIB_REAL *)NULL;  
  

  if( qss_data_arrays->scratch1  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->scratch1 = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);

  if( qss_data_arrays->scratch2  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->scratch2 = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
      
  if( qss_data_arrays->connectivity  == QSSLIB_SERIAL_dummy_pointer_int )  
      qss_data_arrays->connectivity = (int*) malloc(grid->num_gridpts*ISZ);
      
  if( qss_data_arrays->phi_bin  == QSSLIB_SERIAL_dummy_pointer_uchar )  
      qss_data_arrays->phi_bin = (unsigned char*) malloc(grid->num_gridpts*UCSZ);
      
  if( qss_data_arrays->theta  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->theta = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);

  if( qss_data_arrays->overlap  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->overlap = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);

  if( qss_data_arrays->mask_w  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->mask_w = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);

  if( qss_data_arrays->mask_nw  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->mask_nw = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
      
  if( qss_data_arrays->mask_disconn_init  == QSSLIB_SERIAL_dummy_pointer )  
      qss_data_arrays->mask_disconn_init = (QSSLIB_REAL*) malloc(grid->num_gridpts*DSZ);
}     

/*
    Writes QSS Data Arrays into binary files.
*/
void writeDataArrayQSS(QSSLIB_REAL *data, Grid *grid, char *file_name,int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write grid dimensions */
   fwrite(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 

   /* write data array */
   fwrite(data, DSZ, grid->num_gridpts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}

/* same as writeDataArray, but without dimensions, so file is readable by Paraview */
void writeDataArrayRaw(QSSLIB_REAL *data, Grid *grid, char *file_name,int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write data array */
   fwrite(data, DSZ, grid->num_gridpts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}

/* 
    Writes data array which is an unsigned binary character array, instead of the usual
    floating point array of a level set 
*/
void writeDataArrayUchar(unsigned char *data, Grid *grid, char *file_name,int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write grid dimensions */
   fwrite(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 

   /* write data array */
   fwrite(data, sizeof(unsigned char), grid->num_gridpts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}

/* 
    Writes data array which is an unsigned binary character array, instead of the usual
    floating point array of a level set. Additionally, it only writes out the inner fill box,
    ignoring the ghost cells.
*/
void writeDataArrayUcharFB(unsigned char *data_fb, Grid *g, char *file_name, int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   int num_pts = g->grid_dims[0]*g->grid_dims[1]*g->grid_dims[2];
   fwrite(data_fb, sizeof(unsigned char), num_pts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);

}

/* 
    Writes data array which is an int array, instead of the usual
    floating point array of a level set 
*/
void writeDataArrayInt(int *data, Grid *grid, char *file_name,int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write grid dimensions */
   fwrite(grid->grid_dims_ghostbox, sizeof(int), 3, fp); 

   /* write data array */
   fwrite(data, sizeof(int), grid->num_gridpts, fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}

/* 
    Reads data array from a binary file as an int array, instead of the usual
    floating point array of a level set 
*/
int *readDataArrayInt(int *grid_dims_ghostbox,char *file_name)
{
   FILE    *fp;
   int     zip_status;
   int     num_gridpts;
   int    *data = NULL;
   char    *file_base;
   
   checkUnzipFile(file_name,&zip_status,&file_base);
   
   fp = fopen(file_base,"r");

   if( fp != NULL)
   {
     /* read grid dimensions */
     fread(grid_dims_ghostbox, sizeof(int), 3, fp); 
  
     /* allocate memory for data array */ 
     num_gridpts = grid_dims_ghostbox[0] * grid_dims_ghostbox[1]
               * grid_dims_ghostbox[2];
     data = (int *) malloc(num_gridpts*ISZ);

     /* read data array */ 
     fread(data, ISZ, num_gridpts, fp);

     fclose(fp);
     
     zipFile(file_base,zip_status);
   }
   else
   {
      printf("\nCould not open file %s",file_name);
   }
   free(file_base);
   return data;
}

/* 
    Reads data array from a binary file.
*/
QSSLIB_REAL *readDataArrayQSS(int *grid_dims_ghostbox,char *file_name)
{
   FILE    *fp;
   int     zip_status;
   int     num_gridpts;
   QSSLIB_REAL    *data = NULL;
   char    *file_base;
   
   checkUnzipFile(file_name,&zip_status,&file_base);
   
   fp = fopen(file_base,"r");

   if( fp != NULL)
   {
     /* read grid dimensions */
     fread(grid_dims_ghostbox, sizeof(int), 3, fp); 
  
     /* allocate memory for data array */ 
     num_gridpts = grid_dims_ghostbox[0] * grid_dims_ghostbox[1]
               * grid_dims_ghostbox[2];
     data = (QSSLIB_REAL *) malloc(num_gridpts*DSZ);

     /* read data array */ 
     fread(data, DSZ, num_gridpts, fp);

     fclose(fp);
     
     zipFile(file_base,zip_status);
   }
   else
   {
      printf("\nCould not open file %s",file_name);
   }
   free(file_base);
   return data;
}



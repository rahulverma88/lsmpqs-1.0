/*
 * File:        qss_data_arrays.h
 * Function headers for creating and manipulating QSS data arrays
 */

#ifndef INCLUDED_QSS_DATA_ARRAYS_H
#define INCLUDED_QSS_DATA_ARRAYS_H

#include "QSSLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif



#include "qss_grid.h"
#include "qss_file.h"

/*!
 * Structure 'QSS_DataArrays' stores pointers for all arrays needed in a
 * typical QSS computation.
 */
typedef struct _QSS_DataArrays
{
  /* level set function at different time integration steps*/
  QSSLIB_REAL  *phi, *phi_next, *phi_extra;
  
  /* extra storage in reinitialization */
  //QSSLIB_REAL  *phi0;      

  /* extra storage for previous step functions */
  QSSLIB_REAL  *phi_prev;  
   
  /* mask is the level set function that defines restricted domains */
  QSSLIB_REAL  *mask;

   /* LS Equation right hand side */
  QSSLIB_REAL  *lse_rhs; 
  
  /* normal velocity */
  QSSLIB_REAL *normal_velocity; 
  
  /* curvature coefficient */
  QSSLIB_REAL *curvature_coeff; 
  
   /* external velocity field */
  QSSLIB_REAL *external_velocity_x;
  QSSLIB_REAL *external_velocity_y;
  QSSLIB_REAL *external_velocity_z;
   
 /* Arrays storing gradients of the mask */
  QSSLIB_REAL *mask_x, *mask_y, *mask_z;
  
  /* Scratch arrays for reinitialization */
  QSSLIB_REAL *scratch1, *scratch2;
  
  /* Connectivitiy array */
  int *connectivity;
  
  /* Binary array - 0s and 1s */
  unsigned char *phi_bin;
  
  /* Disconnected component masks */
  QSSLIB_REAL *mask_w, *mask_nw, *mask_disconn_init;
  
  /* Contact angle array */
  QSSLIB_REAL *theta;
  
  /*Overlap array */
  QSSLIB_REAL *overlap;

}  QSS_DataArrays;

/*!
 * allocateQSSDataArrays() allocates a QSS_DataArrays data structure 
 * and initializes all of its data pointers to a non-NULL dummy pointer.
 *  
 * Arguments:    none
 *
 * Return value: pointer to QSS_DataArrays structure 
 *
 */
QSS_DataArrays *allocateQSSDataArrays(void);

/*!
 * allocateMemoryForQSSDataArrays() allocates memory for the data 
 * arrays contained within the QSS_DataArrays structure.
 *    
 * Arguments:
 *  - lsm_arrays(in):  pointer to QSS_DataArrays structure
 *  - grid(in):        pointer to Grid 
 *    
 * Return value:       none
 *
 * NOTES: 
 * - The memory for the QSS_DataArrays structure MUST be allocated
 *   (e.g. using the allocateLSMDataArrays() function) before 
 *   allocateMemoryForQSSDataArrays() is called.
 *
 * - Memory is allocated only for array pointers in lsm_data_arrays that 
 *   are not already associated with allocated or that are set to NULL.
 *   If memory has already been allocated for a particular data array or
 *  the data pointer is set to NULL, it will not be reallocated.
 *
 */
void allocateMemoryForQSSDataArrays(
  QSS_DataArrays *qss_data_arrays,
  Grid *grid);

/*!
 * destroyQSSDataArrays() frees ALL memory allocated for the data 
 * arrays contained within the QSS_DataArrays structure as well as the structure
 * itself.
 *   
 * Arguments:
 *  - qss_data_arrays(in):  pointer to QSS_DataArrays 
 *   
 * Return value:            none
 *   
 */
void destroyQSSDataArrays(QSS_DataArrays *qss_data_arrays);

/*!
 * writeDataArrayQSS() writes the specified data array out to a binary file.
 *
 * The data is output in the following order:
 * -# grid dimensions 
 * -# values of data array at all grid points.
 *
 * Arguments:
 *  - data (in):       data array to be output to file
 *  - grid (in):       pointer to Grid 
 *  - file_name (in):  name of output file
 *  - zip_status(in):  integer indicating compression of the file 
 *                     (NO_ZIP,GZIP,BZIP2) 
 *   
 * Return value:       none
 *   
 * NOTES: 
 * - writeDataArrayQSS() is used for 2d and 3d data arrays.  For 2d
 *   data, the third grid dimension MUST be set to 1 (which is 
 *   the default behavior of the createGrid() function).
 *
 * - If a file with the specified file_name already exists, it is
 *   overwritten.
 *
 */   
void writeDataArrayQSS(QSSLIB_REAL *data, Grid *grid, char *file_name,
                    int zip_status);
void writeDataArrayRaw(QSSLIB_REAL *data, Grid *grid, char *file_name,int zip_status);

void writeDataArrayUchar(unsigned char *data, Grid *grid, char *file_name,
                    int zip_status);

void writeDataArrayInt(int *data, Grid *grid, char *file_name,
                    int zip_status);
/*!
 * readDataArrayQSS() loads the data from a binary file into a QSSLIB_REAL
 * array and returns it to the user.  
 *   
 * Arguments:
 *  - grid_dims (out):  dimensions of grid (read from file)
 *  - file_name (in):   name of input file 
 *   
 * Return value:        pointer to data array loaded from file
 *   
 * NOTES: 
 * - readDataArray() dynamically allocates memory for the data array 
 *   that is returned.
 *
 * - The memory for grid_dims is assumed to be allocated by the user.
 *
 * - readDataArray() is used for 2d and 3d data arrays.  For 2d
 *   data, the third grid dimension is set to 1.
 *
 * - Function recognizes if the file name contains .gz or .bz2 extention
 *   and uncompresses the file accordingly.
 */   
QSSLIB_REAL *readDataArrayQSS(int *grid_dims, char *file_name);

void writeDataArrayUcharFB(unsigned char *, Grid *, char *, int);

#ifdef __cplusplus
}
#endif

#endif


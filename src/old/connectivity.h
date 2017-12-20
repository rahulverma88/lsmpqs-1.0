/******************************************************************************
 *
 *   Author:   Masa Prodanovic
 *   Copyright (c) 2009, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
#ifndef INCLUDED_CONNECTIVITY_H
#define INCLUDED_CONNECTIVITY_H

#include "QSSLIB_config.h"
#include "qss_data_arrays.h"
#include "qss_grid.h"
#include "qss_options.h"
#include "qss_macros.h"

/*! \struct TrappedPhase
   structure for recording disconnected(trapped) fluid components(blobs) 
*/
typedef struct {
    /*! level set function describing trapped phase */
    QSSLIB_REAL *phi_trapped;
    
    /*! unsigned array of labels for trapped phase  */        
    unsigned char *label_trapped;
       
    /*! number of voxels in trapped phase */
    int  num_vox;
    
    /*! connectivity type for trapped phase voxels */
    int connectivity;
    
    /*! number of connected components (assuming above connectivity */
    int num_comp;
    
    /*! minimum component size (to be considered at all) */
    int min_comp_size;
    
    /*! phase of the main level set function 
        -1 if the phase is described by phi < 0; 1 otherwise */
    QSSLIB_REAL phase;
    
    /*! scratch space for merge/rechecking candidates, used internally */
    int num_merge_check, num_merge_check_alloc, *merge_check;
    
    /*! scratch space for number of merged voxels, used  internally*/
    int num_merged_vox;
                  
} TrappedPhase;


#define INT_UNBURN -1
#define INT_BURN   -2

#define UNBURN  1
#define TMPBURN 2
#define BURN    4

#define NONTRAP     1
#define TMPTRAP     2
#define TRAP        4


/*! Various utilities for finding voxels on volume boundary */

/*! \fn 
*   find_max_conn_comp2d() calculates number of pixels in the maximal/minimal connected 
*   component of pore space (phi < 0) in specified row(column) 'pos'
*   
*   Arguments
*     dir(in)   - desired direction ('x' or 'y')
*     pos (in)  - integer position of the row/column
*     phi(int)  - level set function data
*     g(in)     - pointer to Grid data
*     max_conn_sz(out)  - (integer) pixel count of the max. connected component
*     min_conn_sz(out)  - (integer) pixel count of the min. connected component
*/
void find_max_conn_comp2d(
     char    dir,
     int     pos,
     QSSLIB_REAL   *phi,
     Grid   *g,
     int    *max_conn_sz,
     int    *min_conn_sz);
     
/*! \fn  
*    find_max_conn_comp3d() finds number of voxels in min/max connected component of the 
*   set where phi < 0 in the specified plane (slice through data)
*   Arguments
*     dir(in)   - desired direction ('x', 'y' or 'z')
*     pos (in)  - integer position of the row/column/slice
*     phi(int)  - level set function data
*     g(in)     - pointer to Grid data
*     max_conn_sz(out)  - (integer) pixel count of the max. connected component
*     min_conn_sz(out)  - (integer) pixel count of the min. connected component
*/
void find_max_conn_comp3d(
     char    dim,
     int     pos,
     QSSLIB_REAL *phi,
     Grid   *g,
     int    *max_conn_sz,
     int    *min_conn_sz);    

/*! \fn 
*   find_opp_side2d() finds voxel indices (and returns them in the array 'vox_ind')
*   that satisfy (mask < 0) and (phi > 0) are are located on the volume edges. 
*   In terms of LSMPQS convention for fluids in the pore space this is
*   equivalent to finding pore space voxels (mask < 0) that are  
*   (wetting) defending fluid (phi > 0).
*   Dimension of data assumed to be 2.
*
*   Arguments:
*      phi(in) - level set function
*     mask(in) - masking level set function (control volume)
*     g(in)    - Grid structure
*     pvox_ind - pointer to array of voxel of the opposite phase, memory allocated
*                within the function
*
*   Returns:
*     number of voxels in array pvox_ind
*/   
int find_opp_side2d(
    QSSLIB_REAL    *phi,
    QSSLIB_REAL    *mask,
    Grid           *g,
    int            **pvox_ind);
    
/*! \fn find_opp_side3d()
   Finds voxel indices (and returns them in the array 'vox_ind')
   that are pore space voxels (mask < 0) on the volume boundary that are not
   initialized as invading fluid (i.e. phi > 0)
   dim == 3
*/   
int find_opp_side3d(
    QSSLIB_REAL *phi,
    QSSLIB_REAL *mask,
    Grid *g,
    int **pvox_ind);
      
/* ! \fn
*   find_opp_side3d_x() finds voxel indices (and returns them in the array 'vox_ind')
*   that satisfy (mask < 0) and (phi > 0). In terms of fluid in the pore space this is
*   equivalent to pore space voxels (mask < 0) that are  (wetting) defending fluid (phi > 0).
*   Dimension of data assumed to be 3.
*   ONLY x-direction volume boundary planes are checked.
*
*/ 
int  find_opp_side3d_x(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int **pvox_ind);   

/*! \fn find_yz_plane_voxels()
   Finds opposite fluid voxels (like above) in a specific yz-plane 
   Used for debugging
*/
int  find_yz_plane_voxels(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int position,
     int **pvox_ind);

/*! \fn find_multiple_yz_plane_voxels()
   Finds opposite fluid voxels (like above) in the a number of consecutive yz-planes 
   Used for debugging
*/
int  find_multiple_yz_plane_voxels(
     QSSLIB_REAL *phi,
     QSSLIB_REAL *mask,
     Grid *g,
     int position_start,
     int position_end,
     int **pvox_ind);

/*! Functions setting stencils for finding neighbors in 2d/3d packed arrays */

/*! \fn set_*_stencil() functions return stencil of shifts to reach nearest neighbors
   for voxels in a packed array.
   The number refers to connectivity - 4,8 in 2D, 6,26 in 3D.
   
   Arguments 
   n(in)        - n[0],n[1],n[2] is number of voxels in each dimension
   sten(in/out) - stencil of shifts in packed array
*/
void set_4_stencil(
	int	*n,
	int	*sten);
	
void set_8_stencil(
	int	*n,
	int	*sten);

void set_6_stencil(
	int	*n,
	int	*sten);

void set_26_stencil(
	int	*n,
	int	*sten);

/*! \fn check_if_*_neighbor_in_volume() functions check whether a particular neighbor
  (as defined by the stencil shift) belongs to the volume. If position in digitized
  space is (i,j,k), then the neighbor (i,j,k+1) will not be in the volume for the 
  largest value of k, etc...
  
   The number refers to connectivity - 4,8 in 2D, 6,26 in 3D.
   
   Arguments 
   pos(in) - packed array position whose neighbors need to be checked
   n*(in)  - array dimensions in x,y,z, directions
             (for 2d, enter nz = 1)
   nxy(in) - nx*ny, provided to function for speed
   is_in_volume(out) - array whose elements correspond to those in relevant stencil
                       1 for the corresponding neighbor in volume, 0 for not	     
*/
void check_if_4_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume);
	
void check_if_8_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume);
	
void check_if_6_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume);
	
void check_if_26_neighbor_in_volume(
	int		pos,
	int		nx,
	int		ny,
	int		nz,
	int		nxy,
	int		*is_in_volume);

/*! \def SET_STENCIL_FUNCTIONS()
   Macro that shortens writing and makes it easy to create functions that work
   for all types of stencil/connectivity/dimension.
   Assumes 'stencil_size' exists and sets function pointers
   'set_stencil_func' and 'check_if_neighbor_in_volume_func' to appropriate 
   functions.

   See findPhaseConnectedComponents() for example use.
*/ 
#define SET_STENCIL_FUNCTIONS()\
  switch( stencil_size )\
  {\
     case  4: set_stencil_func = set_4_stencil;\
	      check_if_neighbor_in_volume_func = check_if_4_neighbor_in_volume;\
	      break;\
     case  8: set_stencil_func = set_8_stencil;\
	      check_if_neighbor_in_volume_func = check_if_8_neighbor_in_volume;\
	      break;\
     case  6: set_stencil_func = set_6_stencil;\
	      check_if_neighbor_in_volume_func = check_if_6_neighbor_in_volume;\
	      break;\
     case 26: set_stencil_func = set_26_stencil;\
	      check_if_neighbor_in_volume_func = check_if_26_neighbor_in_volume;\
	      break;\
     default: set_stencil_func = set_4_stencil;\
	      check_if_neighbor_in_volume_func = check_if_4_neighbor_in_volume;\
	      break;\
  }
    	
/*! Functions for finding and labeling connected components of sets 
    of voxels defined by level set functions
*/
						
/*! \fn  findPhaseConnectedComponents()
*   Find the number of connected components of the set is described by either
*   positive or negative part of level set function.
*   Arguments
*   phi(in)           - level set function
*   phase(in)         - 1.0 or -1.0 depending on whether positive or negative level
*                       set function values should be considered
*   g(in)             - pointer to Grid structure
*   stencil_size(in)  - voxel connectivity type (6 or 26), i.e. number of neighborhood
*                       voxels to be examined for each voxel   
*   pcomponent_size(out)   - array with sizes of the connected components 
*
*   Returns
*     number of connected components    
*/
int findPhaseConnectedComponents(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size);
     

/*! \fn   findPhaseConnectedComponentsControlVolume()
*   Find the number of connected components of the set that is described by either
*   positive or negative phase of a level set function.
*   Arguments
*   phi(in)           - level set function
*   phase(in)         - 1.0 or -1.0 depending on whether positive or negative level
*                       set function values should be considered
*   control_vol(in)   - control volume (mask) level set function
*   control_vol_phase(in) - phase of control volume to consider - only points in this 
*                       phase will be taken into account
*   g(in)             - pointer to Grid structure
*   stencil_size(in)  - voxel connectivity type (6 or 26), i.e. number of neighborhood
*                       voxels to be examined for each voxel   
*   pcomponent_size(out)   - array with sizes of all connected components 
*
*   Returns
*     number of connected components    
*/ 
int findPhaseConnectedComponentsControlVolume(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     QSSLIB_REAL *control_vol,
     QSSLIB_REAL control_vol_phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size);

/*! \fn  processPhaseConnectedComponents()
*   Finds connected components of the set is described by either
*   positive or negative part of level set function. The components that
*   have less than min_size voxels are set to the opposite phase.
*   The largest connected component is set to the opposite phase as well.
*  (This is done with viewing residual wetting phase in mind, otherwise watch out - 
*   data is lost!).
*
*   Arguments
*   phi(in) - level set function
*   phase(in) - 1.0 or -1.0 depending on whether to consider positive level
*               set function values or negative
*   g(in)   - grid structure
*   stencil_size(in) - voxel connectivity type (6 or 26), i.e. number of neighborhood
*                 voxels to be examined for each voxel   
*   pcomponent_size(out) - array with sizes of the maximal connected component 
*
*   Returns
*     number of connected components    
*/
int processPhaseConnectedComponents(
     QSSLIB_REAL *phi,
     QSSLIB_REAL phase,
     Grid   *g,
     int     stencil_size,
     int     min_size);
     
     
/*! \fn  labelPhaseConnectedComponents()
*   Find the number of connected components of the set is described by positive
*   or negative part of level set function and label them in a separate integer array.
*
*   Arguments
*   phi(in) - level set function
*   phase(in) - 1.0 or -1.0 depending on whether to consider positive level
*               set function values or negative
*   g(in)   - grid structure
*   stencil_size(in) - voxel connectivity type (6 or 26), i.e. number of neighborhood
*                 voxels to be examined for each voxel   
*   pcomponent_size(out) - array with sizes of the connected components  
*   pcomponent_size(out) - array with sizes of the connected components    
*   plabel_phi(out) - pointer to integer array whose positive integer values are
*                labels of connected components of the phase, or 0 if voxel
*                does not belong to the phase, component numbering starts at min_label
*   Returns
*     number of connected components    
*/
int labelPhaseConnectedComponents(
     QSSLIB_REAL  *phi,
     QSSLIB_REAL  phase,
     Grid   *g,
     int     stencil_size,
     int     **pcomponent_size,
     int     **pcomp_start_vox,
     int     **plabel_phi,
     int     min_label);
 
 /*! \fn  labelPhaseConnectedComponentsControlVolume()
*   Find the number of connected components of the set is described by positive
*   or negative part of level set function and label them in a separate integer array.
*
*   Arguments
*   phi(in) - level set function
*   phase(in) - 1.0 or -1.0 depending on whether to consider positive level
*               set function values or negative
*   control_vol(in)   - control volume (mask) level set function
*   control_vol_phase(in) - phase of control volume to consider - only points in this phase will be taken into
*                       account
*   g(in)   - grid structure
*   stencil_size(in) - voxel connectivity type (6 or 26), i.e. number of neighborhood
*                 voxels to be examined for each voxel   
*   pcomponent_size(out) - array with sizes of the connected components  
*   pcomponent_size(out) - array with sizes of the connected components    
*   plabel_phi(out) - pointer to integer array whose positive integer values are
*                labels of connected components of the phase, or 0 if voxel
*                does not belong to the phase, component numbering starts at min_label
*   Returns
*     number of connected components    
*/
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
     int     min_label);
     
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
    char    *basename);

  
/*! \fn numberLabeledComponentsInPhaseNeighborhood() checks  neighborhood of the level set function
*  phase and returns number of different labeled components in its neighborhood.
*
*  Arguments
*    phi - level set function array
*    phase - phase of the level set function to examine (positive or negative)
*    labels_data - array with integer labels in [1,max_num_labels]
*    max_num_labels - max number of integer labels
*    plabel_found - pointer to newly allocated array of max_num_labels elements
*                   *plabel_found[i] = 1 if label(component) i+1 is identified 
*                   in the neighborhood of the phase
*    g - Grid structure
*    stencil_size - stencil size, i.e. pixel/voxel connectivity (4,8,6 or 26) 
*    
*/
  
int numberLabeledComponentsInPhaseNeighborhood(
    QSSLIB_REAL  *phi,
    QSSLIB_REAL   phase,
    int    *labels_data,
    int    max_num_labels,
    unsigned char **plabel_found,
    Grid   *g,
    int    stencil_size);

  
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
    int    stencil_size);


/*! \fn findPhaseCenterOfMass()
*/
int  findPhaseCenterOfMass(
   QSSLIB_REAL   *phi,
   QSSLIB_REAL   phase,
   Grid    *g,
   QSSLIB_REAL *comx,
   QSSLIB_REAL *comy,
   QSSLIB_REAL *comz);
   

/*! Functions for connectivity and phase trapping: specific to motion in confined pore spaces */

/*! \fn allocateAndInitializeTrappedPhase initializes and allocates memory for
   trapped NW or W phase.
      
   Arguments:
   options(in) - Options structure
   g(in)       - grid structure
   fp_out(in)  - output file (printing comments on initialization, assumed open)
   phi(in)     - main level set function
   phase(in)   - phase < 0 indicates NW phase, phase > 0 indicates W phase (LSMPQS convention)
   mask(in)    - masking function (control volume)
   mask_phase(in) -relevant phase of the masking function
   
   Returns:
     pointer to the TrappedPhase structure
*/   
    
TrappedPhase *allocateAndInitializeTrappedPhase(
         Options *options,
	 Grid *g,
	 FILE *fp_out,
	 QSSLIB_REAL *phi,
	 QSSLIB_REAL phase,
	 QSSLIB_REAL *mask,
	 QSSLIB_REAL  mask_phase);
	 
/*! \fn  findPhaseConnectedComponentsAndRecordTrappedControlVolume()
   Find the number of connected components of the set is described by either
   positive or negative part of level set function. Trapped (disconnected)
   components are recorded in a separate array.

   Arguments
   phi(in)           - level set function
   phase(in)         - 1.0 or -1.0 depending on whether positive or negative level
                       set function values should be considered
   control_vol(in)   - control volume (mask) level set function
   control_vol_phase(in) - phase of control volume to consider - only points in this phase 
                        will be taken into account
   phi_trapped(in/out)   - array to record trapped phase
                          the trapped phase voxels will be set to the same 'phase'
			  therefore initialize array with the opposite values
   data_trapped(in/out)  - unsigned char array; trapped phase voxels will be
                       marked TRAP; also used for grass-fire search
   g(in)             - pointer to Grid structure
   stencil_size(in)  - voxel connectivity type (4,8 (2D),6 or 26(3D)), i.e. number of 
                       neighborhood voxels to be examined for each voxel
   min_comp_size     - minimal number of component's voxels (otherwise it's
                       disregarded)

   Returns
     number of connected components
     
   Note: function works for both 2D and 3D, as long as connectivity (stencil_size) is chosen
         properly     
*/
int findPhaseConnectedComponentsAndRecordTrappedControlVolume(
     QSSLIB_REAL   *phi,
     QSSLIB_REAL   *control_vol,
     QSSLIB_REAL   control_vol_phase,
     TrappedPhase  *tp,
     Grid          *g);
     
/*! \fn mergeTrappedComponentIntoMainPhase()
   
   Arguments
   search_seed_idx(in) - index of a voxel in the connected component; the component
                        is assumed marked TRAP in data_trapped array
   phi(in)             - level set function
   phase(in)           - 1.0 or -1.0 depending on whether positive or negative level
                       set function values should be considered
   phi_trapped(in/out)  - array that records trapped phase, assumed allocated/initialized
   data_trapped(in/out) - unsigned char array; trapped phase voxels assumed marked TRAP
   num_trapped_vox(in/out) - overall number of trapped voxels, will be reduced by the number
                        of voxels merged into main phase
   g(in)             - pointer to Grid structure
   stencil_size(in)  - voxel connectivity type (6 or 26), i.e. number of neighborhood
                       voxels to be examined for each voxel

   Returns		       
     number of merged voxels
     
   Note: function works for both 2D and 3D, as long as connectivity (stencil_size) is chosen
         properly    		       
*/
void mergeTrappedComponentIntoMainPhase(
     QSSLIB_REAL   *phi,
     TrappedPhase  *tp,
     Grid          *g);
	 
/* \fn checkForTrappedBlobs */
int checkForTrappedBlobs(
     TrappedPhase     *tp,
     QSS_DataArrays   *p,
     Grid             *g,
     QSSLIB_REAL      mask_sign,
     FILE             *fp_out);

/*! \fn destroyTrappedPhase() frees memory for TrappedPhase structure */
void destroyTrappedPhase(TrappedPhase *tp);     					 
			     
#endif

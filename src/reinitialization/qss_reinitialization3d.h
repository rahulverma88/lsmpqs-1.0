/*
 * File:        lsm_reinitialization3d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Header file for Fortran 77 3D reinitialization routines
 */

#ifndef INCLUDED_QSS_REINITIALIZATION_3D_H
#define INCLUDED_QSS_REINITIALIZATION_3D_H

#include "QSSLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_reinitialization3d.h
 *
 * \brief
 * @ref lsm_reinitialization3d.h provides support for computing the right-hand
 * side of the reinitialization and orthogonalization equations in three
 * space dimensions.
 *
 */


/* Link between C/C++ and Fortran function names 
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define COMPUTE_REINITIALIZATION_EQN_RHS_3D                              \
                                     computereinitializationeqnrhs3d_
#define QSS3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS                             \
                                     qss3dcomputeorthogonalizationeqnrhs_
#define QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX                              \
                                       qss3dcomputedistanceforsubcellfix_	
#define QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1           \
                      qss3dcomputereinitializationeqnrhssubcellfixorder1_

void qss_reinitialize_mask(QSS_DataArrays *, Grid *, Options *);

/*!
 * LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS() computes the right-hand side 
 * of the reinitialization equation using a Godunov scheme to select the
 * numerical discretization of the \f$ sgn(\phi) |\nabla \phi| \f$ term.
 * Forward (plus) and backward (minus) spatial derivatives used in
 * the Godunov calculation must be supplied by the user.
 *
 * Arguments:
 *  - reinit_rhs (out):       right-hand side of reinitialization
 *                            equation
 *  - phi (in):               level set function at current iteration
 *                            of reinitialization process
 *  - phi0 (in):              level set function at initial iteration
 *                            iteration of reinitialization process
 *  - phi_*_plus (in):        forward spatial derivatives for 
 *                            \f$ \nabla \phi \f$
 *  - phi_*_minus (in):       backward spatial derivatives for 
 *                            \f$ \nabla \phi \f$
 *  - use_phi0_for_sgn (in):  flag to specify whether \f$ \phi_0 \f$ should 
 *                            be used in the computation of \f$ sgn(\phi) \f$.
 *                              0 = use \f$ \phi \f$ (do NOT use 
 *                                  \f$ \phi_0 \f$ );
 *                              1 = use \f$ \phi_0 \f$.
 *  - *_gb (in):              index range for ghostbox
 *  - *_fb (in):              index range for fillbox
 *
 * NOTES:
 * - if use_phi0_for_sgn is not equal to 0 or 1, the default
 *   behavior of LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS() is to use
 *   \f$ \phi \f$ (i.e. equivalent to setting use_phi0_for_sgn to 0)
 *
 */
void COMPUTE_REINITIALIZATION_EQN_RHS_3D(
  QSSLIB_REAL *reinit_rhs,
  const QSSLIB_REAL* phi,
  const QSSLIB_REAL* phi0,
  const int *ilo_phi0_gb, 
  const int *ihi_phi0_gb, 
  const int *jlo_phi0_gb, 
  const int *jhi_phi0_gb,
  const int *klo_phi0_gb, 
  const int *khi_phi0_gb,
  const int *ilo_fb,   
  const int *ihi_fb,
  const int *jlo_fb,   
  const int *jhi_fb,
  const int *klo_fb,   
  const int *khi_fb,
  const QSSLIB_REAL *dx, 
  const QSSLIB_REAL *dy,
  const QSSLIB_REAL *dz,
  const int *use_phi0_for_sgn);


/*!
*  LSM3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX  determines whether each 
*  gridpoint is close to the interface (within one cell) and if so, 
*  computes the distance that is used in subcell fix to reinitialization. 
*  Following Russo, Smereka: "A Remark on Computing Distance Functions", 
*  J Comp Phys, 2000.
*
*  This ought to be done only once for phi0 and provided to routine
*  LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1
*  Implementation note: still not sure if it's worth storing another array
*  for this or recompute every time... 
*
*  Distance of functions away from the interface is set to a large negative 
*  number.
*
*  Arguments
*    distance0(out):         distance computed for phi0
*    phi0 (in):              level set function at initial iteration
*                            iteration of reinitialization process
*    *_gb (in):              index range for ghostbox
*    *_fb (in):              index range for fillbox
*    d* (in):                grid spacing
*/
void QSS3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX(
  QSSLIB_REAL* distance0, 
  const QSSLIB_REAL* phi0,
  const int *ilo_phi0_gb, 
  const int *ihi_phi0_gb, 
  const int *jlo_phi0_gb, 
  const int *jhi_phi0_gb,
  const int *klo_phi0_gb, 
  const int *khi_phi0_gb,
  const int *ilo_fb,   
  const int *ihi_fb,
  const int *jlo_fb,   
  const int *jhi_fb,
  const int *klo_fb,   
  const int *khi_fb,
  const QSSLIB_REAL *dx, 
  const QSSLIB_REAL *dy, 
  const QSSLIB_REAL *dz);
  
/*!
*  LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1() computes the right-hand 
*  side of the reinitialization equation.  Following 
*  Russo, Smereka: "A Remark on Computing Distance Functions", J Comp Phys,
*  2000.
*  1) Away from the interface 
*     It is using a Godunov scheme to select the 
*     numerical discretization of the sgn(phi) |grad(phi)| term.
*     Forward (plus) and backward (minus) spatial derivatives used in 
*     the Godunov calculation must be supplied by the user.
*  2) Near interface a correction (of first order) is made in order
*     not to use the information from the opposite side of the interface.
* 
*
*  Arguments:
*    reinit_rhs (out):       right-hand side of reinitialization 
*                            equation
*    phi (in):               level set function at current iteration
*                            of reinitialization process
*    phi0 (in):              level set function at initial iteration
*                            iteration of reinitialization process
*    distance0(in):           distance computed by lsm3DComputeDistanceSubcellFix
*                            at the beginning of reinitialization process
*    phi_*_plus (in):        forward spatial derivatives for grad(phi)
*    phi_*_minus (in):       backward spatial derivatives for grad(phi)
*    *_gb (in):              index range for ghostbox
*    *_fb (in):              index range for fillbox
*    d* (in):                grid spacing
*/  
void QSS3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1(
  QSSLIB_REAL *reinit_rhs,
  const QSSLIB_REAL* phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb, 
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const QSSLIB_REAL* phi0,
  const QSSLIB_REAL* distance0,
  const int *ilo_fb,   
  const int *ihi_fb,
  const int *jlo_fb,   
  const int *jhi_fb,
  const int *klo_fb,   
  const int *khi_fb,
  const QSSLIB_REAL *dx, 
  const QSSLIB_REAL *dy,
  const QSSLIB_REAL *dz);

#ifdef __cplusplus
}
#endif

#endif

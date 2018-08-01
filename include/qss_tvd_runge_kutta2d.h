#ifndef INCLUDED_QSS_TVD_RUNGE_KUTTA_2D_H
#define INCLUDED_QSS_TVD_RUNGE_KUTTA_2D_H

#include "QSSLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file qss_tvd_runge_kutta3d.h
 *
 * \brief
 * @ref qss_tvd_runge_kutta3d.h provides support for time integration of
 * partial differential equations in three space dimensions via 
 * total-variation diminishing Runge-Kutta methods.  Support is provided 
 * for first-, second-, and third-order time integration.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define QSS2D_RK1_STEP                      qss2drk1step_
#define QSS2D_TVD_RK2_STAGE2                qss2dtvdrk2stage2_
#define QSS2D_TVD_RK3_STAGE3                qss2dtvdrk3stage3_

/*!
 * QSS3D_RK1_STEP() takes a single first-order Runge-Kutta (i.e. Forward
 * Euler) step.
 *
 * Arguments:
 *  - u_next (out):  u(t_cur+dt)
 *  - u_cur (in):    u(t_cur)
 *  - rhs (in):      right-hand side of time evolution equation
 *  - dt (in):       step size
 *  - *_gb (in):     index range for ghostbox
 *  - *_fb (in):     index range for fillbox
 *
 * Return value:     none
 */
void QSS2D_RK1_STEP(
  QSSLIB_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const QSSLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const QSSLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const QSSLIB_REAL *dt);

/*
    For second stage of TVD Runge-Kutta.
*/
void QSS2D_TVD_RK2_STAGE2(
  QSSLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  QSSLIB_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const QSSLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const QSSLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const QSSLIB_REAL *dt);

void QSS2D_TVD_RK3_STAGE3(
  QSSLIB_REAL *u_stage3,
  const int *ilo_u_stage3_gb,
  const int *ihi_u_stage3_gb,
  const int *jlo_u_stage3_gb,
  const int *jhi_u_stage3_gb,
  QSSLIB_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const QSSLIB_REAL *u_cur,
  const int *ilo_u_cur_gb,
  const int *ihi_u_cur_gb,
  const int *jlo_u_cur_gb,
  const int *jhi_u_cur_gb,
  const QSSLIB_REAL *rhs, 
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const QSSLIB_REAL *dt);


#endif

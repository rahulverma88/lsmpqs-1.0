#ifndef included_qss_util_h
#define included_qss_util_h

#include "QSSLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
 
#define IMPOSE_MASK_PAR            imposemaskparallel_
#define IMPOSE_MASK_PAR_VAR        imposemaskparallelvar_
#define COPY_DATA_PAR_3D           copydataparallel3d_
#define QSS3D_MAX_NORM_DIFF            qss3dmaxnormdiff_
#define QSS3D_MAX_NORM_DIFF_LOCAL            qss3dmaxnormdifflocal_
#define QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO                               \
                                       qss3dvolumeregionphilessthanzero_
#define QSS3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO                            \
                                       qss3dvolumeregionphigreaterthanzero_
#define IMPOSE_TRAP_VEL_3D          imposetrapvelocity3d_
                                       
void IMPOSE_MASK_PAR(
    QSSLIB_REAL *phi_masked,
    QSSLIB_REAL *mask,
    QSSLIB_REAL *phi,
    QSSLIB_REAL *overlap,
    const int *ilo_gb,
    const int *ihi_gb, 
    const int *jlo_gb,
    const int *jhi_gb,
    const int *klo_gb,
    const int *khi_gb,
    const int *cur_klo_gb,
    const int *cur_khi_gb
);

void IMPOSE_MASK_PAR_VAR(
    QSSLIB_REAL *phi_masked,
    QSSLIB_REAL *mask,
    QSSLIB_REAL *phi,
    QSSLIB_REAL *overlap,
    const int *ilo_gb,
    const int *ihi_gb, 
    const int *jlo_gb,
    const int *jhi_gb,
    const int *klo_gb,
    const int *khi_gb,
    const int *cur_klo_gb,
    const int *cur_khi_gb
);

void COPY_DATA_PAR_3D(
    QSSLIB_REAL *phi_copy,
    QSSLIB_REAL *phi_orig,
    const int *ilo_gb,
    const int *ihi_gb, 
    const int *jlo_gb,
    const int *jhi_gb,
    const int *klo_gb,
    const int *khi_gb,
    const int *cur_klo_gb,
    const int *cur_khi_gb
);

/*!
 * QSS3D_MAX_NORM_DIFF() computes the max norm of the difference
 * between the two specified scalar fields.
 *      
 * Arguments:
 *  - max_norm_diff (out):   max norm of the difference between the fields
 *  - field1 (in):           scalar field 1
 *  - field2 (in):           scalar field 2
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for box to include in norm
 *                           calculation
 *
 * Return value:             none
 *
 */
void QSS3D_MAX_NORM_DIFF(
  QSSLIB_REAL *max_norm_diff,
  const QSSLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const int *klo_field1_gb, 
  const int *khi_field1_gb,
  const QSSLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const int *klo_field2_gb, 
  const int *khi_field2_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib);
  
void QSS3D_MAX_NORM_DIFF_LOCAL(
  QSSLIB_REAL *max_norm_diff,
  const QSSLIB_REAL *field1,
  const int *ilo_field1_gb, 
  const int *ihi_field1_gb,
  const int *jlo_field1_gb, 
  const int *jhi_field1_gb,
  const int *klo_field1_gb, 
  const int *khi_field1_gb,
  const QSSLIB_REAL *field2,
  const int *ilo_field2_gb, 
  const int *ihi_field2_gb,
  const int *jlo_field2_gb, 
  const int *jhi_field2_gb,
  const int *klo_field2_gb, 
  const int *khi_field2_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  QSSLIB_REAL *local_zone);

/*!
 * QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO() computes the volume of the
 * region where the level set function is less than 0.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi < 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void QSS3D_VOLUME_REGION_PHI_LESS_THAN_ZERO(
  QSSLIB_REAL *volume,
  const QSSLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const QSSLIB_REAL *dx,
  const QSSLIB_REAL *dy,
  const QSSLIB_REAL *dz,
  const QSSLIB_REAL *epsilon);


/*! 
 * QSS3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO() computes the volume of 
 * the region where the level set function is greater than 0.
 *
 * Arguments:
 *  - volume (out):          volume of the region where \f$ \phi > 0 \f$
 *  - phi (in):              level set function
 *  - dx, dy, dz (in):       grid spacing
 *  - epsilon (in):          width of numerical smoothing to use for 
 *                           Heaviside function
 *  - *_gb (in):             index range for ghostbox
 *  - *_ib (in):             index range for interior box
 *
 * Return value:             none
 *
 */
void QSS3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO(
  QSSLIB_REAL *volume,
  const QSSLIB_REAL *phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const QSSLIB_REAL *dx,
  const QSSLIB_REAL *dy,
  const QSSLIB_REAL *dz,
  const QSSLIB_REAL *epsilon);

void IMPOSE_TRAP_VEL_3D(
  QSSLIB_REAL   *phi_w,
  QSSLIB_REAL   *phi_nw,
  QSSLIB_REAL   *normal_velocity,
  QSSLIB_REAL   *curvature_coeff,
  QSSLIB_REAL   *vel_x,
  QSSLIB_REAL   *vel_y,
  QSSLIB_REAL   *vel_z,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb,
  const int *klo_phi_gb, 
  const int *khi_phi_gb,
  const int *ilo_ib, 
  const int *ihi_ib,
  const int *jlo_ib, 
  const int *jhi_ib,
  const int *klo_ib, 
  const int *khi_ib,
  const QSSLIB_REAL *dx);
  
#ifdef __cplusplus
}
#endif

#endif
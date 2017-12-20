#ifndef INCLUDED_QSS_GENERAL_UTIL_H
#define INCLUDED_QSS_GENERAL_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "qss_data_arrays.h"
#include "qss_options.h"
#include "qss_grid.h"

void reinitialize3d_subcell_fix_qss(QSS_DataArrays *, Grid *, Options *);
void reinitialize2d_subcell_fix_qss(QSS_DataArrays *, Grid *, Options *);
void qss_reinitialize_mask(QSS_DataArrays *, Grid *, Options *);
void qss_reinitialize_mask_no_bc(QSS_DataArrays *, Grid *, Options *);
void qss_reinitialize_mask2d(QSS_DataArrays *, Grid *, Options *);
void signedLinearExtrapolationBCqss(QSSLIB_REAL *, Grid *, int);
void initializeDisconnectedMasks(QSSLIB_REAL *, Grid *);
void qss_reinitializeDisconnectedMask2d(QSS_DataArrays *,Grid *,Options *);
void createReservoirInlet3d(QSS_DataArrays *, Grid *);
void readSpheresBinaryFile(char *,int  *,QSSLIB_REAL *,QSSLIB_REAL *,
    QSSLIB_REAL **, QSSLIB_REAL **, QSSLIB_REAL **, QSSLIB_REAL **);
void writeSpheresBinaryFile(char *, int ,int , QSSLIB_REAL *, QSSLIB_REAL *,
    QSSLIB_REAL *, QSSLIB_REAL *, QSSLIB_REAL *, QSSLIB_REAL *);
void setSpaceDerivFunc(Grid *g, Options *options);
void setThetaOverlap(QSS_DataArrays *, Grid *, Options *);
int return_1();
int return_0();
double return_0_double();

#endif
/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file constCurvModel3d.h

    Function headers for performing level set computations to reach constant curvature solutions
    in three dimensions.
    
    Both drainage and imbibition simulations call these functions.

    Sets up time-stepping parameters based on options structure, and then sets up parallelization.
    Domain decomposition for parallelization is done based on a simple bread-slicing 
    plan in the z-direction. 
    
    There are two functions - constCurvModel3d, and constCurvModel3dNoVar.
    The latter is for zero contact angle cases, where "a" and "b" parameters in the level
    set equation do not need to vary, and the convective "V" term does not exist. So not creating
    matrices for these terms results in significant memory and computational time savings,
    especially for large geometries in 3 dimensions.
             
*/

#ifndef INCLUDED_CONST_CURV_MODEL3D_H
#define INCLUDED_CONST_CURV_MODEL3D_H

QSSLIB_REAL  constCurvModel3d(Options *,QSS_DataArrays  *,Grid  *,FILE *);
QSSLIB_REAL  constCurvModel3dNoVar(Options *,QSS_DataArrays  *,Grid  *,FILE *, QSSLIB_REAL );

#endif
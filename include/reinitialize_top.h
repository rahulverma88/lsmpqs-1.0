/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file reinitialize_top.h
 *
 * Description: Function headers for reinitialization algorithms. Both 2D and 3D.
 *      reinitializeMedium* implements the reinitialization algorithm utilizing the 
 *      Godunov discretization. The SubcellFix variant is supposed to be more accurate,
 *      with less mass losses.
 */
#ifndef INCLUDED_REINIT_TOP_H
#define INCLUDED_REINIT_TOP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "qss_data_arrays.h"
#include "qss_options.h"
#include "qss_grid.h"

void reinitializeMedium2d(QSSLIB_REAL *, Grid *, Options *);
void reinitializeMedium3d(QSSLIB_REAL *, Grid *, Options *);
void reinitializeSubcellFix2d(QSSLIB_REAL *, Grid *, Options *);
void reinitializeSubcellFix3d(QSSLIB_REAL *, Grid *, Options *);

#ifdef __cplusplus
}
#endif

#endif
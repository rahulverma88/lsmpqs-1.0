/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file connectivity.h

    Header definitions for connectivity functions.
    
*/

#ifndef INCLUDED_CONNECTIVITY_H
#define INCLUDED_CONNECTIVITY_H

#include "qss_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void findConnectivity2d(QSS_DataArrays *, Grid *);
void findConnectivity3d(QSS_DataArrays *, Grid *);
void getMainInd(Options *o, QSS_DataArrays *p, Grid *g);
void getMainInd_new(Options *o, QSS_DataArrays *p, Grid *g);
void trapComponents(QSS_DataArrays *p, Grid *g, Options *o, QSSLIB_REAL val);
void trapComponents_mask(QSS_DataArrays *p, Grid *g, Options *o);

#ifdef __cplusplus
}
#endif



#endif
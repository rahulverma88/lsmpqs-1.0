#ifndef INCLUDED_CONNECTIVITY_H
#define INCLUDED_CONNECTIVITY_H

#include "qss_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void findConnectivity2d(QSS_DataArrays *, Grid *);
void findConnectivity3d(QSS_DataArrays *, Grid *);
void getMainInd(Options *o, QSS_DataArrays *p, Grid *g);
void trapComponents(QSS_DataArrays *p, Grid *g, Options *o);

#ifdef __cplusplus
}
#endif



#endif
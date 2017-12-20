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
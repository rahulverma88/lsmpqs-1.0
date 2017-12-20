#ifndef INCLUDED_DRAIN_TOP_H
#define INCLUDED_DRAIN_TOP_H

#define TOL       1e-16
#define EMAX_STOP 0.05
#define TPLOT     0.1

#include "qss_options.h"
#include "qss_data_arrays.h"

int    drainTop(Options *,char *,char *,char *);

#endif

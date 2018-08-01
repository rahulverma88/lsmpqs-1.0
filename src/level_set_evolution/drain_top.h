/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file drain_top.h

    Top level function header for calling the drainage constant curvature model in both 2D
    and 3D. 
    
    Reads in input data, makes sure everything is consistent. Essentially works as an error
    check for input data. Also writes out final data step.
             
*/


#ifndef INCLUDED_DRAIN_TOP_H
#define INCLUDED_DRAIN_TOP_H

#define TOL       1e-16
#define EMAX_STOP 0.05
#define TPLOT     0.1

#include "qss_options.h"
#include "qss_data_arrays.h"

int    drainTop(Options *,char *,char *,char *);

#endif

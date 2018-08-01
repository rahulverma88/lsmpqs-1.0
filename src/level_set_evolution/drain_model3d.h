/******************************************************************************
 *
 *   Author:   Rahul Verma
 *   Copyright (c) 2018, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
/*! \file drain_model3d.h

    Top level function header for calling the drainage constant curvature model in 3D. 
    
    This function does several things:
    
    1. Initializes and re-initializes all masks (solid space, disconnected components).
    2. Sets up the contact angle model - for constant theta, it sets up the variational a
        and b matrices. For spatially varying theta, it reads in a "theta.gz" file which 
        provides values of theta at different points in the domain. Based on the theta
        information, this function then determines overlap values at different points in 
        the domain between the pore space and the solid space.
    3. Calculates gradients of the solid space mask- this will be required in the level set
        computations, and needs to be only computed once since our mask doesn't change.
    4. For zero theta, it calls the less-memory intensive constCurvModel2dNoVar, while for
        all other cases, constCurvModel2d is called.
    5. It updates the "a" term for the next step.
    6. It writes out the data steps for each step.
             
*/

#ifndef INCLUDED_DRAIN_MODEL3D_H
#define INCLUDED_DRAIN_MODEL3D_H

void  drain_model3d(Options *,QSS_DataArrays  *,Grid  *,FILE *);

#endif
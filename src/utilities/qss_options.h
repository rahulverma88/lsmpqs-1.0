
#ifndef INCLUDED_QSS_OPTIONS
#define INCLUDED_QSS_OPTIONS

#ifdef __cplusplus
extern "C" {
#endif


/*! \file qss_options.h
 *
 * \brief
 * QSSLIB provides support for Options structure.
 *
 */

#include "QSSLIB_config.h"
#include "qss_grid.h"

/*! 
 * Structure 'Options' stores various input and running options
 * that can be helpful when coding.
 * Options can be set by user via reading an input file, and each
 * option that is not mentioned in input file will be set to default.
 *   
 * Modify this structure and accompanying functions 
 * - setDefaultOptions
 * - copyOptions
 * - createOptionsFromInputFile
 * - printOptions
 * to add any additional options you need.
 *
 */
typedef struct _Options
{
      /* Main part of the structure */
   char   outfile[256]; /* output file name */
   char   path[256];    /* output file path, used internally */
     
   QSSLIB_REAL dx;           /* grid spacing - assumed same in all dimensions */
   QSSLIB_REAL tmax;         /* max running time */
   
   /*! max running time for reinitalization*/
   QSSLIB_REAL   tmax_r;
   
    /*! time interval for evaluating max.abs.error, checkpointing etc. */
   QSSLIB_REAL  tplot;
   
   char   accuracy[20]; /* accuracy options: low, medium, high, very_high */
   QSSLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy_id;  
                        /* internal accuracy identifier */
   
   /* Assuming phi_t + a |grad_phi| = b kappa |grad_phi|
      where kappa is mean curvature
   */   
   QSSLIB_REAL a;   /* normal velocity */
   QSSLIB_REAL b;   /* mean curvature term */
   
   int    save_data;    /* save all data (1) or not (0); data saved in
                           the same directory as the  ouput file */
   int    do_reinit;    /* reinitialize periodically (1) or not (0) */
   int    do_mask;      /* impose mask (restrict domain of movement) (1) 
                           or not (0); if 1 then set 'mask' array in
			   QSS_DataArrays structure */

   QSSLIB_REAL eps_stop;
   
   /*! use subcell fix for reinitalization functions */ 
   int    subcell_fix;
   
   /*! reinitalization function, set internally */      			   
   void  (*reinit_func)();
   
   /*! curvature increment(decrement) */
   QSSLIB_REAL   dc;
   
   /*! maximum curvature for simulation */
   QSSLIB_REAL   cmax;
   
   /*! maximum amax for simulation */
   QSSLIB_REAL   amax;
   
   /*! minimum amin for simulation */
   QSSLIB_REAL   amin;
     
   /*! maximum NW phase volume fraction for simulation */
   QSSLIB_REAL   vol_frac_max;
   
    /*! minimum NW phase volume fraction for simulation */
   QSSLIB_REAL   vol_frac_min;
   
   /*! reinitialize mask (1) or not (0) */           
   int    mask_reinit;
   	   
   int    narrow_band;      /* use narrow banding or no */	
   
   /* First step of the simulation, if re-starting a stopped simulation */
   int    init_step;

   /* For non-zero, uniform theta, this sets the overlap between the pore 
        space and grain */
   QSSLIB_REAL overlap;
   
   /*! stop the simulation when opp. bdry. is touched - specific to drain or imbibe */
   int    stop_touch;	

   /*! eps (the  width of smearing for Heaviside function, where applicable), 
       is set to eps_coefficient*dx  */
   QSSLIB_REAL   eps_coefficient;
   
   /*! boundary condition on volume edges - 0 (signedLinear),
      1 (linear), 2 (copy), 3 (zeroLevelSet) */	        
   int    extrapol_id;
   
   /*! volume edge boundary condition function, set internally */      			   
   void  (*extrapol_func)();
   
   /* Function pointer to set the spatial derivative function */
   void  (*space_deriv_func)();
   
   /* Variables for storing order of spatial and temporal accuracy */
   int order_space_accur;
   int order_time_accur;
   
   /*! if 1, main level set function data is output every TPLOT */          
   int    checkpoint;	
   	   
   /*! Reintroduce fluid phase (non-wetting) at inlet - specific to 
       drain and mimics a non-wetting phase reservoir at inlet. */
   int    reservoir_inlet;    
    
   /* Contact angle in degrees */
   QSSLIB_REAL theta;
   
   /* Constant C used in variational method */
   QSSLIB_REAL C;
   
   /* Variables to store max(b)/dx and max(U)/dx, used for computing time step 
        Not really options. */
   QSSLIB_REAL b_max_over_dx, max_U_over_dx;
   
   /* Check connectivity or not */
   int check_connectivity;
   
   /* Variables to store fixed indices for wetting and non-wetting phase reservoirs
        at outlet/inlet, respectively */
   int phi_w_ind;
   int phi_nw_ind;
   
   /* Variable to store how far away from the zero level set to compute max_error.
        Useful for achieving convergence faster. */
   QSSLIB_REAL err_check_zone;
   
   /* If 1, then looks for 'theta.gz' in directory, for heterogeneous-wet case */
   int use_var_theta;
   
   /* Used for generating 'theta.gz' for cases where only two theta values need to be 
        imposed. */
   QSSLIB_REAL theta_1;
   QSSLIB_REAL theta_2;
   
   int    print_details;    /* whether to print details (1) or not (0) */ 
   int    use_satn_stop;    /* Use saturation as a stopping criteria: 
                                very small changes in satn stop the simulation */
   int    conserve_imbibe;  /* Conserve wetting phase during imbibition (1) or not (0) */
   int    center_inlet;     /* Mainly used for 2D Juanes micromodel case where inlet is at center
                            of domain, and outlet is all four sides. */
} Options;


/*!
 * Accuracy_settings structure should make it easier to switch between 
 * numerical discretization of different accuracy.
 * 'low'       <-->  HJ ENO1  in space, TVD Runge Kutta 1 in time
 * 'medium'    <-->  HJ ENO2  in space, TVD Runge Kutta 2 in time
 * 'high'      <-->  HJ ENO3  in space, TVD Runge Kutta 3 in time
 * 'very_high' <-->  HJ WENO5 in space, TVD Runge Kutta 3 in time
 */
typedef struct _Accuracy_settings_item {
   char *name;
   int  num_ghostcells;      /* number of ghostcells to be added to each
                                volume side */
} Accuracy_settings_item;


/*!
 * Accuracy_settings_menu connects the accuracy (set by user) with 
 * appropriate number of ghostcells to be added and the pointer
 * to the function that sets Grid for the appropriate numerical scheme.
 */
static Accuracy_settings_item Accuracy_settings_menu[] =
{
   {"low",       2 },
   {"medium",    3 },
   {"high",      5 },
   {"very_high", 4 }
};


/*================= Options structure manipulation ==================*/


/*!
 * createOptionsDefault() allocates a new Options structure and sets all 
 *  variables in Options structure to default values.
 *  
 * Arguments:          none
 *  
 * Return value:       pointer to an Options structure
 *
 * NOTES: Change this function if more Options structure elements are added.
 */
Options *createOptionsDefault(void);


/*!
 *  copyOptions() allocates a new Options structure and copies all 
 *  variables in Options structure from the provided Options structure.
 *  
 * Arguments:          
 *   - options_src(in):    Options structure to be copied
 *  
 * Return value:       pointer to an Options structure
 *
 * NOTES: Change this function if more Options structure elements are added.
 */
Options *copyOptions(Options *options_src);



/*! 'createOptionsFromInputFile' allocated memory for a new Options structure and sets its
 *         variables according to input file. Variables not set by input file
 *         are set to their default values.
 *
 *   Arguments:
 *   - filename (in):   input file name
 *   
 *   Return value:      pointer to a new Options structure 
 *
 *   NOTES: 
 *   - Input file is assumed to have a parameter entered in each line.
 *   - Appropriate keyword (equal to variable name from structure 'Options' )
 *     has to be at the beginnning of the line.
 *   - Change this function if more Options structure elements are added.
 */
Options *createOptionsFromInputFile(char *options);


/*!
 * 'printOptions' prints options structure variables to an output file
 *
 *   Arguments:
 *   - options(in):    Options structure to be printed
 *   - fp(in):         output file pointer (file opened for writing beforehand)
 *
 *   Return value:    none
 *
 *   NOTES: Change this function if more Options structure elements are added.
 */ 
void     printOptions(Options *options,FILE *fp);


/*!
 * destroyOptions() frees memory for an Options structure.
 *
 * Arguments:     
 *  - Options (in):   pointer to Options
 *
 * Return value:  none
 *
 */
void destroyOptions(Options *options);


/* User additions */

/* end User additions */

#ifdef __cplusplus
}
#endif

#endif

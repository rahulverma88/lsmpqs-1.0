#ifndef included_QSSLIB_config_h
#define included_QSSLIB_config_h

/* Debugging macro */
#ifndef QSSLIB_DEBUG_NO_INLINE
/* #undef LSMLIB_DEBUG_NO_INLINE  */
#endif

/* Macro defined if double precision library is being built. */
#ifndef QSSLIB_DOUBLE_PRECISION
/* #undef LSMLIB_DOUBLE_PRECISION */
#endif

/* Floating-point precision for LSMLIB_REAL */
#ifndef QSSLIB_REAL
#define QSSLIB_REAL float
#endif

#ifndef LSMLIB_REAL
#define LSMLIB_REAL float
#endif

/* Zero tolerance */
#ifndef QSSLIB_ZERO_TOL
#define QSSLIB_ZERO_TOL 1.e-5
#endif

/* Maximum value for LSMLIB_REAL */
#ifndef QSSLIB_REAL_MAX
#define QSSLIB_REAL_MAX FLT_MAX
#endif

/* Minimum value for LSMLIB_REAL */
#ifndef QSSLIB_REAL_MIN
#define QSSLIB_REAL_MIN FLT_MIN
#endif

/* Machine epsilon value for LSMLIB_REAL */
#ifndef QSSLIB_REAL_EPSILON
#define QSSLIB_REAL_EPSILON FLT_EPSILON
#endif

#define GB_DIMS  &(g->ilo_gb),&(g->ihi_gb),&(g->jlo_gb),&(g->jhi_gb),&(g->klo_gb),&(g->khi_gb)
#define FB_DIMS  &(g->ilo_fb),&(g->ihi_fb),&(g->jlo_fb),&(g->jhi_fb),&(g->klo_fb),&(g->khi_fb)
#define FB_DIMS_PAR  &(g->ilo_fb),&(g->ihi_fb),&(g->jlo_fb),&(g->jhi_fb),&(cur_klo_fb),&(cur_khi_fb)

#define GB_DIMS_2D  &(g->ilo_gb),&(g->ihi_gb),&(g->jlo_gb),&(g->jhi_gb)
#define FB_DIMS_2D  &(g->ilo_fb),&(g->ihi_fb),&(g->jlo_fb),&(g->jhi_fb)
#define FB_DIMS_PAR_2D  &(g->ilo_fb),&(g->ihi_fb),&(cur_jlo_fb),&(cur_jhi_fb)

#define DSZ  sizeof(QSSLIB_REAL)
#define ISZ  sizeof(int)
#define UCSZ sizeof(unsigned char)


#if defined(_OPENMP)
#else
    #define omp_get_thread_num return_0
    #define omp_get_num_threads return_1
    #define omp_get_wtime return_0_double
#endif

#endif


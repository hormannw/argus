
/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

#ifdef __GNUC__
#  define ATTRIBUTE__UNUSED  __attribute__ ((unused))
#else
#  define ATTRIBUTE__UNUSED
#endif

/*---------------------------------------------------------------------------*/

SEXP rargusRoU(SEXP sexp_n, SEXP sexp_chi);
/*---------------------------------------------------------------------------*/
/* Draw sample from argus distribution.                                      */
/* Wrapper for do_rargusRoU():                                              */
/*   GetRNGstate(); do_rargusRoU(...); PutRNGstate();                       */
/*---------------------------------------------------------------------------*/

SEXP do_rargusRoU(int n, double *chi, int nchi);
/*---------------------------------------------------------------------------*/
/* Draw sample from argus distribution                                       */
/* without calling GetRNGstate() ... PutRNGstate()                           */
/*---------------------------------------------------------------------------*/

SEXP dargus(SEXP sexp_x, SEXP sexp_chi, SEXP sexp_logvalue);
/*---------------------------------------------------------------------------*/
/* evaluate pdf of argus distribution                                        */
/*---------------------------------------------------------------------------*/

SEXP pargus(SEXP sexp_x, SEXP sexp_chi, SEXP sexp_lower, SEXP sexp_logvalue);
/*---------------------------------------------------------------------------*/
/* evaluate CDF of argus distribution                                        */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

#ifdef __GNUC__
#  define ATTRIBUTE__UNUSED  __attribute__ ((unused))
#else
#  define ATTRIBUTE__UNUSED
#endif

/*---------------------------------------------------------------------------*/

SEXP rargus(SEXP sexp_n, SEXP sexp_chi);
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                        */
/* Wrapper for do_rgig():                                                    */
/*   GetRNGstate(); do_rgig(...); PutRNGstate();                             */
/*---------------------------------------------------------------------------*/

SEXP do_rargus(int n, double lambda, double chi);
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution                                         */
/* without calling GetRNGstate() ... PutRNGstate()                           */
/*---------------------------------------------------------------------------*/

SEXP dargus(SEXP sexp_x, SEXP sexp_chi, SEXP sexp_logvalue);
/*---------------------------------------------------------------------------*/
/* evaluate pdf of GIG distribution                                          */
/*---------------------------------------------------------------------------*/


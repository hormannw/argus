
/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "argus.h"




/*****************************************************************************/
/* API                                                                       */
/*****************************************************************************/

SEXP dargus(SEXP sexp_x, SEXP sexp_chi, SEXP sexp_logvalue)
/*---------------------------------------------------------------------------*/
/* Evaluate density of argus distribution                                    */
/*                                                                           */
/* Parameters:                                                               */
/*   x ....... argument(s) for density                                       */
/*   chi   ... parameter vector with same length as x                        */
/*   logvalue  if TRUE the logarithm of the density will be returned         */
/*                                                                           */
/* Return:                                                                   */
/*   density at x                                                            */
/*---------------------------------------------------------------------------*/
{ 
  /* Arguments */
  double *x;                 /* array of points */
  int nx;                    /* length of x */
  double *chi;                /* parameter */
  int nchi;                    /* length of chi */
  int logvalue;
  double NORMCONSTANT;       /* normalization constant for density */
  int i;
  SEXP sexp_res;             /* results of computations */
  double *res;
/*  double res_err;            * result in case error */
/*  int err;                   * indicate invalid input */

  /* get array of points for which the density has to be evaluated */
  nx = length(sexp_x);
  PROTECT(sexp_x = AS_NUMERIC(sexp_x));
  x = REAL(sexp_x);

  /* compute logarithm? */
  logvalue = *(LOGICAL( AS_LOGICAL(sexp_logvalue) ));

  /* extract chi vector */
  nchi = length(sexp_chi);
  PROTECT(sexp_chi = AS_NUMERIC(sexp_chi));
  chi = REAL(sexp_chi);
  
  /* allocate array for results of density evaluations */
  PROTECT(sexp_res = NEW_NUMERIC(nx));
  res = REAL(sexp_res);

  /* evaluate density */
  if(nchi!=1){/* case nchi > 1 */
   for(i=0; i<nx; i++) {
     if (ISNAN(x[i])) {
       res[i] = x[i];
     }
     else if (!R_FINITE(x[i]) || x[i] <= 0. || x[i]>= 1) {
       res[i] =  (logvalue) ? R_NegInf : 0.;
     }
     else {
   /* compute normalization constant */
       NORMCONSTANT = chi[i]*chi[i]*chi[i]/(sqrt(2*M_PI)*(pnorm(chi[i],0.,1.,1,0)-chi[i]*dnorm(chi[i],0.,1.,0)-0.5));
       res[i] = NORMCONSTANT*x[i]*sqrt(1-x[i]*x[i])*exp(-0.5*chi[i]*chi[i]*(1-x[i]*x[i]));
       if (logvalue) res[i] = log(res[i]);
     }
   }
  }else{/* case nchi=1 , same chi value used for all x values*/
   /* compute normalization constant */
       NORMCONSTANT = chi[0]*chi[0]*chi[0]/(sqrt(2*M_PI)*(pnorm(chi[0],0.,1.,1,0)-chi[0]*dnorm(chi[0],0.,1.,0)-0.5));
   for(i=0; i<nx; i++) {
     if (ISNAN(x[i])) {
       res[i] = x[i];
     }
     else if (!R_FINITE(x[i]) || x[i] <= 0. || x[i]>= 1) {
       res[i] =  (logvalue) ? R_NegInf : 0.;
     }
     else { 
       res[i] = NORMCONSTANT*x[i]*sqrt(1-x[i]*x[i])*exp(-0.5*chi[0]*chi[0]*(1-x[i]*x[i]));
       if (logvalue) res[i] = log(res[i]);
     }
   }
  }

  /* return result to R */
  UNPROTECT(3);
  return sexp_res;
  
} /* end of dargus() */

/*---------------------------------------------------------------------------*/

/*
Rcode für pdf

parg<- function(x,chi){
psi <- function(chi){pnorm(chi)-chi*dnorm(chi)-0.5}	
1-psi(chi*sqrt(1-x*x))/psi(chi)	
}
parg((1:10)/10,1:10)

*/

double psi(double chi)
{ return( pnorm(chi,0.,1.,1,0)-chi*dnorm(chi,0.,1.,0)-0.5);}	

SEXP pargus(SEXP sexp_x, SEXP sexp_chi, SEXP sexp_lower, SEXP sexp_logvalue)
/*---------------------------------------------------------------------------*/
/* Evaluate CDF of argus distribution.                                       */
/*                                                                           */
/* Parameters:                                                               */
/*   x ....... argument(s) for density                                       */
/*   chi   ... parameter vector with same length as x                        */
/*   logvalue  if TRUE the logarithm of the density will be returned         */
/*                                                                           */
/* Return:                                                                   */
/*   density at x                                                            */
/*---------------------------------------------------------------------------*/
{ 
  /* Arguments */
  double *x;                 /* array of points */
  int nx;                    /* length of x */
  double *chi;               /* parameter */
  int nchi;                  /* length of chi */
  int lower;
  int logvalue;
  double PSICHI; /* \Psi(\chi) */
/*  double chiTsqrt1mx2;  chi*sqrt(1-chi^2) */
  int i;
  SEXP sexp_res;             /* results of computations */
  double *res;
/*  double res_err;            * result in case error */
/*  int err;                   * indicate invalid input */

  /* get array of points for which the density has to be evaluated */
  nx = length(sexp_x);
  PROTECT(sexp_x = AS_NUMERIC(sexp_x));
  x = REAL(sexp_x);

  /* compute lower tail? */
  lower = *(LOGICAL( AS_LOGICAL(sexp_lower) ));

  /* compute logarithm? */
  logvalue = *(LOGICAL( AS_LOGICAL(sexp_logvalue) ));

  /* extract chi vector */
  nchi = length(sexp_chi);
  PROTECT(sexp_chi = AS_NUMERIC(sexp_chi));
  chi = REAL(sexp_chi);
  
  /* allocate array for results of density evaluations */
  PROTECT(sexp_res = NEW_NUMERIC(nx));
  res = REAL(sexp_res);

  /* evaluate density */
  if(nchi!=1){/* case nchi > 1 */
   for(i=0; i<nx; i++) {
     if (ISNAN(x[i])) {
       res[i] = x[i];
     }
     else if (!R_FINITE(x[i]) || x[i] <= 0. || x[i]>= 1.) {
      if(x[i] <= 0.) res[i] =  (logvalue) ? R_NegInf : 0.;
      else res[i] =  (logvalue) ? 0. : 1.;
     }
     else {
   /* compute constant */
	 PSICHI = psi(chi[i]);
       if (!lower) x[i] = 1-x[i]; /* as the argus distribution has domain (0,1) */
	   res[i] = 1- psi(chi[i]*sqrt(1-x[i]*x[i]))/PSICHI;
       if (logvalue) res[i] = log(res[i]);
     }
   }
  }else{/* case nchi=1 , same chi value used for all x values*/
   /* compute constant */
	 PSICHI = psi(chi[0]);
   for(i=0; i<nx; i++) {
     if (ISNAN(x[i])) {
       res[i] = x[i];
     }
     else if (!R_FINITE(x[i]) || x[i] <= 0. || x[i]>= 1. ) {
       if(x[i] <= 0.) res[i] =  (logvalue) ? R_NegInf : 0.;
       else res[i] =  (logvalue) ? 0. : 1.;
     }
     else { 
       if (!lower) x[i] = 1-x[i]; /* as the argus distribution has domain (0,1) */
	   res[i] = 1- psi(chi[0]*sqrt(1-x[i]*x[i]))/PSICHI;
       if (logvalue) res[i] = log(res[i]);
     }
   }
  }

  /* return result to R */
  UNPROTECT(3);
  return sexp_res;
  
} /* end of pargus() */

/*---------------------------------------------------------------------------*/




SEXP rargusRoU(SEXP sexp_n, SEXP sexp_chi)
/*---------------------------------------------------------------------------*/
/* Draw sample from argus distribution.                                      */
/* Wrapper for do_rargusRoU() with GetRNGstate() ... PutRNGstate()          */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   chi   ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n*length(chi)'                                   */
/*---------------------------------------------------------------------------*/
{
  int n;          /* sample size */
  double *chi;                /* parameter */
  int nchi;                    /* length of chi */
  SEXP sexp_res;  /* results */


  /* extract sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));

  /* extract chi vector */
  nchi = length(sexp_chi);
  PROTECT(sexp_chi = AS_NUMERIC(sexp_chi));
  chi = REAL(sexp_chi);

  /* Get state of R uniform PRNG */
  GetRNGstate();

  /* run generator */
  PROTECT(sexp_res = do_rargusRoU(n, chi, nchi));

  /* Return state of PRNG to R */
  PutRNGstate();

  /* return result to R */
  UNPROTECT(2);
  return (sexp_res);

} /* end of rargusRoU() */

/*---------------------------------------------------------------------------*/

SEXP do_rargusRoU(int n, double *chi, int nchi)
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                        */
/* without calling GetRNGstate() ... PutRNGstate()                           */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   chi  ... parameter vector                                               */
/*   nchi ... length ofparameter vector                                      */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n*nchi                                           */
/*---------------------------------------------------------------------------*/
{
  SEXP sexp_res;           /* results */
  double *res;
  int i,j;

  /* check sample size */
  if (n<0) {
    error("sample size 'n' must be non-negative integer.");
  }

  /* allocate array for random sample */
  PROTECT(sexp_res = NEW_NUMERIC(n*nchi));
  res = REAL(sexp_res);

/* setup forRoU */    
  double m,b,ap,am,xp,xm,U,V,Y;
  for (j=0; j<nchi; j++) {
    if(chi[j]<=1){
      m = chi[j]*chi[j]*0.5;
      b = sqrt(sqrt(0.5)*chi[j])*exp(-chi[j]*chi[j]*0.25);
/*    b = 2^-0.25*sqrt(chi[j])*exp(-chi[j]*chi[j]/4) */
      ap = 0.;
      xm = (chi[j]*chi[j]+5.-sqrt(chi[j]*chi[j]*(chi[j]*chi[j]+6.)+25.))*0.25;
      am = (xm-chi[j]*chi[j]*0.5)*sqrt(sqrt(xm)*exp(-xm));
    }else{
      m = 0.5;
      b = exp(-(1.+log(2.))/4.);
/*  b <- 2^-0.25*exp(-1/4) */
      xp = chi[j]*chi[j]*0.5;
      if(xp > 1.5+sqrt(2.))xp = 1.5+sqrt(2.); 
/*  xp <- min(1.5 + sqrt(2),chi*chi*0.5) */
      xm = 1.5 - sqrt(2.);
      ap = (xp - 0.5)*sqrt(sqrt(xp)*exp(-xp));
      am = (xm - 0.5)*sqrt(sqrt(xm)*exp(-xm));
    }
    /* run generator */
    for (i=0; i<n; i++) {
      do {
/*      ++count; */
        U = b * unif_rand();      
        V = am+unif_rand()*(ap-am);
        Y = V/U + m ;
      }                              /* Acceptance/Rejection */
      while ( Y<0. || Y > chi[j]*chi[j]*0.5 || U*U > sqrt(Y)*exp(-Y)  );
      /* store random point */
      res[n*j+i] = sqrt(1.-2.*Y/(chi[j]*chi[j]));
    }
  }


  /* return result */
  UNPROTECT(1);
  return sexp_res;

} /* end of do_rargusRoU() */


#include<R.h>
#include<Rmath.h>

double mu = 0, sigma = 1;
int give_log = 0, lower_tail = 1;

//---------------------------------
double d_pgumbel();
double d_pgumbel2();
double d_pAO();
double d_plgamma();

double d_dgumbel();
double d_dgumbel2();
double d_dAO();
double d_dlgamma();

double d_gnorm();
double d_glogis();
double d_gcauchy();
double d_ggumbel();
double d_ggumbel2();
double d_glgamma();

//---

double d_pfun();
double d_dfun();
double d_gfun();
double d_gAO();

//--- negative log-likelihood:

double d_nll();

//--- Utilities:

double mmax();
double maxAbs();
void Trace();
//---------------------------------

//------------------------------------------------------------------
// CDFs:

double d_pgumbel(double q, double loc, double scale, int lower_tail)
{
  if(ISNAN(q)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1.;
  else if(q == R_NegInf)
    q = 0.;
  else {
    q = (q - loc) / scale;
    q = exp( -exp( -q));
  }
  return !lower_tail ? 1 - q : q;
}

void pgumbel(double *q, int *nq, double *loc, double *scale,
	     int *lower_tail)
{
// pgumbel()
    int i;
// How can I handle if loc and scale are not of unit length?
    for(i = 0; i < *nq; i++)
	q[i] = d_pgumbel(q[i], *loc, *scale, *lower_tail);
}

double d_pgumbel2(double q, double loc, double scale, int lower_tail)
// this is (partly) redundant since d_pgumbel2(q) = 1 - d_pgumbel(-q)
{
  if(ISNAN(q)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1;
  else if(q == R_NegInf)
    q = 0;
  else {
    q = (-q - loc) / scale;
    q = exp(-exp(-q));
  }
  return !lower_tail ? q : 1 - q;
}

void pgumbel2(double *q, int *nq, double *loc, double *scale,
	      int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = 1 - d_pgumbel(-q[i], *loc, *scale, *lower_tail);
}

double d_pAO(double q, double lambda, int lower_tail)
{
  if(ISNAN(q) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1;
  else if(q == R_NegInf)
    q = 0;
  else {
    if(lambda < 1.0e-6)
      error("'lambda' has to be positive. lambda = %e was supplied\n",
	    lambda);
    q = 1 - R_pow(lambda * exp(q) + 1, -1/lambda);
  }
  return !lower_tail ? 1 - q : q;
}

void pAO(double *q, int *nq, double *lambda, int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = d_pAO(q[i], *lambda, *lower_tail);
}

double d_plgamma(double eta, double lambda, int lower_tail)
{

  double v; 
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf) 
    v = 1;
  else if(eta == R_NegInf)
    v = 0;
  else {
    v = R_pow_di(lambda, -2) * exp(lambda * eta);
    if(lambda < 1.0e-6)
      v = 1 - pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		     lower_tail, give_log);
    if(lambda > -1.0e-6)
      v = pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		 lower_tail, give_log);
    if(lambda >= -1.0e-6 && lambda <= 1.0e-6)
      v = pnorm(eta, mu, sigma, lower_tail, give_log);
  }
  return !lower_tail ? 1 - v : v;
}

void plgamma(double *q, int *nq, double *lambda, int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = d_plgamma(q[i], *lambda, *lower_tail);
}

//------------------------------------------------------------------
// PDFs:

double d_dgumbel(double x, double loc, double scale, int give_log)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    // if(x == INFINITE || x == -INFINITE) // seems to work as well.
    return 0; // this special case needs to be handled separately 
  x = (x - loc) / scale;
  x = -exp(-x) - x - log(scale);
  return give_log ? x : exp(x);
}

void dgumbel(double *x, int *nx, double *loc, double *scale,
	     int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dgumbel(x[i], *loc, *scale, *give_log);
}

double d_dgumbel2(double x, double loc, double scale, int give_log)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0;
  x = (-x - loc) / scale;
  x = -exp(-x) - x - log(scale);
  return give_log ? x : exp(x);
}

void dgumbel2(double *x, int *nx, double *loc, double *scale,
	     int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dgumbel2(x[i], *loc, *scale, *give_log);
}

double d_dAO(double eta, double lambda, int give_log)
{
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf || eta == R_NegInf)
    return 0;
  if(lambda < 1.0e-6)
    error("'lambda' has to be positive. lambda = %e was supplied\n",
	  lambda);
  eta -= (1 + 1 / lambda) * log(lambda * exp(eta) + 1);
  return give_log ? eta : exp(eta);
}

void dAO(double *x, int *nx, double *lambda, int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dAO(x[i], *lambda, *give_log);
}

double d_dlgamma(double x, double lambda, int give_log)
{
  if(ISNAN(x) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 
  if(lambda < 1.0e-5 && lambda > -1.0e-5) // lambda close to zero
    return dnorm(x, mu, sigma, give_log); 

  double q_2 = R_pow_di(lambda, -2);
  x *= lambda;
  x = log(fabs(lambda)) + q_2 * log(q_2) -
    lgammafn(q_2) + q_2 * (x - exp(x));
  return !give_log ? exp(x) : x;
}

void dlgamma(double *x, int *nx, double *lambda, int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dlgamma(x[i], *lambda, *give_log);
}

//------------------------------------------------------------------
// gradients of PDFs:

double d_glogis(double x)
{
  // Gradient of dlogis(x) wrt. x
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    // if(x == INFINITE || x == -INFINITE) // seems to work as well.
    return 0; // this special case needs to be handled separately 

  /* Store the sign of x, compute the gradient for the absolute value
     and restore the sign. This is needed to avoid exp(LARGE) to blow
     up and the function to return NaN.
  */
  int sign = x > 0; //could use fsign() instead...
  x = exp(-fabs(x));
  x = 2 * x * x * R_pow_di(1 + x, -3) - x *
    R_pow_di(1 + x, -2);
  return sign ? x : -x;
}


void glogis(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_glogis(x[i]);
}

double d_gnorm(double x) {
  
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == INFINITY || x == -INFINITY)
    return 0;
  else
    return -x * dnorm(x, mu, sigma, give_log);
}

void gnorm(double *x, int *nx)
{
// Gradient of dnorm(x) wrt. x
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gnorm(x[i]);
}

double d_gcauchy(double x)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 
  
  return x = -2 * x / M_PI * R_pow_di(1 + x * x, -2);
}

void gcauchy(double *x, int *n)
{
// Gradient of dcauchy(x) wrt. x
    int i;
    for(i = 0; i < *n; i++)
	x[i] = d_gcauchy(x[i]);
}

double d_ggumbel(double x)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 

  x = exp(-x);
  if(x == INFINITY)
    return 0;
  
  double eq = exp(-x);
  return -eq * x + eq * x * x;
}

void ggumbel(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_ggumbel(x[i]);
}

double d_ggumbel2(double x)
// redundant function...
{
    return -d_ggumbel(-x);
}

void ggumbel2(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = -d_ggumbel(-x[i]);
    // or x[i] = d_ggumbel2(x[i]);
}

double d_gAO(double eta, double lambda)
{
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf || eta == R_NegInf)
    return 0; 

  double lex = lambda * exp(eta);
  if(lex == R_PosInf || lex == 0)
    return 0.0;
  double y = d_dAO(eta, lambda, give_log);
  
  return y == 0 ? 0.0 : y * (1 - (1 + 1/lambda) * lex / (1 + lex)); 
}

void gAO(double *x, int *nx, double *lambda)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gAO(x[i], *lambda);
}

double d_glgamma(double x, double lambda)
{
  if(ISNAN(x) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 
  if(lambda < 1.0e-5 && lambda > -1.0e-5) // lambda close to zero
    return -x * dnorm(x, mu, sigma, give_log);
  
  double z = exp(lambda * x);
  if(z == R_PosInf || z == 0)
    return 0.0;
  double y = d_dlgamma(x, lambda, give_log);
  return y <= 0 ? 0.0 : y * (1 - exp(lambda * x)) / lambda; 
  // Equivalent to:
  /* if(y <= 0) 
     return 0.0;
     else 
     return y * (1 - exp(lambda * x)) / lambda;
  */
}

void glgamma(double *x, int *nx, double *lambda)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_glgamma(x[i], *lambda);
}

//------------------------------------------------------------------
// link utility functions:

/* Link functions::
  1: logistic
  2: probit
  3: cloglog
  4: loglog
  5: cauchit
  6: Aranda-Ordaz
  7: log-gamma
 */

double d_pfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return plogis(x, mu, sigma, lower_tail, give_log);
    case 2: // probit
	return pnorm(x, mu, sigma, lower_tail, give_log);
    case 3: // cloglog
	return d_pgumbel(x, mu, sigma, lower_tail);
    case 4: // loglog
	return d_pgumbel2(x, mu, sigma, lower_tail);
    case 5: // cauchit
	return pcauchy(x, mu, sigma, lower_tail, give_log);
    case 6: // Aranda-Ordaz
	return d_pAO(x, lambda, lower_tail);
    case 7: // log-gamma
	return d_plgamma(x, lambda, lower_tail);
    default : // all other
	// if(link == 6)
	//     error("the Aranda-Ordaz link is not available");
	// if(link == 7)
	//     error("the log-gamma link is not available");
	// else
	error("link not recognized\n");
	return NA_REAL;
    }
}

void pfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_pfun(x[i], *lambda, *link);
}

double d_dfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return dlogis(x, mu, sigma, give_log);
    case 2: // probit
	return dnorm(x, mu, sigma, give_log);
    case 3: // cloglog
	return d_dgumbel(x, mu, sigma, give_log);
    case 4: // loglog
	return d_dgumbel2(x, mu, sigma, give_log);
    case 5: // cauchit
	return dcauchy(x, mu, sigma, give_log);
    case 6:
	return d_dAO(x, lambda, give_log);
    case 7:
	return d_dlgamma(x, lambda, give_log);
    default : // all other
	error("link not recognized\n");
	return NA_REAL;
    }
}

void dfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dfun(x[i], *lambda, *link);
}


double d_gfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return d_glogis(x);
    case 2: // probit
	return d_gnorm(x);
    case 3: // cloglog
	return d_ggumbel(x);
    case 4: // loglog
	return d_ggumbel2(x);
    case 5: // cauchit
	return d_gcauchy(x);
    case 6:
	return d_gAO(x, lambda);
    case 7:
	return d_glgamma(x, lambda);
    default : // all other
	error("link not recognized\n");
	return NA_REAL;
    }
}

void gfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gfun(x[i], *lambda, *link);
}

//------------------------------------------------------------------
// Gradients and Hessians for update.b in clmm():

void grFacSum(double *x, int *grFac, int *nx, double *u, int *nu)
// compute tapply(x, grFac, sum) + u
{
    int i, j;
    double z = 0;

    for(i = 0; i < *nu; i++) {
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                z = z + x[j];
	}
	u[i] = u[i] + z;
        z = 0;
    }
}

// FIXME: grFacSum such that it can be used by gradC and hessC - this
// should simplify the code

double d_nll(double *u, int nu, int *grFac, double stDev,
	     double *o1, double *o2, int no, double *eta1,
	     double *eta2, double *eta1Fix, double *eta2Fix,
	     double *sigma, double *pr, double *weights,
	     double lambda, int *link)
/*
  Returns:          nll
  Updates:          eta1, eta2, pr  given the new value of u
  Leaves unchanged: u, grFac, stDev, o1, o2, eta1Fix, eta2Fix, sigma, weights
 */
{
    int i, j;
    double o, nll = 0.0;

    for(i = 0; i < no; i++) {
	o = u[grFac[i] - 1] * stDev;
	eta1[i] = (eta1Fix[i] + o1[i] - o) / sigma[i];
	eta2[i] = (eta2Fix[i] + o2[i] - o) / sigma[i];
	pr[i] = d_pfun(eta1[i], lambda, *link) -
	    d_pfun(eta2[i], lambda, *link);
	if(!R_FINITE(pr[i]) || pr[i] <= 0.0) {
	    return INFINITY;
	}
	nll -= weights[i] * log(pr[i]);
    }
    for(j = 0; j < nu; j++)
	nll -= dnorm(u[j], 0.0, 1.0, 1);
    return nll;
}

void nll(double *u, int *nu, int *grFac, double *stDev,
	 double *o1, double *o2, int *no, double *eta1,
	 double *eta2, double *eta1Fix, double *eta2Fix,
	 double *sigma, double *pr, double *weights,
	 double *lambda, int *link, double *nll)
{
    *nll = d_nll(u, *nu, grFac, *stDev, o1, o2, *no, eta1, eta2,
		 eta1Fix, eta2Fix, sigma, pr, weights, *lambda, link);
}

void grad(double *stDev, double *p1, double *p2, double *pr,
	  double *weights, double *sigma, double *wtprSig,
	  double *eta1, double *eta2, double *gradValues,
	  double *u, int *grFac, int *nx, int *ngv,
	  double *lambda, int *link)
/*
  Returns:          void
  Updates:          gradValues, p1, p2, wtprSig  given the new values of eta1, eta2
  Leaves unchanged: grFac, stDev, eta1, eta2, pr, sigma, weights, link, nx, ngv
  Assumes:
  nx: length of grFac, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2
  ngv: length of gradValues
 */
{
    int i, j;
    // double tmp[*nx], z = 0;

    // update p1, p2, wtprSig:
    for(i = 0; i < *nx; i++) {
	p1[i] = d_dfun(eta1[i], *lambda, *link);
	p2[i] = d_dfun(eta2[i], *lambda, *link);
	wtprSig[i] = weights[i] / pr[i] / sigma[i];
    }

    // sum for each level of the grouping factor:
    for(i = 0; i < *ngv; i++) {
	gradValues[i] = 0; // Could set these to
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
		gradValues[i] += *stDev * wtprSig[j] *
		    (p1[j] - p2[j]);
	}
	gradValues[i] += u[i];
    }
}

void gradC(double *stDev, double *p1, double *p2, double *wtprSig,
	   int *grFac, int *nx, double *u, int *nu)
{
    // gradient for update.b
    int i, j;
    double z = 0;

    for(i = 0; i < *nx; i++) {
	wtprSig[i] = *stDev * wtprSig[i] * (p1[i] - p2[i]);
    }

    for(i = 0; i < *nu; i++) {
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                z += wtprSig[j];
	}
	u[i] += z;
        z = 0;
    }
}

void hess(double *stDev, double *p1, double *p2, double *pr,
	  double *wtprSig, double *eta1, double *eta2, int *link,
	  int *grFac, int *nx, double *hessValues, double *lambda,
	  int *nhv)
/*
  Returns:          void
  Updates:          hessValues  given the new values of eta1, eta2
  Leaves unchanged: grFac, stDev, eta1, eta2, p1, p2, pr, sigma, weights, link, nx, ngv
  Assumes:
  nx: length of grFac, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2
  nhv: length of hessValues
 */
{
    int i, j;

    // sum for each level of the grouping factor:
    for(i = 0; i < *nhv; i++) {
	hessValues[i] = 0;
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                hessValues[i] +=
		    (R_pow_di(p1[j] - p2[j], 2) / pr[j] -
		     d_gfun(eta1[j], *lambda, *link) +
		     d_gfun(eta2[j], *lambda, *link)) * wtprSig[j];
	}
	hessValues[i] = (hessValues[i] * *stDev * *stDev) + 1;
	// hessValues[i] *= *stDev * *stDev;
	// ++hessValues[i];
    }
}

void hessC(double *stDev, double *p1, double *p2, double *pr,
	   double *g1, double *g2, double *wtprSig,
	   int *grFac, int *nx, double *z, int *nz)
{
    // hessian for update.b
    int i, j;
    double sigma2;

    sigma2 = R_pow_di(*stDev, 2);

    for(i = 0; i < *nx; i++)
	pr[i] = (R_pow_di(p1[i] - p2[i], 2) / pr[i] -
		 g1[i] + g2[i]) * wtprSig[i];

    for(i = 0; i < *nz; i++) {
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                z[i] = z[i] + pr[j];
	}
	z[i] = z[i] * sigma2 + 1;
    }
}

//------------------------------------------------------------------
// Trace function:

void Trace(int iter, double stepFactor, double val, double maxGrad,
	   double *par, int npar, int first)
{
    int i;

    if(first)
	Rprintf("iter:  step factor:     Value:     max|grad|:   Parameters:\n");
    Rprintf(" %3d:    %1.3e:   %.3f:     %1.3e:  ", iter, stepFactor, val, maxGrad);
    for(i = 0; i < npar; i++)
	Rprintf(" %.4f", par[i]);
    Rprintf("\n");
}

//------------------------------------------------------------------

void NRalg(int *trace, int *maxIter, double *gradTol,
	   int *maxLineIter, int *grFac,
	   double *stDev, double *o1, double *o2,
	   double *eta1Fix, double *eta2Fix, double *eta1,
	   double *eta2, double *sigma, int *link,
	   double *weights, double *u,
	   double *pr, double *funValue,
	   double *gradValues, double *hessValues,
	   int *nx, int *nu, double *maxGrad, int *conv,
	   double *p1, double *p2, double *wtprSig,
	   double *lambda, int *Niter)
{
/*
  nx: length(pr)
  r:  length(start) = length(u)

  updates: u, funValue, gradValues, hessValues, maxGrad,

  correct vector input:
  eta1, eta2, pr, funValue (grad is called before d_nll), u = 0,
  grFac, o1, o2, eta1Fix, eta2Fix, sigma, weights

  arbitrary input:
  p1, p2, wtprSig, gradValues, hessValues,

  needed output:
  u, funValue, gradValues, hessValues, conv, Niter,
*/
    int lineIter, innerIter = 0, i, j;
    double stepFactor = 1, funValueTry, step[*nu];

    *funValue = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1, eta2,
		      eta1Fix, eta2Fix, sigma, pr, weights,
		      *lambda, link);
    if(!R_FINITE(*funValue)) {
	*conv = 0;
	return ;
    }
    grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	 gradValues, u, grFac, nx, nu, lambda, link);
    *maxGrad = maxAbs(gradValues, *nu);
    *conv = -1; // Convergence flag
    if(*trace)
	Trace(0, stepFactor, *funValue, *maxGrad, u, *nu, 1);

    // Newton-Raphson algorithm:
    for(i = 0; i < *maxIter; i++) {
        if(*maxGrad < *gradTol) {
            *conv = 1;
            return ;
	}
	hess(stDev, p1, p2, pr, wtprSig, eta1, eta2, link,
	     grFac, nx, hessValues, lambda, nu);
	for(j = 0; j < *nu; j++) {
	    step[j] = gradValues[j] / hessValues[j];
	    u[j] -= stepFactor * step[j];
	}
	funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
			    eta2, eta1Fix, eta2Fix, sigma, pr,
			    weights, *lambda, link);
	lineIter = 0;
	//  simple line search, i.e. step halfing:
	while(funValueTry > *funValue) {
	    stepFactor *= 0.5;
	    for(j = 0; j < *nu; j++)
		u[j] += stepFactor * step[j];
	    funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
				eta2, eta1Fix, eta2Fix, sigma, pr,
				weights, *lambda, link);
	    lineIter++;
	    if(*trace)
		Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad,
		      u, *nu, 0);
	    if(lineIter > *maxLineIter){
		*conv = -2;
		return ;
	    }
	    innerIter++;
        }
        *funValue = funValueTry;
	grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	     gradValues, u, grFac, nx, nu, lambda, link);
	*maxGrad = maxAbs(gradValues, *nu);
	if(*trace)
	    Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad, u,
		  *nu, 0);
	stepFactor = fmin2(1.0, stepFactor * 2.0);
	(*Niter)++;
    }
}

void NRalgv3(int *trace, int *maxIter, double *gradTol,
	     int *maxLineIter, int *grFac,
	     double *stDev, double *o1, double *o2,
	     double *eta1Fix, double *eta2Fix, double *sigma,
	     int *link, double *weights, double *u,
	     double *pr, double *funValue,
	     double *gradValues, double *hessValues,
	     int *nx, int *nu, double *maxGrad, int *conv,
	     double *lambda, int *Niter)
// Less input and slightly faster than NRalg().
{
/*
  control arguments from clmm - see ?clmm.control:
  trace, maxIter, gradTol, maxLineIter all of length 1

  length = nx: grFac, o1, o2, eta1Fix, eta2Fix, sigma, weights
  length = 1: stDev, funValue, nx, nu, maxGrad, conv, lambda, Niter
  length = nu: gradValues, hessValues, u

  updates: u, funValue, gradValues, hessValues, maxGrad, conv, Niter,
  pr,

  correct vector input:
  eta1, eta2, pr, u = 0, grFac, o1, o2, eta1Fix, eta2Fix, sigma,
  weights

  arbitrary input:
  gradValues, hessValues,

  needed output:
  u, funValue, gradValues, hessValues, conv, Niter,
*/
    int lineIter, innerIter = 0, i, j;
    double stepFactor = 1, funValueTry, step[*nu];
    double eta1[*nx], eta2[*nx], p1[*nx], p2[*nx], wtprSig[*nx];

    *funValue = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1, eta2,
		      eta1Fix, eta2Fix, sigma, pr, weights,
		      *lambda, link);
    if(!R_FINITE(*funValue)) {
	*conv = 0;
	return ;
    }
    grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	 gradValues, u, grFac, nx, nu, lambda, link);
    *maxGrad = maxAbs(gradValues, *nu);
    *conv = -1; // Convergence flag
    if(*trace)
	Trace(0, stepFactor, *funValue, *maxGrad, u, *nu, 1);

    // Newton-Raphson algorithm:
    for(i = 0; i < *maxIter; i++) {
        if(*maxGrad < *gradTol) {
            *conv = 1;
            return ;
	}
	hess(stDev, p1, p2, pr, wtprSig, eta1, eta2, link, grFac, nx,
	     hessValues, lambda, nu);
	for(j = 0; j < *nu; j++) {
	    /* Actually there is no need to store 'step' since
	       'gradValues' could hold the step values (maintained
	       here for code clarity) */
	    step[j] = gradValues[j] / hessValues[j];
	    u[j] -= stepFactor * step[j];
	}
	funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
			    eta2, eta1Fix, eta2Fix, sigma, pr,
			    weights, *lambda, link);
	lineIter = 0;
	//  simple line search, i.e. step halfing:
	while(funValueTry > *funValue) {
	    stepFactor *= 0.5;
	    for(j = 0; j < *nu; j++)
		u[j] += stepFactor * step[j];
	    funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
				eta2, eta1Fix, eta2Fix, sigma, pr,
				weights, *lambda, link);
	    lineIter++;
	    if(*trace)
		Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad,
		      u, *nu, 0);
	    if(lineIter > *maxLineIter){
		*conv = -2;
		return ;
	    }
	    innerIter++;
        }
        *funValue = funValueTry;
	grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	     gradValues, u, grFac, nx, nu, lambda, link);
	*maxGrad = maxAbs(gradValues, *nu);
	if(*trace)
	    Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad, u,
		  *nu, 0);
	stepFactor = fmin2(1.0, stepFactor * 2.0);
	(*Niter)++;
    }
}

//------------------------------------------------------------------

void getNGHQ(double *nll, int *grFac, double *stDev,
	     double *eta1Fix, double *eta2Fix, double *o1, double *o2,
	     double *Sigma, double *weights, int *nx, int *nu,
	     double *ghqns, /* double *ghqws,*/ double *lghqws,
	     int *nGHQ, int *link, double *ns, double *lambda)
{
    int i, j, h;
    double SS = 0, SS1 = 0, SS2 = 0, eta1tmp, eta2tmp;

    for(i = 0; i < *nu; i++) {
	for(h = 0; h < *nGHQ; h++) {
	    for(j = 0; j < *nx; j++) {
		if(grFac[j] == i + 1) {
		    eta1tmp = (eta1Fix[j] + o1[j] - ns[h]) / Sigma[j];
		    eta2tmp = (eta2Fix[j] + o2[j] - ns[h]) / Sigma[j];
		    SS1 += weights[j] *
			      log(d_pfun(eta1tmp, *lambda, *link) -
				  d_pfun(eta2tmp, *lambda, *link));
		}
	    }
	    // SS2 += exp(SS1) * ghqws[h];
	    // SS2 += exp(SS1 + log(ghqws[h]));
	    SS2 += exp(SS1 + lghqws[h]);
	    SS1 = 0;
	}
	SS += log(SS2);
	SS2 = 0;
    }
    *nll = -SS + *nu * log(M_PI * 2) * 0.5;
}

void getNAGQ(double *nll, int *grFac, double *stDev,
	     double *eta1Fix, double *eta2Fix, double *o1, double *o2,
	     double *Sigma, double *weights, int *nx, int *nu,
	     double *ghqns, double *lghqws, /* double *lghqws, */
	     double *ghqns2, double *u, double *D,
	     int *nAGQ, int *link, double *lambda)
/*
  nll: negative log-likelihood (return value)

  length = nx: grFac, o1, o2, eta1Fix, eta2Fix, Sigma, weights
  length = 1: stDev, nll, nx, nu, nAGQ, lambda, link
  length = nu: D, u
  length = nAGQ: ghqns, lghqws (log ghqws) / ghqws
 */
{
    int i, j, h;
    double SS1 = 0, SS2 = 0, eta1tmp, eta2tmp, K, ranNew;
    *nll = 0;

    for(i = 0; i < *nu; i++) {
	K = sqrt(2.0 / D[i]);
	for(h = 0; h < *nAGQ; h++) {
	    for(j = 0; j < *nx; j++) {
		if(grFac[j] == i + 1) {
		    ranNew = *stDev * (u[i] + K * ghqns[h]);
		    eta1tmp = (eta1Fix[j] + o1[j] - ranNew) / Sigma[j];
		    eta2tmp = (eta2Fix[j] + o2[j] - ranNew) / Sigma[j];
		    SS1 += weights[j] *
			      log(d_pfun(eta1tmp, *lambda, *link) -
				  d_pfun(eta2tmp, *lambda, *link));
		}
	    }
	    // SS2 += exp(SS1) * K * ghqws[h] *
	    // 	dnorm(u[i] + K * ghqns[h], mu, sigma, give_log);
//  	    SS2 += exp(SS1 + lghqws[h] + ghqns2[h] - //R_pow_di(ghqns[h], 2) +
//  		       0.5 * R_pow_di(u[i] + K * ghqns[h], 2)) * K;
 	    SS2 += exp(SS1 + lghqws[h] + ghqns2[h] - //R_pow_di(ghqns[h], 2) +
 		       0.5 * R_pow_di(u[i] + K * ghqns[h], 2));
	    SS1 = 0;
	}
	// *nll -= log(SS2);
	*nll -= log(SS2) + log(K);
	SS2 = 0;
    }
    *nll += *nu * log(M_PI * 2) * 0.5;
}


//------------------------------------------------------------------

double mmax(double *x, int nx)
/*
   Return the maximum of the elements in x
   nx: length of x ( >= 1)
 */
{
    int i;
    double cmax; // current max

    cmax = x[0];
    if(nx == 1)
	return cmax;
    for(i = 1; i < nx; i++) {
	if(x[i] > cmax)
	    cmax = x[i];
    }
    return cmax;
}

double maxAbs(double *x, int nx)
/*
  Return max(abs(x))
  nx: length of x ( >= 1 )
 */
{
    int i;
    double cmax; // current max

    cmax = fabs(x[0]);
    if(nx == 1)
	return cmax;
    for(i = 1; i < nx; i++) {
	if(fabs(x[i]) > cmax)
	    cmax = fabs(x[i]);
    }
    return cmax;
}


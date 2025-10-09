/**
 * @file RbMathFunctions
 * This file contains some useful (standard) math functions.
 *
 * @brief Implementation of standard math functions.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream> // IWYU pragma: keep

#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"


using namespace RevBayesCore;


// Math Functions


/*!
 * This function calculates the beta function.
 *
 * \brief Beta function.
 * \param a is an argument. 
 * \param b is an argument. 
 * \return Returns the value of the beta function. 
 * \throws Does not throw an error.
 */
double RbMath::beta(double a, double b) {
    
	return ( exp(RbMath::lnGamma(a) + RbMath::lnGamma(b) - RbMath::lnGamma(a + b)) );
}

/*!
 *  DESCRIPTION
 *	Evaluates the "deviance part"
 *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
 *		  =  x * log(x/M) + M - x
 *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
 *
 *	in a manner that should be stable (with small relative error)
 *	for all x and M=np. In particular for x/np close to 1, direct
 *	evaluation fails, and evaluation is based on the Taylor series
 *	of log((1+v)/(1-v)) with v = (x-np)/(x+np).
 */
double RbMath::binomialDeviance(double x, double np)
{
    double ej, s, s1, v;
    int j;
    
    if (!isFinite(x) || !isFinite(np) || np == 0.0)
    {
        std::ostringstream ss;
        ss << "Cannot compute the binomial deviance function for x = " << x << " and np = " << np;
        throw RbException(ss.str());
    }
    
    if (fabs(x-np) < 0.1*(x+np))
    {
        v = (x-np)/(x+np);
        s = (x-np)*v;/* s using v -- change by MM */
        ej = 2*x*v;
        v = v*v;
        for (j=1; ; j++) { /* Taylor series */
            ej *= v;
            s1 = s+ej/((j<<1)+1);
            if (s1==s) /* last term was effectively 0 */
                return(s1);
            s = s1;
        }
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}


/**
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */
int RbMath::chebyshev_init(double *dos, int nos, double eta)
{
    int i, ii;
    double err;
    
    if (nos < 1)
        return 0;
    
    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
        i = nos - ii;
        err += fabs(dos[i]);
        if (err > eta) {
            return i;
        }
    }
    return i;
}


double RbMath::chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;
    
    if (n < 1 || n > 1000) 
    {
        std::ostringstream s;
        s << "Cannot compute chebyshev function for n = " << n;
        throw RbException(s.str());
    }
    
    if (x < -1.1 || x > 1.1) 
    {
        std::ostringstream s;
        s << "Cannot compute chebyshev function for x = " << x;
        throw RbException(s.str());
    }
    
    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++)
    {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

/**
 * C++ version of the expm1 function. We provide our own, since this
 * function is not available in the Microsoft cmath header
 *
 * Compute exp(x) - 1 without loss of precision for small values of x.
 */
double RbMath::expm1(double x)
{
    if (fabs(x) < 1e-5)
        return x + 0.5*x*x;
    else
        return exp(x) - 1.0;
}

/*!
 * This function calculates the gamma function for real x.
 *
 * \brief Gamma function.
 * \param x is the argument. 
 * \return Returns the value of the gamma function. 
 * \throws Does not throw an error.
 */
double RbMath::gamma(double x) 
{
    const static double gamcs[42] = {
        +.8571195590989331421920062399942e-2,
        +.4415381324841006757191315771652e-2,
        +.5685043681599363378632664588789e-1,
        -.4219835396418560501012500186624e-2,
        +.1326808181212460220584006796352e-2,
        -.1893024529798880432523947023886e-3,
        +.3606925327441245256578082217225e-4,
        -.6056761904460864218485548290365e-5,
        +.1055829546302283344731823509093e-5,
        -.1811967365542384048291855891166e-6,
        +.3117724964715322277790254593169e-7,
        -.5354219639019687140874081024347e-8,
        +.9193275519859588946887786825940e-9,
        -.1577941280288339761767423273953e-9,
        +.2707980622934954543266540433089e-10,
        -.4646818653825730144081661058933e-11,
        +.7973350192007419656460767175359e-12,
        -.1368078209830916025799499172309e-12,
        +.2347319486563800657233471771688e-13,
        -.4027432614949066932766570534699e-14,
        +.6910051747372100912138336975257e-15,
        -.1185584500221992907052387126192e-15,
        +.2034148542496373955201026051932e-16,
        -.3490054341717405849274012949108e-17,
        +.5987993856485305567135051066026e-18,
        -.1027378057872228074490069778431e-18,
        +.1762702816060529824942759660748e-19,
        -.3024320653735306260958772112042e-20,
        +.5188914660218397839717833550506e-21,
        -.8902770842456576692449251601066e-22,
        +.1527474068493342602274596891306e-22,
        -.2620731256187362900257328332799e-23,
        +.4496464047830538670331046570666e-24,
        -.7714712731336877911703901525333e-25,
        +.1323635453126044036486572714666e-25,
        -.2270999412942928816702313813333e-26,
        +.3896418998003991449320816639999e-27,
        -.6685198115125953327792127999999e-28,
        +.1146998663140024384347613866666e-28,
        -.1967938586345134677295103999999e-29,
        +.3376448816585338090334890666666e-30,
        -.5793070335782135784625493333333e-31
    };
    
    int i, n;
    double y;
    double sinpiy, value;
    
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     * (xmin, xmax) are non-trivial, see ./gammalims.c
     * xsml = exp(.01)*DBL_MIN
     * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
     */
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
    
    if (isNan(x)) return x;
    
    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == (std::int64_t)x)) 
    {
        std::ostringstream s;
        s << "Cannot compute gamma function for x = " << x;
        throw RbException(s.str());
    }
    
    y = fabs(x);
    
    if (y <= 10) {
        
        /* Compute gamma(x) for -10 <= x <= 10
         * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
         * first of all. */
        
        n = int(x);
        if (x < 0) --n;
        y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
        --n;
        value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
        if (n == 0)
            return value;/* x = 1.dddd = 1+y */
        
        if (n < 0) {
            /* compute gamma(x) for -10 <= x < 1 */
            
            /* exact 0 or "-n" checked already above */
            
            /* The answer is less than half precision */
            /* because x too near a negative integer. */
            if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
                std::ostringstream s;
                s << "Cannot compute gamma function for x = " << x;
                throw RbException(s.str());
            }
            
            /* The argument is so close to 0 that the result would overflow. */
            if (y < xsml) {
                std::ostringstream s;
                s << "Cannot compute gamma function for x = " << x;
                throw RbException(s.str());
            }
            
            n = -n;
            
            for (i = 0; i < n; i++) {
                value /= (x + i);
            }
            return value;
        }
        else {
            /* gamma(x) for 2 <= x <= 10 */
            
            for (i = 1; i <= n; i++) {
                value *= (y + i);
            }
            return value;
        }
    }
    else {
        /* gamma(x) for	 y = |x| > 10. */
        
        if (x > xmax) {			/* Overflow */
            std::ostringstream s;
            s << "Cannot compute gamma function for x = " << x;
            throw RbException(s.str());
        }
        
        if (x < xmin) {			/* Underflow */
            std::ostringstream s;
            s << "Cannot compute gamma function for x = " << x;
            throw RbException(s.str());
        }
        
        if (y <= 50 && y == (int)y) { /* compute (n - 1)! */
            value = 1.;
            for (i = 2; i < y; i++) value *= i;
        }
        else { /* normal case */
            value = exp((y - 0.5) * log(y) - y + RbConstants::LN_SQRT_2PI +
                        ((2*y == (int)2*y)? stirlerr(y) : lnGammacor(y)));
        }
        if (x > 0)
            return value;
        
        if (fabs((x - (int)(x - 0.5))/x) < dxrel){
            
            /* The answer is less than half precision because */
            /* the argument is too near a negative integer. */
            
            std::ostringstream s;
            s << "Cannot compute gamma function for x = " << x;
            throw RbException(s.str());
        }
        
        sinpiy = sin(RbConstants::PI * y);
        if (sinpiy == 0) {		/* Negative integer arg - overflow */
            std::ostringstream s;
            s << "Cannot compute gamma function for x = " << x;
            throw RbException(s.str());
        }
        
        return -RbConstants::PI / (y * sinpiy * value);
    }

}

double RbMath::gamma_old(double x) {
    
	double c[7] = { -1.910444077728E-03, 
        8.4171387781295E-04, 
        -5.952379913043012E-04, 
        7.93650793500350248E-04, 
        -2.777777777777681622553E-03, 
        8.333333333333333331554247E-02, 
        5.7083835261E-03 };
	double fact;
	int i;
	int n;
	double p[8] = { -1.71618513886549492533811, 
        2.47656508055759199108314E+01, 
        -3.79804256470945635097577E+02, 
        6.29331155312818442661052E+02, 
        8.66966202790413211295064E+02, 
        -3.14512729688483675254357E+04, 
        -3.61444134186911729807069E+04, 
        6.64561438202405440627855E+04 };
	double q[8] = { -3.08402300119738975254353E+01, 
        3.15350626979604161529144E+02,
        -1.01515636749021914166146E+03,
        -3.10777167157231109440444E+03, 
        2.25381184209801510330112E+04, 
        4.75584627752788110767815E+03, 
        -1.34659959864969306392456E+05, 
        -1.15132259675553483497211E+05 };
	double sqrtpi = 0.9189385332046727417803297;
	double sum2;
	double value;
	double xbig = 35.040;
	double xden;
	double xminin = 1.18E-38;
	double xnum;
	double y;
	double y1;
	double ysq;
	double z;
    
	bool parity = false;
	fact = 1.0;
	n = 0;
	y = x;
    
	if ( y <= 0.0 )
    {
		/* argument negative */
		y = -x;
		y1 = ( double ) ( ( int ) ( y ) );
		value = y - y1;
        
		if ( value != 0.0 )
        {
			if ( y1 != ( double ) ( ( int ) ( y1 * 0.5 ) ) * 2.0 )
				parity = true;
			fact = -RbConstants::PI / std::sin(RbConstants::PI * value);
			y = y + 1.0;
        }
		else
        {
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
        }
    }
	if ( y < rbEpsilon() )
    {
		/* argument < EPS */
		if ( xminin <= y )
        {
			value = 1.0 / y;
        }
		else
        {
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
        }
    }
	else if ( y < 12.0 )
    {
		y1 = y;
		/* 0.0 < argument < 1.0 */
		if ( y < 1.0 )
        {
			z = y;
			y = y + 1.0;
        }
		/* 1.0 < argument < 12.0, reduce argument if necessary */
		else
        {
			n = int ( y ) - 1;
			y = y - ( double ) ( n );
			z = y - 1.0;
        }
		/* evaluate approximation for 1.0 < argument < 2.0 */
		xnum = 0.0;
		xden = 1.0;
		for ( i = 0; i < 8; i++ )
        {
			xnum = ( xnum + p[i] ) * z;
			xden = xden * z + q[i];
        }
        
		value = xnum / xden + 1.0;
		/* adjust result for case  0.0 < argument < 1.0 */
		if ( y1 < y )
        {
			value = value / y1;
        }
		/* adjust result for case  2.0 < argument < 12.0 */
		else if ( y < y1 )
        {
			for ( i = 1; i <= n; i++ )
            {
				value = value * y;
				y = y + 1.0;
            }
        }
    }
	else
    {
		/* evaluate for 12 <= argument */
		if ( y <= xbig )
        {
			ysq = y * y;
			sum2 = c[6];
			for ( i = 0; i < 6; i++ )
            {
				sum2 = sum2 / ysq + c[i];
            }
			sum2 = sum2 / y - y + sqrtpi;
			sum2 = sum2 + ( y - 0.5 ) * log ( y );
			value = exp ( sum2 );
        }
		else
        {
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
        }
        
    }
	/* final adjustments and return */
	if ( parity )
    {
		value = -value;
    }
	if ( fact != 1.0 )
    {
		value = fact / value;
    }
    
	return value;
}


/*!
 * This function returns the incomplete beta function, which is
 *
 * BI(a,b,x) = Integral(0 <= T <= X) T**(A-1) (1-T)**(B-1) dt /
 *             Integral(0 <= T <= 1) T**(A-1) (1-T)**(B-1) dt
 *
 * \brief Incomplete beta function.
 * \param a is a beta parameter. 
 * \param b is a beta parameter. 
 * \param x is the upper limit of integration. 
 * \return Returns the incomplete beta function. 
 * \throws Does not throw an error.
 * \see Majumder & Bhattacharjee. 1973. Algorithm AS63. Applied 
 *      Statistics, 22.
 */
double RbMath::incompleteBeta(double a, double b, double x) {
    
	double tol = 1.0E-07;
    
	double value;
	if ( x <= 0.0 )
    {
		value = 0.0;
		return value;
    }
	else if ( 1.0 <= x )
    {
		value = 1.0;
		return value;
    }
    
	/* change tail if necessary and determine S */
	double psq = a + b;
    
	double xx, cx, pp, qq;
	bool indx;
	if ( a < (a + b) * x )
    {
		xx = 1.0 - x;
		cx = x;
		pp = b;
		qq = a;
		indx = true;
    }
	else
    {
		xx = x;
		cx = 1.0 - x;
		pp = a;
		qq = b;
		indx = false;
    }
    
	double term = 1.0;
	int i = 1;
	value = 1.0;
	int ns = (int)(qq + cx * (a + b));
    
	/* use Soper's reduction formulas */
	double rx = xx / cx;
    
	double temp = qq - (double)i;
	if ( ns == 0 )
		rx = xx;
    
	int it = 0;
	int it_max = 1000;
	for (;;)
    {
		it++;
		if ( it_max < it )
        {
			return RbConstants::Double::inf;
        }
		term = term * temp * rx / ( pp + ( double ) ( i ) );
		value = value + term;
		temp = fabs(term);
		if ( temp <= tol && temp <= tol * value )
			break;
		i++;
		ns--;
		if ( 0 <= ns )
        {
			temp = qq - (double)i;
			if ( ns == 0 )
				rx = xx;
        }
		else
        {
			temp = psq;
			psq = psq + 1.0;
        }
    }
    
	/* finish calculation */
	value = value * exp(pp * log(xx) + (qq - 1.0) * log(cx) - lnBeta(a, b) - log(pp));
	if ( indx )
		value = 1.0 - value;
	return value;
}


/*!
* This function returns the incomplete gamma function
 * where x is the upper or lower limit of the integration
 * and alpha is the shape parameter.
 *
 * \brief Incomplete gamma function.
 * \param alpha is the shape parameter of the gamma. 
 * \param x is the upper limit of integration. 
 * \param regularized indicates the regularized incomplete gamma
 * \param lower indicates x is the lower limit of integration, upper otherwise
 * \return Returns the incomplete gamma function.
 * \throws RbException on error
 * \see Bhattacharjee, G. P. 1970. The incomplete gamma integral. Applied 
 *      Statistics, 19:285-287.
 */


double RbMath::incompleteGamma(double x, double alpha, bool regularized, bool lower)
{

    // (1) series expansion     if (alpha>x || x<=1)
    // (2) continued fraction   otherwise
    // RATNEST FORTRAN by
    // Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
    // 19: 285-287 (ASA032)
    
    double accurate = 1e-8, overflow = 1e30;
    double factor, gin, rn, a, b, an, dif, term;
    double pn0, pn1, pn2, pn3, pn4, pn5;
    
    if (x == 0.0)
    {
        return lower ? 0.0 : regularized ? 1.0 : RbMath::gamma(alpha);
    }
    if ( alpha <= 0.0 )
    {
        std::ostringstream s;
        s << "Cannot compute incomplete gamma function for alpha = " << alpha;
        throw RbException(s.str());
    }
    if ( x < 0 && RbMath::isInt(alpha) == false )
    {
        std::ostringstream s;
        s << "Cannot compute incomplete gamma function for x = " << x << " and alpha = " << alpha;
        throw RbException(s.str());
    }

    double scale = RbMath::lnGamma(alpha);
    factor = exp(alpha * log(fabs(x)) - x - scale) * pow(sgn(x), alpha);
    
    if (x > 1 && x >= alpha) {
        // continued fraction
        a = 1 - alpha;
        b = a + x + 1;
        term = 0;
        pn0 = 1;
        pn1 = x;
        pn2 = x + 1;
        pn3 = x * b;
        gin = pn2 / pn3;
        
        do {
            a++;
            b += 2;
            term++;
            an = a * term;
            pn4 = b * pn2 - an * pn0;
            pn5 = b * pn3 - an * pn1;
            
            if (pn5 != 0) {
                rn = pn4 / pn5;
                dif = fabs(gin - rn);
                if (dif <= accurate) {
                    if (dif <= accurate * rn) {
                        break;
                    }
                }
                
                gin = rn;
            }
            pn0 = pn2;
            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            if (fabs(pn4) >= overflow) {
                pn0 /= overflow;
                pn1 /= overflow;
                pn2 /= overflow;
                pn3 /= overflow;
            }
        } while (true);
        gin = 1.0 - factor * gin;
    } else {
        // series expansion
        gin = 1;
        term = 1;
        rn = alpha;
        do {
            rn++;
            term *= x / rn;
            gin += term;
        }
        while (fabs(term) > accurate);
        gin *= factor / alpha;
    }

    gin = lower ? gin : 1.0 - gin;
    gin = regularized ? gin : exp(log(fabs(gin)) + scale) * sgn(gin);
	
    return gin;
}


/*!
 * This function calculates the log of the beta function. The code was adapted from R.
 *
 * \brief Natural log of the beta function.
 * \param a is the alpha parameter of the beta function. 
 * \param b is the beta parameter of the beta function. 
 * \return Returns the log of the beta function. 
 * \throws Does not throw an error.
 * \see function lbeta() in R version 2.12.2
 */
double RbMath::lnBeta(double a, double b)
{
    double corr, p, q;
       
    p = q = a;
    if (b < p) p = b;/* := min(a,b) */
    if (b > q) q = b;/* := max(a,b) */
    
    /* both arguments must be >= 0 */
    if (p < 0)
    {
        std::ostringstream s;
        s << "Cannot compute log-beta function for a = " << a << " and b = " << b;
        throw RbException(s.str());
    }
    else if (p == 0) 
    {
        return RbConstants::Double::inf;
    }
    else if (!isFinite(q)) 
    { /* q == +Inf */
        return RbConstants::Double::neginf;
    }
        
    if (p >= 10) {
        /* p and q are big. */
        corr = lnGammacor(p) + lnGammacor(q) - lnGammacor(p + q);
        return log(q) * -0.5 + RbConstants::LN_SQRT_2PI + corr
            + (p - 0.5) * log(p / (p + q)) + q * RbMath::log1p(-p / (p + q));
    }
    else if (q >= 10) {
        /* p is small, but q is big. */
        corr = lnGammacor(q) - lnGammacor(p + q);
        return lnGamma(p) + corr + p - p * log(p + q)
            + (q - 0.5) * RbMath::log1p(-p / (p + q));
    }
    else
        /* p and q are small: p <= q < 10. */
        return log(gamma(p) * (gamma(q) / gamma(p + q)));
    
}

/*!
 * This function calculates the log of the gamma function, which is equal to:
 * Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
 * The result is accurate to 10 decimal places. Stirling's formula is used
 * for the central polynomial part of the procedure.
 *
 * \brief Natural log of the gamma function.
 * \param alp is the parameter of the gamma function. 
 * \return Returns the log of the gamma function. 
 * \throws Does not throw an error.
 * \see Pike, M. C. and I. D. Hill. 1966. Algorithm 291: Logarithm of the gamma
 *      function. Communications of the Association for Computing Machinery, 9:684.
 *
 *    The new adoption from R 2.12.2:
 *
 *    The function lgammafn computes log|gamma(x)|.  The function
 *    lgammafn_sign in addition assigns the sign of the gamma function
 *    to the address in the second argument if this is not NULL.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 */
double RbMath::lnGamma_sign(double x, int *sgn)
{
    double ans, y, sinpiy;
    
#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;
    
    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
        xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
        dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
     dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
     */
#undef xmax
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765696e-8
#endif
    
    if (sgn != NULL) *sgn = 1;
    
    if (x < 0 && fmod(floor(-x), 2.) == 0)
        if (sgn != NULL) *sgn = -1;
    
    if (x <= 0 && x == RbMath::trunc(x)) { /* Negative integer argument */
        std::ostringstream s;
        s << "Cannot compute log-gamma function for x = " << x;
        throw RbException(RbException::MATH_ERROR,s.str());
    }
    
    y = fabs(x);
    
    if (y <= 10)
        return log(fabs(gamma(x)));
    /*
     ELSE  y = |x| > 10 ---------------------- */
    
    if (y > xmax) {
        std::ostringstream s;
        s << "Cannot compute log-gamma function for x = " << x;
        throw RbException(RbException::MATH_ERROR,s.str());
    }
    
    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
        if (x > 1e17)
            return(x*(log(x) - 1.));
        else if (x > 4934720.)
            return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
        else
#endif
            return RbConstants::LN_SQRT_2PI + (x - 0.5) * log(x) - x + lnGammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sin(RbConstants::PI * y));
    
    ans = RbConstants::LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lnGammacor(y);
    
    if (fabs((x - RbMath::trunc(x - 0.5)) * ans / x) < dxrel) {
        
        /* The answer is less than half precision because
         * the argument is too near a negative integer. */
        
        std::ostringstream s;
        s << "Cannot compute log-gamma function for x = " << x;
        throw RbException(RbException::MATH_ERROR,s.str());
    }
    
    return ans;
}

double RbMath::lnGamma(double x)
{
    return lnGamma_sign(x, NULL);
}

double RbMath::lnGamma_old(double a) {
    
	double x = a;
	double f = 0.0;
	double z;
	if (x < 7) 
    {
		f = 1.0;  
		z = x - 1.0;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
    }
	z = 1.0 / (x*x);
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
             (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
              0.083333333333333)/x);  
}

/*!
 *  Compute the log gamma correction factor for x >= 10 so that
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */
double RbMath::lnGammacor(double x)
{
    const static double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
        -.3399615005417721944303330599666e-15,
        +.2683181998482698748957538846666e-17,
        -.2868042435334643284144622399999e-19,
        +.3962837061046434803679306666666e-21,
        -.6831888753985766870111999999999e-23,
        +.1429227355942498147573333333333e-24,
        -.3547598158101070547199999999999e-26,
        +.1025680058010470912000000000000e-27,
        -.3401102254316748799999999999999e-29,
        +.1276642195630062933333333333333e-30
    };
    
    double tmp;
    
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     *   xbig = 2 ^ 26.5
     *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#undef  xmax
#define xmax  3.745194030963158e306
    
    if (x < 10)
    {
        std::ostringstream s;
        s << "Cannot compute log-gammacor function for x = " << x;
        throw RbException(RbException::MATH_ERROR,s.str());
    }
    else if (x >= xmax)
    {
        std::ostringstream s;
        s << "Cannot compute log-gammacor function for x = " << x;
        throw RbException(RbException::MATH_ERROR,s.str());
        /* allow to underflow below */
    }
    else if (x < xbig)
    {
        tmp = 10 / x;
        return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double RbMath::logGamma1p(double x) 
{
    const double eulers_const =	 0.5772156649015328606065120900824024;
    
    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
        0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
        0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
        0.2058080842778454787900092413529198e-1,
        0.7385551028673985266273097291406834e-2,
        0.2890510330741523285752988298486755e-2,
        0.1192753911703260977113935692828109e-2,
        0.5096695247430424223356548135815582e-3,
        0.2231547584535793797614188036013401e-3,
        0.9945751278180853371459589003190170e-4,
        0.4492623673813314170020750240635786e-4,
        0.2050721277567069155316650397830591e-4,
        0.9439488275268395903987425104415055e-5,
        0.4374866789907487804181793223952411e-5,
        0.2039215753801366236781900709670839e-5,
        0.9551412130407419832857179772951265e-6,
        0.4492469198764566043294290331193655e-6,
        0.2120718480555466586923135901077628e-6,
        0.1004322482396809960872083050053344e-6,
        0.4769810169363980565760193417246730e-7,
        0.2271109460894316491031998116062124e-7,
        0.1083865921489695409107491757968159e-7,
        0.5183475041970046655121248647057669e-8,
        0.2483674543802478317185008663991718e-8,
        0.1192140140586091207442548202774640e-8,
        0.5731367241678862013330194857961011e-9,
        0.2759522885124233145178149692816341e-9,
        0.1330476437424448948149715720858008e-9,
        0.6422964563838100022082448087644648e-10,
        0.3104424774732227276239215783404066e-10,
        0.1502138408075414217093301048780668e-10,
        0.7275974480239079662504549924814047e-11,
        0.3527742476575915083615072228655483e-11,
        0.1711991790559617908601084114443031e-11,
        0.8315385841420284819798357793954418e-12,
        0.4042200525289440065536008957032895e-12,
        0.1966475631096616490411045679010286e-12,
        0.9573630387838555763782200936508615e-13,
        0.4664076026428374224576492565974577e-13,
        0.2273736960065972320633279596737272e-13,
        0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };
    
    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;
    double lgam;
    int i;
    
    if (fabs(x) >= 0.5)
        return RbMath::lnGamma (x + 1);
    
    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * logcf(-x / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--)
        lgam = coeffs[i] - x * lgam;
    
    return (x * lgam - eulers_const) * x - log1pmx (x);
} /* lgamma1p */

/**
 * C++ version of the log1p function. We provide our own, since this
 * function is not available in the Microsoft cmath header
 *
 * Compute log(1+x) without losing precision for small values of x
 */
double RbMath::log1p(double x)
{
    if (x <= -1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << x 
           << ") for log1p; must be greater than -1.0";
        throw RbException( os.str() );
    }

    if (fabs(x) > 1e-4)
    {
        // x is large enough that the obvious evaluation is OK
        return log(1.0 + x);
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8

    return (-0.5*x + 1.0)*x;
}




/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double RbMath::log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;
    
    if (x > 1 || x < minLog1Value)
        return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
            * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
            * ---------------------------------------------
            * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
            */
        double r = x / (2 + x), y = r * r;
        if (fabs(x) < 1e-2) {
            static const double two = 2;
            return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
                         two / 3) * y - x);
        } else {
            static const double tol_logcf = 1e-14;
            return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
        }
    }
}

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR


/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */
double RbMath::logcf (double x, double i, double d,
       double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;
    
#if 0
    assert (i > 0);
    assert (d >= 0);
#endif
    
    b2 = c4 * b1 - i * b2;
    
    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
        double c3 = c2*c2*x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;
        
        c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;
        
        if (fabs (b2) > scalefactor) {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        } else if (fabs (b2) < 1 / scalefactor) {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }
    
    return a2 / b2;
}

/*!
 * This function returns the round off unit for floating point arithmetic.
 * The returned value is a number, r, which is a power of 2 with the property
 * that, to the precision of the computer's arithmetic, 1 < 1 + r, but
 * 1 = 1 + r / 2. This function comes from John Burkardt.
 *
 * \brief Round off unity for floating point arithmetic.
 * \return Returns round off unity for floating point arithmetic. 
 * \throws Does not throw an error.
 */
double RbMath::rbEpsilon(void) {
    
	double r = 1.0;
	while ( 1.0 < (double)(1.0 + r)  )
		r = r / 2.0;
	return 2.0 * r;
}

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lnGammacor() which computes almost the same!
 */
double RbMath::stirlerr(double n)
{
    
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
    
    /*
     error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
     */
    const static double sferr_halves[31] = {
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,  /* 0.5 */
        0.0810614667953272582196702,  /* 1.0 */
        0.0548141210519176538961390,  /* 1.5 */
        0.0413406959554092940938221,  /* 2.0 */
        0.03316287351993628748511048, /* 2.5 */
        0.02767792568499833914878929, /* 3.0 */
        0.02374616365629749597132920, /* 3.5 */
        0.02079067210376509311152277, /* 4.0 */
        0.01848845053267318523077934, /* 4.5 */
        0.01664469118982119216319487, /* 5.0 */
        0.01513497322191737887351255, /* 5.5 */
        0.01387612882307074799874573, /* 6.0 */
        0.01281046524292022692424986, /* 6.5 */
        0.01189670994589177009505572, /* 7.0 */
        0.01110455975820691732662991, /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
    };
    double nn;
    
    if (n <= 15.0) {
        nn = n + n;
        if (nn == (int)nn) return(sferr_halves[(int)nn]);
        return(lnGamma(n + 1.) - (n + 0.5)*log(n) + n - RbConstants::LN_SQRT_2PI);
    }
    
    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

/**
 * C++ version of the trunc function. We provide our own, since this
 * function is not available in the Microsoft cmath header
 *
 * Truncate a floating-point value.
 */
double RbMath::trunc(double x)
{
    return double( int( x ) );
}

/**
 *
 * Calculate the log of the sum of the exponentials of a vector of numbers
 */
double RbMath::log_sum_exp(const std::vector<double>& x, double max)
{
    double lse = 0.0;
    for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); it++)
        lse += exp(*it - max);

    return max + log(lse);
}

/**
 *
 * Calculate the log of the sum of the exponentials of a vector of numbers
 */
double RbMath::log_sum_exp(const std::vector<double>& x)
{
    double max = x.front();
    for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); it++)
        max = std::max(max,*it);

    return log_sum_exp(x,max);
}



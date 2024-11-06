/**
 * @file RbMathCombinatorialFunctions
 * This file contains some useful (standard) combinatorial math functions.
 *
 * @brief Implementation of standard combinatorial math functions.
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
#include <cmath>
#include <iostream>
#include <sstream> // IWYU pragma: keep

#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"
#include "RbException.h"

using namespace RevBayesCore;

/**	Binomial coefficients.
*	choose(n, k)   and  lchoose(n,k) := log(abs(choose(n,k))
*
*	These work for the *generalized* binomial theorem,
*	i.e., are also defined for non-integer n  (integer k).
*
*  We use the simple explicit product formula for  k <= k_small_max
*  and also have added statements to make sure that the symmetry
*    (n \\ k ) == (n \\ n-k)  is preserved for non-negative integer n.
*
* Note: These functions are adapted from the R source code.
*/


#define k_small_max 30
/* 30 is somewhat arbitrary: it is on the *safe* side:
 * both speed and precision are clearly improved for k < 30.
 */
double RbMath::choose(double n, double k)
{

    double r, k0 = k;
    k = floor(k + 0.5);
    if (fabs(k - k0) > 1e-7)
    {
        std::ostringstream s;
        s << "'k' (" << k0 << ") must be integer, rounded to " << k;
        throw RbException(s.str());
    }
    if (k < k_small_max)
    {
        int j;
        if (n-k < k && n >= 0 && RbMath::isInt(n)) k = n-k; /* <- Symmetry */
        if (k <	 0) return 0.0;
        if (k == 0) return 1.0;
        /* else: k >= 1 */
        r = n;
        for (j = 2; j <= k; j++)
            r *= (n-j+1)/j;
        return RbMath::isInt(n) ? floor(r + 0.5) : r;
        /* might have got rounding errors */
    }
    /* else: k >= k_small_max */
    if (n < 0)
    {
        r = choose(-n+ k-1, k);
        if (k != 2 * floor(k / 2.)) r = -r;
        return r;
    }
    else if (RbMath::isInt(n))
    {
        if (n < k) return 0.;
        if (n - k < k_small_max) return choose(n, n-k); /* <- Symmetry */
        return floor(exp(RbMath::lfastchoose(n, k)) + 0.5);
    }
    /* else non-integer n >= 0 : */
    if (n < k-1)
    {
        int s_choose;
        r = RbMath::lfastchoose2(n, k, /* -> */ &s_choose);
        return s_choose * exp(r);
    }
    return exp(RbMath::lfastchoose(n, k));
}


/*!
 * Number of ways to choose 2 elements from k elements
 *
 * \brief Return Kc2 (k choose 2)
 * \param x The x value
 * \return K choose 2
 */
int RbMath::kchoose2(int k){

//    if (k==0 | k==1 ){
//        return 0;
//    }
//    else if (k==2){
//        return 1;
//    }
//    else{
        return (k*(k-1)/2);
//    }

}

/*!
 * This function returns the harmonic number of x=> sum_i=1^x(1/i)
 *
 * \brief Return sum_i=1^x(1/i)
 * \param x The x value
 * \return The harmonic number of x
 */
double RbMath::harmonicNumber(size_t n)
{
    
    double hm = 0.0;
    for (int i=1; i<=n; i++)
    {
        hm += 1.0/i;
    }
    
    return (hm);
}

/*!
 * This function returns the factorial of x, x!
 *
 * \brief Return x!
 * \param x The x value
 * \return The factorial x!
 */
double RbMath::factorial(int x)
{

    double fac = 1.0;
    for (int i=1; i<=x; i++)
    {
	fac *= i;
    }
    
    return (fac);
}

#define ODD(_K_) ((_K_) != 2 * floor((_K_) / 2.))
/* matching R_D_nonint() in ./dpq.h : */
#define R_IS_INT(x)	  (fabs((x) - floor((x)+0.5)) <= 1e-7)

double RbMath::lnChoose(double n, double k) {

    double k0 = k;
    k = floor(k + 0.5);
    if (fabs(k - k0) > 1e-7)
    {
        std::ostringstream s;
        s << "'k' (" << k0 << ") must be integer, rounded to " << k;
        throw RbException(s.str());
    }
    if (k < 2)
    {
        if (k <	 0) return RbConstants::Double::neginf;
        if (k == 0) return 0.;
        /* else: k == 1 */
        return log(fabs(n));
    }
    /* else: k >= 2 */
    if (n < 0)
    {
        return lnChoose(-n+ k-1, k);
    }
    else if (R_IS_INT(n))
    {
        if (n < k)
            return RbConstants::Double::neginf;
        /* k <= n :*/
        if (n - k < 2) return lnChoose(n, n-k); /* <- Symmetry */
        /* else: n >= k+2 */
        return RbMath::lfastchoose(n, k);
        }
    /* else non-integer n >= 0 : */
    if (n < k-1)
        {
        int s;
        return RbMath::lfastchoose2(n, k, &s);
        }
    return RbMath::lfastchoose(n, k);
}


double RbMath::lfastchoose(double n, double k) {

    return -log(n + 1.) - RbMath::lnBeta(n - k + 1., k + 1.);
}

/* mathematically the same:
 less stable typically, but useful if n-k+1 < 0 : */
double RbMath::lfastchoose2(double n, double k, int *s_choose) {

    double r = lnGamma_sign(n - k + 1., s_choose);
    return RbMath::lnGamma(n + 1.) - RbMath::lnGamma(k + 1.) - r;
}



/* log factorial ln(n!) */
/*!
 * This function calculates the natural log of the factorial of n using
 * the Stirling approximation.
 *
 * \brief Natural log of the factorial.
 * \param n is the number for the factorial (n!).
 * \return Returns the natural log of the factorial of n.
 * \throws Does not throw an error.
 */
double RbMath::lnFactorial(int n) {

    if (n == 0)
        return 0;
    
    double n1 = n;
    double r  = 1.0 / n1;
    double C0 =  0.918938533204672722;
    double C1 =  1.0/12.0;
    double C3 = -1.0/360.0;
    return (n1 + 0.5) * log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}


std::uint64_t RbMath::stirlingFirst(int n, int k) {

    std::uint64_t r = 0;
    if (n == k)
        return 1;
    if (k == 1)
        return factorial(n-1);
    if (n == 1 && k > 1)
        return 0;
    if (n > 0 && k == 0)
        return 0;
    if (k == n-1)
        return kchoose2(n);
    r = stirlingFirst(n-1, k-1) + ((n-1)*stirlingFirst(n-1, k));
    return r;
}

std::int64_t RbMath::stirlingSecond(int n, int k) {
	
    std::int64_t r = 0;
    if (k > n || k < 0){
	throw RbException("Error: k must be in 0,...,n");
    }
    else{
	if (k == 1 || k == n) {
	    return 1;
	}
	else{
	    std::int64_t v1 = factorial(k);
	    std::int64_t sumv = 0;
	    for (int j=0; j<k+1; j++){
		std::int64_t v2 = 1;
		if ((k-j) % 2 == 1)
		    v2 = -1;
		std::int64_t v3 = v2 * choose(k,j) * pow(j,n);
		sumv += v3;
	    }
	    r = sumv / v1;
	    return r;
	}
    }
    return r;
}

std::int64_t RbMath::bell(int n) {

    std::int64_t r = 0;
    for (int k = 0; k <= n; k++)
    {
        r += stirlingSecond(n, k);
    }
    return r;
}

double RbMath::lnStirlingFirst(int n, int k) {
	
	double r = log(stirlingFirst(n, k));
	return r;
}


std::int64_t RbMath::signedStirlingFirst(int n, int k) {
	int sign = 1;
	if ((n-k) % 2 == 1)
		sign = -1;
	auto v = stirlingFirst(n,k);
	std::cout << v << std::endl;
	return sign * (std::int64_t)v;
}

/*!
 * This function returns the harmonic number of x=> sum_i=1^x(1/i^2)
 *
 * \brief Return sum_i=1^x(1/i^2)
 * \param x The x value
 * \return The harmonic number of x
 */
double RbMath::squaredHarmonicNumber(size_t n)
{
    
    double hm = 0.0;
    for (int i=1; i<=n; i++)
    {
        hm += 1.0/(i*i);
    }
    
    return (hm);
}

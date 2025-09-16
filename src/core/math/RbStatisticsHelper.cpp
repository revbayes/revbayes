//
//  RbStatisticsHelper.cpp
//  RevBayes
//
//  Created by Sebastian Hoehna on 6/7/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <cmath>
#include <cstdint>
#include <optional>
#include <sstream> // IWYU pragma: keep

#include "DistributionNormal.h"
#include "RbStatisticsHelper.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"

using namespace RevBayesCore;

/*!
 * This function is used to calculate the concentration parameter
 * of the DPP from a given number of tables and a given number of 
 * data elements.
 *
 * \brief Subfunction for DPP.
 * \param tables is the given number of tables. 
 * \param num is the given number of data elemnts. 
 * \return Returns an approximate value of alpha. 
 * \throws 
 */
double RbStatistics::Helper::dppConcParamFromNumTables(double tables, double num){
    
    if ( tables < 1.0 ) {
        std::stringstream ss;
        ss << "The Dirichlet Process expects there to be 1 or more clusters. You entered: " << tables << ".";
        
		throw RbException( ss.str() );
    }
	if ( num < 2.0) {
        std::stringstream ss;
		ss << "The Dirichlet Process expects there to be more than 1 data element. You entered: " << num << ".";
		throw RbException( ss.str() );
    }
    if ( tables > num) {
        std::stringstream ss;
        ss << "The Dirichlet Process expects there to be fewer clusters than data elements. You entered: " << tables << " clusters, and " << num << " data elements.";
        throw RbException( ss.str() );
    }
	double a = 0.000001;
	double ea = dppExpectNumTableFromConcParam(a, num);
	bool goUp;
	double target = tables;
	if (target <= 1.0) 
		target = 1.01;
	if (ea < target)
		goUp = true;
	else
		goUp = false;
	double increment = 0.1;
	while ( fabs(ea - target) > 0.000001 ) {
		if (ea < target && goUp == true) {
			a += increment;
		}
		else if (ea > target && goUp == false) {
			a -= increment;
		}
		else if (ea < target && goUp == false) {
			increment /= 2.0;
			goUp = true;
			a += increment;
		}
		else {
			increment /= 2.0;
			goUp = false;
			a -= increment;
		}
		ea = dppExpectNumTableFromConcParam(a, num);
	}
	return a;
}

/*!
 * This function is used to calculate the expected number of tables
 * of the DPP from a given value of alpha and a given number of 
 * data elements.
 *
 * \brief Subfunction for DPP.
 * \param tables is the given number of tables. 
 * \param num is the given number of data elemnts. 
 * \return Returns an approximate value of alpha. 
 * \throws 
 */
double RbStatistics::Helper::dppExpectNumTableFromConcParam(double conp, double num){
    
    if ( conp <= 0.0 )
		throw RbException( "The concentration parameter must be larger than 0.0." );
	if ( num < 2.0)
		throw RbException( "The Dirichlet Process expects there to be more than 1 data element: number > 2." );
	double expectedNum = 0.0;
	for (int i=1; i<=num; i++)
		expectedNum += ( 1.0 / (i - 1.0 + conp) );
	expectedNum *= conp;
	return expectedNum;
}


/**
 * This function is used to calculate the auto-correlation time for MCMC
 * samples, represented here as an arbitrary vector of numeric values.
 *
 * @param values Vector of numeric parameter values
 * @param begin Pointer to start index
 * @param end Pointer to end index
 * @return The effective sample size
 */
double RbStatistics::Helper::getACT(const std::vector<double> values, std::optional<std::int64_t> begin, std::optional<std::int64_t> end)
{
    if (not begin) begin = 0;
    if (not end) end = values.size();
    
    size_t samples = *end - *begin;
    double mean = getRangeMean(values, begin, end);
    double varStat = getVarStat(values, begin, end);
    double gammaStat = 0.0;
    
    for (size_t i = 0; i < samples; i++) {
        double del = values.at( *begin + i ) - mean;
        gammaStat += (del * del);
    }
    
    gammaStat /= ( (double)samples );
    
    double act = varStat / gammaStat;
    return act;
}
    

/**
 * This function is equivalent to the one used in Tracer, and also (except for
 * the choice of the maximum lag value) to convenience:::essTracer().
 *
 * @param values Vector of numeric parameter values
 * @param begin Pointer to start index
 * @param end Pointer to end index
 * @return The effective sample size
 */
double RbStatistics::Helper::getESS(const std::vector<double> values, std::optional<std::int64_t> begin, std::optional<std::int64_t> end)
{
    if (not begin) begin = 0;
    if (not end) end = values.size();
    
    size_t samples = *end - *begin;
    double act = getACT(values, begin, end);
    double ess = samples / act;
    
    return ess;
}


/**
 * This function is used to calculate the mean of an arbitrary contiguous
 * subset of a vector of numeric values: usually, these would represent the
 * samples of numeric parameters from an MCMC simulation.
 *
 * @param values Vector of numeric parameter values
 * @param begin Pointer to start index
 * @param end Pointer to end index
 * @return The effective sample size
 */
double RbStatistics::Helper::getRangeMean(const std::vector<double> values, std::optional<std::int64_t> begin, std::optional<std::int64_t> end)
{
    if (not begin) begin = 0;
    if (not end) end = values.size();
    
    size_t samples = *end - *begin;
    double mean = 0.0;
    
    for (size_t i = 0; i < samples; ++i)
    {
        mean += values.at( *begin + i );
    }
    mean /= samples;
    
    return mean;
}


/**
 * This function is used to calculate the standard error of the mean of an
 * arbitrary contiguous subset of a vector of numeric values: usually, these would
 * represent the samples of numeric parameters from an MCMC simulation.
 *
 * @param values Vector of numeric parameter values
 * @param begin Pointer to start index
 * @param end Pointer to end index
 * @return The effective sample size
 */
double RbStatistics::Helper::getSEM(const std::vector<double> values, std::optional<std::int64_t> begin, std::optional<std::int64_t> end)
{
    if (not begin) begin = 0;
    if (not end) end = values.size();
    
    size_t samples = *end - *begin;
    double varStat = getVarStat(values, begin, end);
    double sem = sqrt(varStat / samples);
    
    return sem;
}
    

double RbStatistics::Helper::getVarStat(const std::vector<double> values, std::optional<std::int64_t> begin, std::optional<std::int64_t> end)
{
    if (not begin) begin = 0;
    if (not end) end = values.size();
    
    size_t MAX_LAG = 1000;
    size_t samples = *end - *begin;
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
    
    // get the mean
    double mean = getRangeMean(values, begin, end);
    
    std::vector<double> gammaStat(maxLag, 0.0);
    double varStat = 0.0;
    
    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values.at( *begin + j ) - mean;
            double del2 = values.at( *begin + j + lag ) - mean;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            } else {
                maxLag = lag;
            }
        }
    }
    
    return varStat;
}


double RbStatistics::Helper::pointChi2(double prob, double v) {
    // Returns z so that Prob{x<z}=prob where x is Chi2 distributed with df
    // = v
    // RATNEST FORTRAN by
    // Best DJ & Roberts DE (1975) The percentage points of the
    // Chi2 distribution. Applied Statistics 24: 385-388. (AS91)
    
    double e = 0.5e-6, aa = 0.6931471805, p = prob;
    double ch, a, q, p1, p2, t, x, b, s1, s2, s3, s4, s5, s6;
    double epsi = .01;
    if ( p < 0.000002 || p > 1 - 0.000002)  {
        epsi = .000001;
    }
    // if (p < 0.000002 || p > 0.999998 || v <= 0) {
    //      throw new IllegalArgumentException("Arguments out of range p" + p + " v " + v);
    //  }
    double g = RbMath::lnGamma(v / 2);
    double xx = v / 2;
    double c = xx - 1;
    if (v < -1.24 * log(p)) {
        ch = pow((p * xx * exp(g + xx * aa)), 1 / xx);
        if (ch - e < 0) {
            return ch;
        }
    } else {
        if (v > 0.32) {
            x = RbStatistics::Normal::quantile(p, 0, 1);
            p1 = 0.222222 / v;
            ch = v * pow((x * sqrt(p1) + 1 - p1), 3.0);
            if (ch > 2.2 * v + 6) {
                ch = -2 * (log(1 - p) - c * log(.5 * ch) + g);
            }
        } else {
            ch = 0.4;
            a = log(1 - p);
            
            
            do {
                q = ch;
                p1 = 1 + ch * (4.67 + ch);
                p2 = ch * (6.73 + ch * (6.66 + ch));
                t = -0.5 + (4.67 + 2 * ch) / p1
                - (6.73 + ch * (13.32 + 3 * ch)) / p2;
                ch -= (1 - exp(a + g + .5 * ch + c * aa) * p2 / p1)
                / t;
            } while (fabs(q / ch - 1) - epsi > 0);
        }
    }
    do {
        q = ch;
        p1 = 0.5 * ch;
        if ((t = RbMath::incompleteGamma(p1, xx)) < 0) {
            throw new RbException("Arguments out of range: t < 0");
        }
        p2 = p - t;
        t = p2 * exp(xx * aa + g + p1 - c * log(ch));
        b = t / ch;
        a = 0.5 * t - b * c;
        
        s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
        s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
        s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
        s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
        s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
        s6 = (120 + c * (346 + 127 * c)) / 5040;
        ch += t
        * (1 + 0.5 * t * s1 - b
           * c
           * (s1 - b
              * (s2 - b
                 * (s3 - b
                    * (s4 - b * (s5 - b * s6))))));
    } while (fabs(q / ch - 1) > e);
    
    return (ch);
}

/*!
 * This function generates a Poisson-distributed random variable using
 * inversion by the chop down method. 
 *
 * \brief Poisson random variables using the chop down method.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 */
int RbStatistics::Helper::poissonInver(double lambda, RandomNumberGenerator& rng) {
    
	const int bound = 130;
	static double p_L_last = -1.0;
	static double p_f0;
	int x;
    
	if (lambda != p_L_last) {
		p_L_last = lambda;
		p_f0 = exp(-lambda);
    } 
    
	while (1) {  
		double r = rng.uniform01();  
		x = 0;  
		double f = p_f0;
		do {
			r -= f;
			if (r <= 0.0) 
				return x;
			x++;
			f *= lambda;
			r *= x;
        } while (x <= bound);
    }
}

/*!
 * This function generates a Poisson-distributed random variable for
 * small values of lambda. The method is a simple calculation of the 
 * probabilities of x = 1 and x = 2. Larger values are ignored.
 *
 * \brief Poisson random variables for small lambda.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 */
int RbStatistics::Helper::poissonLow(double lambda, RandomNumberGenerator& rng) {
    
	double d = sqrt(lambda);
	if (rng.uniform01() >= d) 
		return 0;
	double r = rng.uniform01() * d;
	if (r > lambda * (1.0 - lambda)) 
		return 0;
	if (r > 0.5 * lambda * lambda * (1.0 - lambda)) 
		return 1;
	return 2;
}

/*!
 * This function generates a Poisson-distributed random variable using
 * the ratio-of-uniforms rejectin method. 
 *
 * \brief Poisson random variables using the ratio-of-uniforms method.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 * \see Stadlober, E. 1990. The ratio of uniforms approach for generating
 *      discrete random variates. Journal of Computational and Applied 
 *      Mathematics 31:181-189.
 */
int RbStatistics::Helper::poissonRatioUniforms(double lambda, RandomNumberGenerator& rng) {
    
	static double p_L_last = -1.0;  /* previous L */
	static double p_a;              /* hat center */
	static double p_h;              /* hat width */
	static double p_g;              /* ln(L) */
	static double p_q;              /* value at mode */
	static int p_bound;             /* upper bound */
	int mode;                       /* mode */
	double u;                       /* uniform random */
	double lf;                      /* ln(f(x)) */
	double x;                       /* real sample */
	int k;                          /* integer sample */
    
	if (p_L_last != lambda) {
		p_L_last = lambda;
		p_a = lambda + 0.5;
		mode = (int)lambda;
		p_g  = log(lambda);
		p_q = mode * p_g - RbMath::lnFactorial(mode);
		p_h = sqrt(2.943035529371538573 * (lambda + 0.5)) + 0.8989161620588987408;
		p_bound = (int)(p_a + 6.0 * p_h);
    }
    
	while(1) {
		u = rng.uniform01();
		if (u == 0.0) 
			continue;
		x = p_a + p_h * (rng.uniform01() - 0.5) / u;
		if (x < 0 || x >= p_bound) 
			continue;
		k = (int)(x);
		lf = k * p_g - RbMath::lnFactorial(k) - p_q;
		if (lf >= u * (4.0 - u) - 3.0) 
			break;
		if (u * (u - lf) > 1.0) 
			continue;
		if (2.0 * log(u) <= lf) 
			break;
    }
	return(k);
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double RbStatistics::Helper::rndGamma(double s, RandomNumberGenerator& rng) {

//    double r = rndGamma1(s, rng); // this algorithm is broken
//    double r = rndGamma2(s, rng); // this algorithm is broken
//    double r = rndGamma3(s, rng); // this algorithm is unstable for a < 1
    double r = rndGamma4(s, rng);
    return (r);
    
//	double r = 0.0;
//	if (s <= 0.0)      
//		std::cout << "Gamma parameter less than zero" << std::endl;
//	else if (s < 1.0)  
//		r = rndGamma1(s, rng);
//	else if (s > 1.0)  
//		r = rndGamma2(s, rng);
//	else           
//		r = - log(rng.uniform01());
//	return (r);
    
//    // Adapted code from R (rgamma)
//    /* Constants : */
//    const static double sqrt32 = 5.656854;
//    const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
//    
//    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
//     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
//     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
//     */
//    const static double q1 = 0.04166669;
//    const static double q2 = 0.02083148;
//    const static double q3 = 0.00801191;
//    const static double q4 = 0.00144121;
//    const static double q5 = -7.388e-5;
//    const static double q6 = 2.4511e-4;
//    const static double q7 = 2.424e-4;
//    
//    const static double a1 = 0.3333333;
//    const static double a2 = -0.250003;
//    const static double a3 = 0.2000062;
//    const static double a4 = -0.1662921;
//    const static double a5 = 0.1423657;
//    const static double a6 = -0.1367177;
//    const static double a7 = 0.1233795;
//    
//    /* State variables [FIXME for threading!] :*/
//    static double aa = 0.;
//    static double aaa = 0.;
//    static double s, s2, d;    /* no. 1 (step 1) */
//    static double q0, b, si, c;/* no. 2 (step 4) */
//    
//    double e, p, q, r, t, u, v, w, x, ret_val;
//    
//    if (!RbMath::isFinite(a) || a < 0.0 ) {
//        throw RbException("Infinite parameters for rgamma.");
//    }
//    
//    if (a < 1.) { /* GS algorithm for parameters a < 1 */
//        if (a == 0)
//            return 0.;
//        e = 1.0 + exp_m1 * a;
//        while ( true ) {
//            p = e * rng.uniform01();
//            if (p >= 1.0) {
//                x = -log((e - p) / a);
//                if (exp( rng.uniform01() ) >= (1.0 - a) * log(x))
//                    break;
//            } else {
//                x = exp(log(p) / a);
//                if (exp( rng.uniform01() ) >= x)
//                    break;
//            }
//        }
//        return x;
//    }
//    
//    /* --- a >= 1 : GD algorithm --- */
//    
//    /* Step 1: Recalculations of s2, s, d if a has changed */
//    if (a != aa) {
//        aa = a;
//        s2 = a - 0.5;
//        s = sqrt(s2);
//        d = sqrt32 - s * 12.0;
//    }
//    /* Step 2: t = standard normal deviate,
//     x = (s,1/2) -normal deviate. */
//    
//    /* immediate acceptance (i) */
//    t = RbStatistics::Normal::rv(rng);
//    x = s + 0.5 * t;
//    ret_val = x * x;
//    if (t >= 0.0)
//        return ret_val;
//    
//    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
//    u = rng.uniform01();
//    if (d * u <= t * t * t)
//        return ret_val;
//    
//    /* Step 4: recalculations of q0, b, si, c if necessary */
//    
//    if (a != aaa) {
//        aaa = a;
//        r = 1.0 / a;
//        q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
//               + q2) * r + q1) * r;
//        
//        /* Approximation depending on size of parameter a */
//        /* The constants in the expressions for b, si and c */
//        /* were established by numerical experiments */
//        
//        if (a <= 3.686) {
//            b = 0.463 + s + 0.178 * s2;
//            si = 1.235;
//            c = 0.195 / s - 0.079 + 0.16 * s;
//        } else if (a <= 13.022) {
//            b = 1.654 + 0.0076 * s2;
//            si = 1.68 / s + 0.275;
//            c = 0.062 / s + 0.024;
//        } else {
//            b = 1.77;
//            si = 0.75;
//            c = 0.1515 / s;
//        }
//    }
//    /* Step 5: no quotient test if x not positive */
//    
//    if (x > 0.0) {
//        /* Step 6: calculation of v and quotient q */
//        v = t / (s + s);
//        if (fabs(v) <= 0.25)
//            q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
//                                      + a3) * v + a2) * v + a1) * v;
//        else
//            q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
//        
//        
//        /* Step 7: quotient acceptance (q) */
//        if (log(1.0 - u) <= q)
//            return ret_val;
//    }
//    
//    while (true) {
//        /* Step 8: e = standard exponential deviate
//         *	u =  0,1 -uniform deviate
//         *	t = (b,si)-double exponential (laplace) sample */
//        e = exp( rng.uniform01() );
//        u = rng.uniform01();
//        u = u + u - 1.0;
//        if (u < 0.0)
//            t = b - si * e;
//        else
//            t = b + si * e;
//        /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
//        if (t >= -0.71874483771719) {
//            /* Step 10:	 calculation of v and quotient q */
//            v = t / (s + s);
//            if (fabs(v) <= 0.25)
//                q = q0 + 0.5 * t * t *
//                ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
//                  + a2) * v + a1) * v;
//            else
//                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
//            /* Step 11:	 hat acceptance (h) */
//            /* (if q not positive go to step 8) */
//            if (q > 0.0) {
//                w = expm1(q);
//                /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
//                /* if t is rejected sample again at step 8 */
//                if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
//                    break;
//            }
//        }
//    } /* repeat .. until  `t' is accepted */
//    x = s + 0.5 * t;
//    return x * x;

    
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double RbStatistics::Helper::rndGamma1(double s, RandomNumberGenerator& rng)
{
    
    double r, x = 0.0, small = 1e-37, w;
    static double   a, p, uf, ss = 10.0, d;

    if (s != ss) {
        a  = 1.0 - s;
        p  = a / (a + s * exp(-a));
        uf = p * pow(small / a, s);
        d  = a * log(a);
        ss = s;
    }
    for (;;) {
        r = rng.uniform01();
        if (r > p)
        {
            x = a - log((1.0 - r) / (1.0 - p));
            w = a * log(x) - d;
        }
        else if (r>uf)
        {
            x = a * pow(r / p, 1.0 / s);
            w = x;
        }
        else
            return (0.0);
        r = rng.uniform01();
        if (1.0 - r <= w && r > 0.0)
            if (r*(w + 1.0) >= 1.0 || -log(r) <= w)  
                continue;
        break;
    }
    
    return (x);
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double RbStatistics::Helper::rndGamma2(double s, RandomNumberGenerator& rng)
{
    
    double              r, d, f, g, x;
    static double       b, h, ss = 0.0;

    if (s != ss) {
        b  = s - 1.0;
        h  = sqrt(3.0 * s - 0.75);
        ss = s;
    }
    for (;;) {
        r = rng.uniform01();
        g = r - r * r;
        f = (r - 0.5) * h / sqrt(g);
        x = b + f;
        if (x <= 0.0)
            continue;
        r = rng.uniform01();
        d = 64.0 * r * r * g * g * g;
        if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))  
            break;
    }
    
    return x;
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param a is the shape parameter of the gamma.
 * \return Returns a gamma-distributed random variable.
 * \throws Does not throw an error.
 */
double RbStatistics::Helper::rndGamma3(double a, RandomNumberGenerator& rng) {
    
    // Adapted code from R (rgamma)
    /* Constants : */
    const static double sqrt32 = 5.656854;
    const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
    
    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
     */
    const static double q1 = 0.04166669;
    const static double q2 = 0.02083148;
    const static double q3 = 0.00801191;
    const static double q4 = 0.00144121;
    const static double q5 = -7.388e-5;
    const static double q6 = 2.4511e-4;
    const static double q7 = 2.424e-4;
    
    const static double a1 = 0.3333333;
    const static double a2 = -0.250003;
    const static double a3 = 0.2000062;
    const static double a4 = -0.1662921;
    const static double a5 = 0.1423657;
    const static double a6 = -0.1367177;
    const static double a7 = 0.1233795;
    
    /* State variables [FIXME for threading!] :*/
    static double aa = 0.;
    static double aaa = 0.;
    static double s, s2, d;    /* no. 1 (step 1) */
    static double q0, b, si, c;/* no. 2 (step 4) */
    
    double e, p, q, r, t, u, v, w, x, ret_val;
    
    if (!RbMath::isFinite(a) || a < 0.0 ) {
        throw RbException("NaN produced in rgamma.");
    }
    
    if (a < 1.) { /* GS algorithm for parameters a < 1 */
        if (a == 0)
            return 0.;
        e = 1.0 + exp_m1 * a;
        while ( true ) {
            p = e * rng.uniform01();
            if (p >= 1.0) {
                x = -log((e - p) / a);
                if (-log( 1.0 - rng.uniform01() ) >= (1.0 - a) * log(x))
                    break;
            } else {
                x = exp(log(p) / a);
                if (-log( 1.0 - rng.uniform01() ) >= x)
                    break;
            }
        }
        return x;
    }
    
    /* --- a >= 1 : GD algorithm --- */
    
    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a != aa) {
        aa = a;
        s2 = a - 0.5;
        s = sqrt(s2);
        d = sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
     x = (s,1/2) -normal deviate. */
    
    /* immediate acceptance (i) */
    t = RbStatistics::Normal::rv(rng);
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0)
        return ret_val;
    
    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = rng.uniform01();
    if (d * u <= t * t * t)
        return ret_val;
    
    /* Step 4: recalculations of q0, b, si, c if necessary */
    
    if (a != aaa) {
        aaa = a;
        r = 1.0 / a;
        q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
               + q2) * r + q1) * r;
        
        /* Approximation depending on size of parameter a */
        /* The constants in the expressions for b, si and c */
        /* were established by numerical experiments */
        
        if (a <= 3.686) {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
        } else if (a <= 13.022) {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
        } else {
            b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
        }
    }
    /* Step 5: no quotient test if x not positive */
    
    if (x > 0.0) {
        /* Step 6: calculation of v and quotient q */
        v = t / (s + s);
        if (fabs(v) <= 0.25)
            q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                      + a3) * v + a2) * v + a1) * v;
        else
            q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        
        
        /* Step 7: quotient acceptance (q) */
        if (log(1.0 - u) <= q)
            return ret_val;
    }
    
    while (true) {
        /* Step 8: e = standard exponential deviate
         *	u =  0,1 -uniform deviate
         *	t = (b,si)-double exponential (laplace) sample */
        e = -log( 1.0 - rng.uniform01() );
        u = rng.uniform01();
        u = u + u - 1.0;
        if (u < 0.0)
            t = b - si * e;
        else
            t = b + si * e;
        /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
        if (t >= -0.71874483771719) {
            /* Step 10:	 calculation of v and quotient q */
            v = t / (s + s);
            if (fabs(v) <= 0.25)
                q = q0 + 0.5 * t * t *
                ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
                  + a2) * v + a1) * v;
            else
                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
            /* Step 11:	 hat acceptance (h) */
            /* (if q not positive go to step 8) */
            if (q > 0.0) {
                w = expm1(q);
                /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                /* if t is rejected sample again at step 8 */
                if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
                    break;
            }
        }
    } /* repeat .. until  `t' is accepted */
    x = s + 0.5 * t;
    return x * x;

}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma.
 * \return Returns a gamma-distributed random variable.
 * \throws Does not throw an error.
 */
double RbStatistics::Helper::rndGamma4(double s, RandomNumberGenerator& rng) {
    
    // This algorithm is based on:
    // G. Marsaglia and W. Tsang. A simple method for generating gamma variables. ACM Transactions on Mathematical Software, 26(3):363-372, 2000
    
    double alpha        = s;
    bool   alpha_scaled = false;
    
    // check that alpha is finite
    if ( !RbMath::isFinite(s) || s < 0.0 )
    {
        throw RbException("NaN produced in rgamma.");
    }
    
    // check if alpha is less than 1
    if ( s < 1.0 )
    {
        alpha = s + 1;
        alpha_scaled = true;
    }
    
    double d = alpha - 1.0 / 3.0;
    double c = 1 / sqrt(9 * d);
    
    double x;
    while(true) {
        
        double z = RbStatistics::Normal::rv(rng);
        while(z <= -1 / c)
        {
            z = RbStatistics::Normal::rv(rng);
        }
        double v = pow(1 + c * z, 3);
        
        double u = rng.uniform01();
        if ( u < 1 - 0.0331 * pow(z, 4) )
        {
            x = d * v;
            break;
        }
        
        double log_u = log(u);
        double rhs   = 0.5 * pow(z, 2) + d * (1 - v + log(v));
        if ( log_u < rhs )
        {
            x = d * v;
            break;
        }
        
    }
    
    // check if we need to transform the RV
    if ( alpha_scaled == true ) {
        double u = rng.uniform01();
        x *= pow(u, 1 / s);
        
//        // make sure that x cannot be exactly 0
//        if (x == 0) {
//            x = 1e-300;
//        }
        
    }
    
    return x;
    
}








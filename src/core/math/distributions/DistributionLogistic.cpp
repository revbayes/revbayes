/**
 * @file DistributionLogistic
 * This file contains the functions of the logistic distribution.
 *
 * @brief Implementation of the logistic distribution.
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

#include <cmath>
#include "DistributionLogistic.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbStatisticsHelper.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

double RbStatistics::Logistic::pdf(double x){

    return RbStatistics::Logistic::pdf(0,1, x);
}

/*!
 * This function calculates the probability density
 * for a logistic-distributed random variable.
 *
 * \brief logistic probability density.
 * \param location
 * \param scale
 * \param x
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
//double RbStatistics::Dirichlet::pdf(const std::vector<double> &a, const std::vector<double> &z) {



double RbStatistics::Logistic::pdf(double location, double scale, double x){

    double e, f;
    if (scale <= 0.0){
        std::ostringstream s;
        s << "Fatal error in Logistic PDF";
        throw RbException(s.str());
    }
////	ML_ERR_return_NAN;

    x = fabs((x - location) / scale);
    e = exp(-x);
    f = 1.0 + e;
    double pdf =  e / (scale * f * f);
    return pdf;

}

/*!
 * This function calculates the natural log of the probability density
 * for a Logistic-distributed random variable.
 *
 * \brief Natural log of Logistic probability density.
 * \param location: location parameter for logistic distribution
 * \param scale:    scale parameter for logistic distribution
 * \param x:        random variable for logistic distirbution
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */

double RbStatistics::Logistic::lnPdf(double location, double scale, double x){


    double e, f;
    if (scale <= 0.0){
        std::ostringstream s;
        s << "Fatal error in Logistic PDF";
        throw RbException(s.str());
    }

    x = fabs((x - location) / scale);
    e = exp(-x);
    f = 1.0 + e;
    double lnP =  -(x + log(scale * f * f)) ;


	return lnP;
}


double RbStatistics::Logistic::cdf(double location, double scale, double x)
{


    if (scale <= 0.0)
    {
        std::ostringstream s;
        s << "Error in Logistic CDF";
        throw RbException(s.str());
    }

    x = (x - location) / scale;

    if (RbMath::isNan(x))
    {
        std::ostringstream s;
        s << "Error in Logistic CDF";
        throw RbException(s.str());
    }

    x = exp(-x );
    return (1/(1 + x));

}

double RbStatistics::Logistic::quantile( double p)
{

    return log((p / (1 - p)));

}

double RbStatistics::Logistic::quantile(double location, double scale, double p)
{

    if (scale < 0.0){
        std::ostringstream s;
        s << "error in Logistic quatile function";
        throw RbException(s.str());
    }
    if (scale == 0)	{
        return location;
    }

	p = Logistic::quantile(p);

    return location + scale * p;

}
/*!
 * This function generates a logistic-distributed random variable.
 *
 * \brief logistic random variable.
 * \param a is a reference to a vector of doubles containing the logistic parameters.
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the logistic random variable.
 * \throws Does not throw an error.
 */

double RbStatistics::Logistic::rv(double location, double scale, RandomNumberGenerator& rng) {

    if (scale == 0){
        return (location);

    }
    else{

        double u = rng.uniform01();
        return location + scale * log(u / (1. - u));
    }
}

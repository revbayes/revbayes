//
//  DistributionWishart.cpp
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-28.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//



#include <cstddef>
#include <cmath>
#include <vector>

#include "DistributionMultivariateNormal.h"
#include "DistributionWishart.h"
#include "RbConstants.h"
#include "DistributionNormal.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a Wishart-distributed random variable.
 *
 * \brief Wishart probability density.
 * \param omega0 is a reference to the scale matrix
 * \param df is the number of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::Wishart::pdf(const MatrixReal &omega0, size_t df, const MatrixReal &z) {
	
    return exp(lnPdf(omega0,df,z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a Wishart-distributed random variable.
 *
 * \brief Natural log of Wishart probability density.
 * \param omega0 is a reference to the scale matrix (a precision matrix)
 * \param df is the numger of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Wishart::lnPdf(const MatrixReal &omega0, size_t df, const MatrixReal &z)
{
    
//    omega0.update();
//    z.update();
    
    if (! z.isPositiveDefinite() )
    {
         return RbConstants::Double::neginf;
    }
    
    double ret = 0;
    ret -= 0.5 * df * omega0.getLogDet();
    ret += 0.5 * (df - omega0.getDim() - 1) * z.getLogDet();
    
    double trace = 0;
    
    MatrixReal inv = omega0.computeInverse();
    
    for (size_t i=0; i<omega0.getDim(); i++)
    {
        for (size_t j=0; j<omega0.getDim(); j++)
        {
            trace += inv[i][j] * z[j][i];
        }
    }
    
    ret -= 0.5 * trace;
    
    return ret;

}



/*!
 * This function generates a Wishart-distributed random variable.
 *
 * \brief Wishart random variable.
 * \param omega0 is a reference to the scale matrix (a precision matrix)
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the Wishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::Wishart::rv(const MatrixReal &omega0, size_t df, RandomNumberGenerator& rng)
{

//    omega0.update();
    
    size_t dim = omega0.getDim();
    MatrixReal z(dim);
    std::vector<double> mean = std::vector<double>(dim, 0.0);
    
    for (size_t k=0; k<df; k++)
    {
        std::vector<double> tmp = RbStatistics::MultivariateNormal::rvCovariance(mean, omega0, rng);
        for (size_t i=0; i<dim; i++)
        {
            for (size_t j=0; j<dim; j++)
            {
                z[i][j] += tmp[i] * tmp[j];
            }
        }
    }

    return z;
}


/*!
 * This function calculates the probability density
 * for a Wishart-distributed random variable.
 *
 * \brief Wishart probability density.
 * \param kappa is a reference to a double; the scale matrix is then omega0 = kappa * I
 * \param df is the number of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::Wishart::pdf(double kappa, size_t df, const MatrixReal &z) {
	
    return exp(lnPdf(kappa,df,z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a Wishart-distributed random variable.
 *
 * \brief Natural log of Wishart probability density.
 * \param kappa is a reference to a double; the scale matrix is then omega0 = kappa * I
 * \param df is the numger of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */

// this log density is only up to a normalization factor that *does* depend on df
// df is therefore assumed to be constant throughout

double RbStatistics::Wishart::lnPdf(double kappa, size_t df, const MatrixReal &z)
{

    
//    z.update();
    int dim = int(z.getDim());
    
    if ( !z.isPositiveDefinite() )
    {
        return RbConstants::Double::neginf;
    }
    
    double ret = 0;
    ret -= 0.5 * df * dim * log(kappa);
    double logDet = z.getLogDet();
    ret += 0.5 * (int(df) - dim - 1) * logDet;
    
    double trace = 0;
    
    for (size_t i=0; i<dim; i++)
    {
        trace += z[i][i];
    }
    trace /= kappa;
    
    ret -= 0.5 * trace;
    
    return ret;

}

/*!
 * This function generates a Wishart-distributed random variable.
 *
 * \brief Wishart random variable.
 * \param kappa is a reference to a double; the scale matrix is then omega0 = kappa * I
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the Wishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::Wishart::rv(double kappa, size_t dim, size_t df, RandomNumberGenerator& rng)
{
        
    MatrixReal z(dim);
    std::vector<double> tmp(dim);

    double sk = sqrt(kappa);
    for (size_t k=0; k<df; k++)
    {
        
        for (size_t i=0; i<dim; i++)
        {
            tmp[i] = RbStatistics::Normal::rv(0, sk, rng);
        }

        for (size_t i=0; i<dim; i++)
        {
            for (size_t j=0; j<dim; j++)
            {
                z[i][j] += tmp[i] * tmp[j];
            }
        }
    }
    
    return z;
}


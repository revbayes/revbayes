/* 
 * File:   DistributionInverseInverseWishart.cpp
 * Author: nl
 * 
 * Created on 15 juillet 2014, 10:13
 */


#include "DistributionInverseWishart.h"
#include "DistributionMultivariateNormal.h"
#include "RbException.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "RbStatisticsHelper.h"
#include "DistributionNormal.h"

#include <cmath>

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart probability density.
 * \param sigma0 is a reference to the scale matrix
 * \param df is the number of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::InverseWishart::pdf(const MatrixReal &sigma0, size_t df, const MatrixReal &z) {
	
    return exp(lnPdf(sigma0,df,z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief Natural log of InverseWishart probability density.
 * \param sigma0 is a reference to the scale matrix (a precision matrix)
 * \param df is the numger of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseWishart::lnPdf(const MatrixReal &sigma0, size_t df, const MatrixReal &z)
{
    
    if ( !z.isPositiveDefinite() )
    {
        return RbConstants::Double::neginf;
    }
    
    double ret = 0;
    ret += 0.5 * df * sigma0.getLogDet();
    ret -= 0.5 * (df + sigma0.getDim() + 1) * z.getLogDet();
    
    double trace = 0;
    
    MatrixReal invz = z.computeInverse();
    
    for (size_t i=0; i<sigma0.getDim(); i++)
    {
        for (size_t j=0; j<sigma0.getDim(); j++)
        {
            trace += sigma0[i][j] * invz[j][i];
        }
    }
    
    ret -= 0.5 * trace;
    
    return ret;

}



/*!
 * This function generates a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart random variable.
 * \param sigma0 is a reference to the scale matrix (a precision matrix)
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the InverseWishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::InverseWishart::rv(const MatrixReal &sigma0, size_t df, RandomNumberGenerator& rng)
{
    
    size_t dim = sigma0.getDim();
        
    MatrixReal z = MatrixReal(dim);
    std::vector<double> mean = std::vector<double>(dim, 0.0);
    
    for (size_t k=0; k<df; k++)
    {
        std::vector<double> tmp = RbStatistics::MultivariateNormal::rvPrecision(mean, sigma0, rng);
        for (size_t i=0; i<dim; i++)
        {
            for (size_t j=0; j<dim; j++)
            {
                z[i][j] += tmp[i] * tmp[j];
            }
        }
    }

    z.setCholesky(true);
    return z.computeInverse();
//    return z;
}

/*!
 * This function generates a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart random variable.
 * \param sigma0 is a reference to the scale matrix (a covariance matrix)
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the InverseWishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::InverseWishart::rvCovariance(const MatrixReal &sigma0, size_t df, RandomNumberGenerator& rng)
{
    
    size_t dim = sigma0.getDim();
    
    MatrixReal z = MatrixReal(dim);
    std::vector<double> mean = std::vector<double>(dim, 0.0);
    
    for (size_t k=0; k<df; k++)
    {
        std::vector<double> tmp = RbStatistics::MultivariateNormal::rvCovariance(mean, sigma0, rng);
        for (size_t i=0; i<dim; i++)
        {
            for (size_t j=0; j<dim; j++)
            {
                z[i][j] += tmp[i] * tmp[j];
            }
        }
    }
    
    z.setCholesky(true);
    return z.computeInverse();
    //    return z;
}


/*!
 * This function calculates the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart probability density.
 * \param kappa is a reference to a vector of double; the scale matrix is then sigma0 = Diagonal(kappa)
 * \param df is the number of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::InverseWishart::pdf(const std::vector<double>& kappa, size_t df, const MatrixReal &z)
{
	
    return exp(lnPdf(kappa,df,z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief Natural log of InverseWishart probability density.
 * \param kappa is a reference to a vector of double; the scale matrix is then sigma0 = Diagonal(kappa)
 * \param df is the numger of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */

// this log density is only up to a normalization factor that *does* depend on df
// df is therefore assumed to be constant throughout

double RbStatistics::InverseWishart::lnPdf(const std::vector<double>& kappa, size_t df, const MatrixReal &z)
{

    
    size_t dim = z.getDim();
    
    if ( z.isPositiveDefinite() == false )
    {
        return RbConstants::Double::neginf;
    }
    
    double ret = 0;
    for (size_t i=0; i<dim; i++)
    {
        ret += 0.5 * df * log(kappa[i]);
    }
    ret -= 0.5 * (df + dim + 1) * z.getLogDet();
   
    double trace = 0;
    
    MatrixReal invz = z.computeInverse();
    for (size_t i=0; i<dim; i++)
    {
        trace += kappa[i] * invz[i][i];
    }
    
    ret -= 0.5 * trace;
    
    return ret;

}

/*!
 * This function generates a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart random variable.
 * \param kappa is a reference to a vector of double; the scale matrix is then sigma0 = Diagonal(kappa)
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the InverseWishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::InverseWishart::rv(const std::vector<double>& kappa, size_t df, RandomNumberGenerator& rng) {
        
    size_t dim = kappa.size();
    MatrixReal z(dim);
    std::vector<double> tmp(dim);

    for (size_t k=0; k<df; k++)
    {

        for (size_t i=0; i<dim; i++)
        {
            double sk = 1.0 / sqrt(kappa[i]);
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
    
    return z.computeInverse();
}



/*!
 * This function calculates the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart probability density.
 * \param kappa is a reference to a double; the scale matrix is then sigma0 = kappa * I
 * \param df is the number of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::InverseWishart::pdf(double kappa, size_t df, const MatrixReal &z)
{
	
    return exp(lnPdf(kappa,df,z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a InverseWishart-distributed random variable.
 *
 * \brief Natural log of InverseWishart probability density.
 * \param kappa is a reference to a double; the scale matrix is then sigma0 = kappa * I
 * \param df is the numger of degrees of freedom
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */

// this log density is only up to a normalization factor that *does* depend on df
// df is therefore assumed to be constant throughout

double RbStatistics::InverseWishart::lnPdf(double kappa, size_t df, const MatrixReal &z)
{

    
    size_t dim = z.getDim();
    
    if ( !z.isPositiveDefinite() )
    {
        return RbConstants::Double::neginf;
    }
    
    double ret = 0;
    ret += 0.5 * df * dim * log(kappa);
    ret -= 0.5 * (df + dim + 1) * z.getLogDet();
   
    double trace = 0;
    
    MatrixReal invz = z.computeInverse();
    for (size_t i=0; i<dim; i++)
    {
        trace += kappa * invz[i][i];
    }
    
    ret -= 0.5 * trace;
    
    return ret;

}

/*!
 * This function generates a InverseWishart-distributed random variable.
 *
 * \brief InverseWishart random variable.
 * \param kappa is a reference to a double; the scale matrix is then sigma0 = kappa * I
 * \param df is the numger of degrees of freedom
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the InverseWishart random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::InverseWishart::rv(double kappa, size_t dim, size_t df, RandomNumberGenerator& rng)
{

    MatrixReal z(dim);
    std::vector<double> tmp = std::vector<double>(dim,0.0);
    
    double sk = 1.0 / sqrt(kappa);
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
    
    return z.computeInverse();
}


//
//  DistributionMultivariateNormal.cpp
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-28.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//


#include "DistributionMultivariateNormal.h"

#include <cstddef>
#include <cmath>

#include "DistributionNormal.h"
#include "CholeskyDecomposition.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a MultivariateNormal-distributed random variable.
 *
 * \brief MultivariateNormal probability density.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param sigma is a reference to a precision matrix containing the covariance
 * \param z is a reference to a vector of doubles containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::MultivariateNormal::pdfCovariance(const std::vector<double>& mu, const MatrixReal& sigma, const std::vector<double> &x, double scale)
{
	
    return exp(lnPdfCovariance(mu,sigma,x,scale));
}


/*!
 * This function calculates the natural log of the probability density
 * for a MultivariateNormal-distributed random variable.
 *
 * \brief Natural log of MultivariateNormal probability density.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param omega0 is a reference to a precision matrix containing the covariance
 * \param x is a reference to a vector of doubles containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::MultivariateNormal::lnPdfCovariance(const std::vector<double>& mu, const MatrixReal& sigma, const std::vector<double> &x, double scale)
{
    // we compute the precision matrix, which is the inverse of the covariance matrix
    // and then simply call the lnPDF for the precision matrix.
    // This simplifies the coding.
    sigma.setCholesky(true);
    MatrixReal omega = sigma.computeInverse();

    return lnPdfPrecision(mu, omega, x, scale);
}


/*!
 * This function calculates the natural log of the probability density
 * for a MultivariateNormal-distributed random variable.
 *
 * \brief Natural log of MultivariateNormal probability density.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param omega0 is a reference to a precision matrix containing the covariance
 * \param x is a reference to a vector of doubles containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::MultivariateNormal::lnPdfCovariance(const std::vector<double>& mu, const MatrixReal& sigma, const std::vector<double> &x, const std::vector<double> & scale)
{
    // we compute the precision matrix, which is the inverse of the covariance matrix
    // and then simply call the lnPDF for the precision matrix.
    // This simplifies the coding.
    MatrixReal sigma_scaled = sigma;
    for (size_t i=0; i<scale.size(); ++i)
    {
        for (size_t j=0; j<scale.size(); ++j)
        {
            sigma_scaled[i][j] = sigma[i][j]*scale[i]*scale[j];
        }
    }
    sigma_scaled.setCholesky(true);
    MatrixReal omega = sigma_scaled.computeInverse();
    
    return lnPdfPrecision(mu, omega, x, 1.0);
}

/*!
 * This function generates a MultivariateNormal-distributed random variable.
 *
 * \brief MultivariateNormal random variable.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param sigma is a reference to a precision matrix containing the covariance
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the MultivariateNormal random variable.
 * \throws Does not throw an error.
 */

std::vector<double> RbStatistics::MultivariateNormal::rvCovariance(const std::vector<double>& mu, const MatrixReal& sigma, RandomNumberGenerator& rng, double scale)
{
    
    sigma.setCholesky(true);
    double sqrtScale = sqrt(scale);
    size_t dim = sigma.getDim();
    
    MatrixReal W(dim, 1, 0.0);
    for (size_t i = 0; i < dim; ++i)
    {
        W[i][0] = RbStatistics::Normal::rv(0, sqrtScale, rng);
    }
    const CholeskyDecomposition& cd = sigma.getCholeskyDecomposition();
    const MatrixReal L = cd.getLowerCholeskyFactor();
    
    MatrixReal V = L * W;
    std::vector<double> v = std::vector<double>(dim, 0.0);
    for (size_t i = 0; i < dim; ++i)
    {
        v[i] = mu[i] + V[i][0];
    }
    
    return v;

//    double sqrtScale = sqrt(scale);
//    size_t dim = sigma.getDim();
//    std::vector<double> v = std::vector<double>(dim, 0.0);
//    std::vector<double> w = std::vector<double>(dim, 0.0);
//
//    // the eigen system of the covariance matrix
//    const EigenSystem& eigensystem = sigma.getEigenSystem();
//    
//    // get the eigenvalues of the *covariance* matrix
//    const std::vector<double>& eigen = eigensystem.getRealEigenvalues();
//    
//    // draw the normal variate in eigen basis
//    for (size_t i=0; i<dim; i++)
//    {
////        w[i] = RbStatistics::Normal::rv(0, sqrtScale, rng);
//        
//        if ( eigen[i] < 0.0 )
//        {
//            throw RbException("Cannot draw random value of multivariate normal distribution because eigenvalues of the covariance matrix are negative.");
//        }
//        
//        w[i] = RbStatistics::Normal::rv(0, sqrtScale * sqrt(eigen[i]), rng);
//    }
//    
//    // get the eigenvector
//    const MatrixReal& eigenvect = eigensystem.getEigenvectors();
//    
//    // change basis
//    for (size_t i = 0; i < dim; ++i)
//    {
//        double tmp = 0;
//        for (size_t j = 0; j < dim; ++j)
//        {
//            tmp += eigenvect[i][j] * w[j];
//        }
//        v[i] = tmp + mu[i];
//        
//    }
//    
//    return v;
    
}

/*!
 * This function calculates the probability density
 * for a MultivariateNormal-distributed random variable.
 *
 * \brief MultivariateNormal probability density.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param omega is a reference to a precision matrix containing the precision
 * \param x is a reference to a vector of doubles containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::MultivariateNormal::pdfPrecision(const std::vector<double>& mu, const MatrixReal& omega, const std::vector<double> &x, double scale)
{
	
    return exp(lnPdfPrecision(mu,omega,x,scale));
}


/*!
 * This function calculates the natural log of the probability density
 * for a MultivariateNormal-distributed random variable.
 *
 * \brief Natural log of MultivariateNormal probability density.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param omega is a reference to a precision matrix containing the precision
 * \param x is a reference to a vector of doubles containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::MultivariateNormal::lnPdfPrecision(const std::vector<double>& mu, const MatrixReal& omega, const std::vector<double> &x, double scale)
{
    
    double logNormalize = -0.5 * log( RbConstants::TwoPI );
    
    double logDet = omega.getLogDet();
    if ( RbMath::isAComputableNumber(logDet) == false )
    {
        return logDet;
    }
   
    // check positive semidefiniteness
    if ( omega.isUsingCholesky() == true )
    {
        if ( omega.getCholeskyDecomposition().checkPositiveSemidefinite() == false )
        {
            return RbConstants::Double::neginf;
        }
    }
    else
    {
        CholeskyDecomposition c = CholeskyDecomposition(&omega);
        if ( c.checkPositiveSemidefinite() == false )
        {
            return RbConstants::Double::neginf;
        }
    }

    size_t dim = x.size();
    std::vector<double> tmp = std::vector<double>(dim,0.0);
    
    double s2 = 0;
    for (size_t i=0; i<dim; i++)
    {
        double tmp = 0;
        for (size_t j=0; j<dim; j++)
        {
            tmp += omega[i][j] * (x[j] - mu[j]);
        }
        s2 += (x[i] - mu[i]) * tmp;
    }
    
    double lnProb = dim * logNormalize + 0.5 * (logDet - dim * log(scale) - s2 / scale);

    return lnProb;
    
}


/*!
 * This function generates a MultivariateNormal-distributed random variable.
 *
 * \brief MultivariateNormal random variable.
 * \param mu is a reference to a vector of doubles containing the mean
 * \param omega is a reference to a precision matrix containing the precision
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the MultivariateNormal random variable.
 * \throws Does not throw an error.
 */
std::vector<double> RbStatistics::MultivariateNormal::rvPrecision(const std::vector<double>& mu, const MatrixReal& omega, RandomNumberGenerator& rng, double scale)
{
 
    omega.setCholesky(true);
    return rvCovariance(mu, omega.computeInverse(), rng, scale);
    
//    double sqrtScale = sqrt(scale);
//    size_t dim = omega.getDim();
//    std::vector<double> v = std::vector<double>(dim, 0.0);
//    std::vector<double> w = std::vector<double>(dim, 0.0);
//
//    const EigenSystem &eigensystem = omega.getEigenSystem();
//    
//    // get the eigenvalues of the *precision* matrix
//    const std::vector<double>& eigen = eigensystem.getRealEigenvalues();
//    
//    // draw the normal variate in eigen basis
//    for (size_t i=0; i<dim; i++)
//    {
////        w[i] = RbStatistics::Normal::rv(0, sqrtScale, rng);
//        if ( eigen[i] < 0.0 )
//        {
//            throw RbException("Cannot draw random value of multivariate normal distribution because eigenvalues of the covariance matrix are negative.");
//        }
//
//        w[i] = RbStatistics::Normal::rv(0, sqrtScale / sqrt(eigen[i]), rng);
//    }
//    
//    // get the eigenvector
//    const MatrixReal& eigenvect = eigensystem.getEigenvectors();
//    
//    // change basis
//    for (size_t i=0; i<dim; i++)
//    {
//        double tmp = 0;
//        for (size_t j=0; j<dim; j++)
//        {
//            tmp += eigenvect[i][j] * w[j];
//        }
//        v[i] = tmp + mu[i];
//    }
//
//    
//    return v;
}

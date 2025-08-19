#include <cstddef>
#include <cmath>

#include "CholeskyDecomposition.h"
#include "DistributionBeta.h"
#include "DistributionLKJ.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a LKJ-distributed random variable.
 *
 * \brief LKJ probability density.
 * \param eta is the parameter of the LKJ distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::LKJ::pdf(double eta, const MatrixReal &z) {
	
    return exp(lnPdf(eta, z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a LKJ-distributed random variable.
 *
 * \brief Natural log of LKJ probability density.
 * \param eta is the parameter of the LKJ distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::LKJ::lnPdf(double eta, const MatrixReal &z)
{
    
    size_t dim = z.getNumberOfRows();
    
    // check that all rows are 1 and all colums are between -1 and 1
    for (size_t r = 0; r < dim; ++r)
    {
        
        if (z[r][r] != 1.0) {
            return RbConstants::Double::neginf;
        }
        
        for (size_t c = r; c < dim; ++c)
        {
            if (c > r)
            {
                if ( z[r][c] > 1.0 or z[r][c] < -1.0 )
                {
                    return RbConstants::Double::neginf;
                }
            }
        }
    }
    
    z.setCholesky(true);

    // check that the matrix is positive definite
    if ( z.getCholeskyDecomposition().checkPositiveSemidefinite() == false )
    {
//        std::cout << "Rejecting non-positive-definite matrix!" << std::endl;
        return RbConstants::Double::neginf;
    }
    
    // compute the normalizing constant
    // equation 16 from Lewandowski, Kurowicka and Joe 2009
    double s = 0.0;
    double p = 0.0;
    for (int k = 1; k <= dim - 1; ++k)
    {
        s += (2 * eta - 2 + dim - k) * (dim - k);
        p += (dim - k) * RbMath::lnBeta(eta + 0.5 * (dim - k - 1), eta + 0.5 * (dim - k - 1));
    }
    
    double c = s * log(2) + p; // the log normalizing constant
    
    return c + (eta - 1) * z.getLogDet();

}



/*!
 * This function generates a LKJ-distributed random variable.
 *
 * \brief LKJ random variable.
 * \param eta is the parameter of the LKJ distribution
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the LKJ random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::LKJ::rv(double eta, size_t dim, RandomNumberGenerator& rng)
{

    MatrixReal P(dim); // this matrix holds the partial correlations
    MatrixReal S(dim); // this matrix holds the product-moment correlation matrix
    
    // S must start as an identity matrix
    for (size_t i = 0; i < dim; ++i)
    {
        S[i][i] = 1.0;
    }
    
    // this algorithm is only safe for eta > 1.
    // if eta < 1, we'll just return the identity matrix.
    // during MCMC, this only means that the initial value will not be from the prior
    if (eta < 1)
    {
        return S;
    }
    
    // initialize beta
    double beta = eta + (dim - 1) / 2;
    
    for (int k = 0; k < dim - 1; ++k)
    {
        
        // decrement beta
        beta -= 0.5;
        
        for (int i = k + 1; i < dim; ++i)
        {
            
            P[k][i] = RbStatistics::Beta::rv(beta, beta, rng) * 2.0 - 1.0; // sample the partial correlation
            
            double p = P[k][i]; // initial value for the product-moment correlation
            for (int l = k - 1; l >= 0; --l)
            {
                p = p * pow( ( 1.0 - pow(P[l][i], 2) ) * ( 1.0 - pow(P[l][k], 2) ), 0.5) + P[l][i] * P[l][k];
            }
            
            S[i][k] = p;
            S[k][i] = p;
            
        }
        
    }
    
    return S;
    
}

/*!
 * This function calculates the probability density
 * for a LKJ-distributed random variable.
 *
 * \brief LKJ probability density.
 * \param eta is the parameter of the LKJ distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::LKJ::pdfPartial(double eta, const MatrixReal &z) {
    
    return exp(lnPdfPartial(eta, z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a LKJ-distributed random variable.
 *
 * \brief Natural log of LKJ probability density.
 * \param eta is the parameter of the LKJ distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::LKJ::lnPdfPartial(double eta, const MatrixReal &z)
{
    
    size_t dim = z.getNumberOfRows();
    
    // compute the correlation matrix
    MatrixReal P(dim);
    
    // P must start as an identity matrix
    for (size_t i = 0; i < dim; ++i)
    {
        P[i][i] = 1.0;
    }

    for (int k = 0; k < dim - 1; ++k)
    {
        for (int i = k + 1; i < dim; ++i)
        {
            double p = z[k][i]; // initial value for the product-moment correlation
            for (int l = k - 1; l >= 0; --l)
            {
                p = p * pow( ( 1.0 - pow(z[l][i], 2) ) * ( 1.0 - pow(z[l][k], 2) ), 0.5) + z[l][i] * z[l][k];
            }
            P[i][k] = p;
            P[k][i] = p;
        }
    }
    
    // compute the Jacobian
    // equation 11 from Lewandowski, Kurowicka and Joe 2009
    double ln_jacobian = 0.0;
    for (size_t k = 0; k < dim - 2; ++k)
    {
        for (size_t i = k + 1; i < dim; ++i)
        {
            ln_jacobian += (dim - k - 2) * log(1.0 - pow(z[k][i], 2.0));
        }
    }
    ln_jacobian *= -0.5;

    // compute the probability of the induced correlation matrix
    double ln_prob = lnPdf(eta, P);
    
    return ln_prob - ln_jacobian;
    
}


/*!
 * This function generates a LKJ-distributed random variable (the partial correlations).
 *
 * \brief LKJ random variable.
 * \param eta is the parameter of the LKJ distribution
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the LKJ random variable (the partial correlations).
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::LKJ::rvPartial(double eta, size_t dim, RandomNumberGenerator& rng)
{
    
    MatrixReal P(dim); // this matrix holds the partial correlations
    for (int i = 0; i < dim; ++i)
    {
        P[i][i] = 1.0;
    }
    
    // I think this algorithm is only safe for eta > 1.
    // if eta < 1, we'll just return the identity matrix.
    // during MCMC, this only means that the initial value will not be from the prior
    if (eta < 1)
    {
        return P;
    }
    
    // initialize beta
    double beta = eta + (dim - 1) / 2;
    
    for (int k = 0; k < dim - 1; ++k)
    {
        // decrement beta
        beta -= 0.5;
        for (int i = k + 1; i < dim; ++i)
        {
            P[k][i] = RbStatistics::Beta::rv(beta, beta, rng) * 2.0 - 1.0; // sample the partial correlation
            P[i][k] = P[k][i];
        }
    }
    
    return P;
    
}

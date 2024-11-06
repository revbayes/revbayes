/*
 * File:   DistributionExponentialError.cpp
 * Author: David Cerny
 *
 * Created on October 8, 2019
 */


#include "DistributionExponentialError.h"

#include <cmath>
#include <cstdint>
#include <numeric>
#include <algorithm>

#include "DistributionDirichlet.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"

using namespace RevBayesCore;



/*!
 * This function calculates the probability density for a distance matrix
 * whose distance from the average distance matrix is exponentially
 * distributed.
 *
 * \brief ExponentialError probability density.
 * \param avgDistMat is a reference to the average distance matrix.
 * \param lambda is the rate parameter of the exponential.
 * \param z is a reference to a distance matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ExponentialError::pdf(const AverageDistanceMatrix &avgDistMat, double lambda, const AverageDistanceMatrix &z)
{
    return exp(lnPdf(avgDistMat, lambda, z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a distance matrix whose distance from the average distance
 * matrix is exponentially distributed.
 *
 * \brief Natural log of ExponentialError probability density.
 * \param avgDistMat is a reference to the average distance matrix.
 * \param lambda is the rate parameter of the exponential.
 * \param z is a reference to a distance matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ExponentialError::lnPdf(const AverageDistanceMatrix &avgDistMat, double lambda, const AverageDistanceMatrix &z)
{
    
    if ( avgDistMat.getSize() != z.getSize() )
    {
        return RbConstants::Double::neginf;
    }
    
    /* AverageDistanceMatrix objects already have their rows and columns alphabetically sorted,
     * so there is no need to further align them against each other by index matching (like we
     * do in the overloaded function below)
     */

    double dist = 0;

    for (size_t i = 0; i < avgDistMat.getSize(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (z.getMask()[i][j])
            {
                double difference = z.getDistanceMatrix()[i][j] - avgDistMat.getDistanceMatrix()[i][j];
                dist += difference * difference;
            }
        }
    }
    
    return (std::log(lambda) - lambda * dist);

}



/*!
 * This function generates an "average" distance matrix of random
 * variables whose distance from the user-specified average distance
 * matrix is exponentially distributed.
 *
 * \brief ExponentialError random variable.
 * \param avgDistMat is a reference to the average distance matrix.
 * \param lambda is the rate parameter of the exponential.
 * \param rng is a pointer to a random number object.
 * \return Returns an average distance matrix of exponential error-distributed random variables.
 * \throws Does not throw an error.
 */
AverageDistanceMatrix RbStatistics::ExponentialError::rv(const AverageDistanceMatrix &avgDistMat, double lambda, RandomNumberGenerator& rng)
{
 
    size_t dim = avgDistMat.getSize();
    
    // first, draw from the exponential of rate lambda:
    double u = rng.uniform01();
    double draw = -(1.0/lambda) * std::log(u);
    
    // next, get a vector of draws from a Dirichlet distribution that has half as many
    // elements as there are off-diagonal elements in the average distance matrix:
    size_t vect_len = (dim * (dim - 1))/2;
    std::vector<double> DirichParams(vect_len, 1);
    std::vector<double> DirichRandomVars = RbStatistics::Dirichlet::rv(DirichParams, rng);
    
    // get an equally std::int64_t vector of signs (-1 or 1):
    std::vector<double> Signs(vect_len);
    for (size_t i=0; i<vect_len; i++)
    {
        double signProb = rng.uniform01();
        if (signProb <= 0.5)
        {
            Signs[i] = -1;
        } else {
            Signs[i] = 1;
        }
    }
    
    // initialize the z and y matrices from avgDistMat using the copy constructor:
    DistanceMatrix z = DistanceMatrix( avgDistMat.getDistanceMatrix() );
    MatrixBoolean y = MatrixBoolean( avgDistMat.getMask() );
    
    // fill it in:
    for (size_t i=0; i<dim; i++)
    {
        for (size_t j=0; j<dim; j++)
        {
            if (i > j) // lower triangular
            {
                size_t index = (j * (2 * dim - j - 1))/2 + (i - j - 1);
                z[i][j] += sqrt(draw * DirichRandomVars[index]) * Signs[index];
            }
            else if (j > i) // upper triangular
            {
                size_t index = (i * (2 * dim - i - 1))/2 + (j - i - 1);
                z[i][j] += sqrt(draw * DirichRandomVars[index]) * Signs[index];
            }
        }
    }
 
    // construct a new AverageDistanceMatrix from the newly calculated distance matrix and the old mask
    AverageDistanceMatrix x = AverageDistanceMatrix( z, y );
    return x;
}



/*!
* This function calculates the probability density for a distance matrix
* whose distance from the another distance matrix is exponentially
* distributed.
*
* \brief ExponentialError probability density.
* \param distMat is a reference to the user-specified distance matrix.
* \param lambda is the rate parameter of the exponential.
* \param z is a reference to a distance matrix containing the random variables.
* \return Returns the probability density.
* \throws Does not throw an error.
*/
double RbStatistics::ExponentialError::pdf(const DistanceMatrix &distMat, double lambda, const AverageDistanceMatrix &z)
{
    return exp(lnPdf(distMat, lambda, z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a distance matrix whose distance from another distance
 * matrix is exponentially distributed.
 *
 * \brief Natural log of ExponentialError probability density.
 * \param distMat is a reference to the user-specified distance matrix.
 * \param lambda is the rate parameter of the exponential.
 * \param z is a reference to a distance matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ExponentialError::lnPdf(const DistanceMatrix &distMat, double lambda, const AverageDistanceMatrix &z)
{
    
    if ( distMat.getSize() != z.getSize() )
    {
        return RbConstants::Double::neginf;
    }

    std::vector<std::string> dmNames( distMat.getSize() );
    for(size_t i = 0; i < dmNames.size(); i++)
    {
        dmNames[i] = distMat.getTaxa()[i].getName();
    }

    std::vector<uint32_t> ix = StringUtilities::stringSortIndices(dmNames);
    
    double dist = 0;
    
    for (size_t i = 0; i < distMat.getSize(); i++)
    {
        uint32_t k = ix[i];
        for (size_t j = 0; j < i; j++)
        {
            uint32_t l = ix[j];
            if (z.getMask()[i][j])
            {
                double difference = z.getDistanceMatrix()[i][j] - distMat[k][l];
                dist += difference * difference;
            }
        }
    }
    
    return (std::log(lambda) - lambda * dist);

}



/*!
 * This function generates a distance matrix of random variables whose
 * distance from the user-specified distance matrix is exponentially
 * distributed.
 *
 * \brief ExponentialError random variable.
 * \param distMat is a reference to the user-specified distance matrix.
 * \param lambda is the rate parameter of the exponential.
 * \param rng is a pointer to a random number object.
 * \return Returns a distance matrix of exponential error-distributed random variables.
 * \throws Does not throw an error.
 */
AverageDistanceMatrix RbStatistics::ExponentialError::rv(const DistanceMatrix &distMat, double lambda, RandomNumberGenerator& rng)
{
 
    size_t dim = distMat.getSize();
    
    // first, draw from the exponential of rate lambda:
    double u = rng.uniform01();
    double draw = -(1.0/lambda) * std::log(u);
    
    // next, get a vector of draws from a Dirichlet distribution that has half as many
    // elements as there are off-diagonal elements in the average distance matrix:
    size_t vect_len = (dim * (dim - 1))/2;
    std::vector<double> DirichParams(vect_len, 1);
    std::vector<double> DirichRandomVars = RbStatistics::Dirichlet::rv(DirichParams, rng);
    
    // get an equally std::int64_t vector of signs (-1 or 1):
    std::vector<double> Signs(vect_len);
    for (size_t i=0; i<vect_len; i++)
    {
        double signProb = rng.uniform01();
        if (signProb <= 0.5)
        {
            Signs[i] = -1;
        } else {
            Signs[i] = 1;
        }
    }
    
    // initialize the z matrix from distMat using the copy constructor:
    DistanceMatrix z = DistanceMatrix( distMat );
    
    // fill it in:
    for (size_t i=0; i<dim; i++)
    {
        for (size_t j=0; j<dim; j++)
        {
            if (i > j) // lower triangular
            {
                size_t index = (j * (2 * dim - j - 1))/2 + (i - j - 1);
                z[i][j] += sqrt(draw * DirichRandomVars[index]) * Signs[index];
            }
            else if (j > i) // upper triangular
            {
                size_t index = (i * (2 * dim - i - 1))/2 + (j - i - 1);
                z[i][j] += sqrt(draw * DirichRandomVars[index]) * Signs[index];
            }
        }
    }
    
    // get a mask of the right dimension, with all bits flipped to 1 (= all distances treated as meaningful)
    MatrixBoolean y = MatrixBoolean( dim ).negate();
    AverageDistanceMatrix x = AverageDistanceMatrix( z, y );
    return x;
    
}

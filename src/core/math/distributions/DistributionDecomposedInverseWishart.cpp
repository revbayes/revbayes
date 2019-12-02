/* 
 * File:   DistributionInverseInverseWishart.cpp
 * Author: nl
 * 
 * Created on 15 juillet 2014, 10:13
 */


#include "DistributionDecomposedInverseWishart.h"
#include "DistributionMultivariateNormal.h"
#include "DistributionUniform.h"
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
double RbStatistics::DecomposedInverseWishart::pdf(double nu, const MatrixReal &r) {
	
    return exp(lnPdf(nu, r));
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

double RbStatistics::DecomposedInverseWishart::lnPdf(double nu, const MatrixReal &r) {

    size_t k = r.getDim();
    
    if ( r.isPositiveDefinite() == false )
        {
        return RbConstants::Double::neginf;
        }
    
    double lnP = (0.5 * (k - 1.0) * (nu - 1.0) - 1.0);
    lnP += r.getLogDet();
    
    MatrixReal submatrix(k-1);
    for (size_t i=0; i<k; i++)
        {
        size_t ai = 0;
        for (size_t a=0; a<k; a++)
            {
            if (a != i)
                {
                size_t bi = 0;
                for (size_t b=0; b<k; b++)
                    {
                    if ( b != i )
                        {
                        submatrix[ai][bi] = r[a][b];
                        bi++;
                        }
                    }
                ai++;
                }
            }

        lnP += submatrix.getLogDet();
//        std::cout << "logdet=" << submatrix.getLogDet() << std::endl;
        }

    return lnP;
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
MatrixReal RbStatistics::DecomposedInverseWishart::rv(size_t k, double nu, RandomNumberGenerator& rng) {
        
    MatrixReal z(k);
    
    for (int i=0; i<k; i++)
        {
        for (int j=0; j<k; j++)
            {
            if (i == j)
                z[i][j] = 1.0;
            else
                z[i][j] = 0.0;
            }
        }
    
    /*for (int i=0; i<k; i++)
        {
        for (int j=0; j<k; j++)
            {
            double r_ij = RbStatistics::Uniform::rv(-1.0, 1.0, rng);
            z[i][j] = r_ij;
            }
        }*/

    return z;
}



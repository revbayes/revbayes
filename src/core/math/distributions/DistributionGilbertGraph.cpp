//
//  DistributionGilbertGraph.cpp
//  revbayes
//
//  Created by Michael R. May on 2017-08-04.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#include <cstddef>
#include <cmath>

#include "DistributionGilbertGraph.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a GilbertGraph-distributed random variable.
 *
 * \brief GilbertGraph probability density.
 * \param eta is the parameter of the GilbertGraph distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::GilbertGraph::pdf(double p, const MatrixReal &z) {
    
    return exp(lnPdf(p, z));
}



/*!
 * This function calculates the natural log of the probability density
 * for a GilbertGraph-distributed random variable.
 *
 * \brief Natural log of GilbertGraph probability density.
 * \param eta is the parameter of the GilbertGraph distribution
 * \param z is a reference to a matrix containing the random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::GilbertGraph::lnPdf(double p, const MatrixReal &z)
{
    
    size_t dim = z.getNumberOfRows();
    size_t num_possible_edges = dim * (dim-1) / 2;
    size_t num_assigned_edges = 0;
    for (size_t i = 0; i < z.size(); i++) {
        for (size_t j = i+1; j < z.size(); j++) {
            num_assigned_edges += z[i][j];
        }
    }
    size_t num_unassigned_edges = num_possible_edges - num_assigned_edges;
    double ret = num_assigned_edges * log(p) + num_unassigned_edges * log(1.0-p);
    return ret;
}



/*!
 * This function generates a GilbertGraph-distributed random variable.
 *
 * \brief GilbertGraph random variable.
 * \param eta is the parameter of the GilbertGraph distribution
 * \param rng is a pointer to a random number object.
 * \return Returns a vector containing the GilbertGraph random variable.
 * \throws Does not throw an error.
 */
MatrixReal RbStatistics::GilbertGraph::rv(double p, size_t dim, RandomNumberGenerator& rng)
{
    
    MatrixReal P(dim, dim, 0);
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = i+1; j < dim; j++)
        {
            if (GLOBAL_RNG->uniform01() < p) {
                P[i][j] = 1;
                P[j][i] = 1;
            }
        }
    }
    
    return P;
    
}

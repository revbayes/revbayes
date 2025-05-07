#include <cstddef>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

#include "HalfCauchyDistribution.h"
#include "GMRFUnevenGridHyperpriorGibbsMove.h"
#include "NormalDistribution.h"
#include "RandomNumberFactory.h"
#include "RbStatisticsHelper.h"
#include "TypedDagNode.h"
#include "AbstractGibbsMove.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
GMRFUnevenGridHyperpriorGibbsMove::GMRFUnevenGridHyperpriorGibbsMove( StochasticNode<double> *g, std::vector< StochasticNode<double> *> n, TypedDagNode< RbVector<double> > *u, double z, double w) : AbstractGibbsMove(w),
    global_scale( g ),
    normals( n ),
    grid( u ),
    zeta( z )
{
    
    // tell the base class to add the node
    assert( not global_scale->isClamped() );
    
    HalfCauchyDistribution* dist = dynamic_cast<HalfCauchyDistribution *>( &global_scale->getDistribution() );
    
    if (dist == NULL)
    {
        throw(RbException("GMRFUnevenGridHyperpriorGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }
    
    if (dist->getScale()->getValue() != 1.0)
    {
        throw(RbException("GMRFUnevenGridHyperpriorGibbsMove requiresglobal scale has Half-Cauchy(0,1) Distribution"));
    }
    
    if (dist->getLocation()->getValue() != 0.0)
    {
        throw(RbException("GMRFUnevenGridHyperpriorGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }
    
    addNode( global_scale );
    
    // tell the base class to add the node
    for (size_t i=0; i<normals.size(); ++i)
    {
        
        NormalDistribution* dist = dynamic_cast<NormalDistribution *>( &normals[i]->getDistribution() );
        
        if (dist == NULL)
        {
            throw(RbException("GMRFUnevenGridHyperpriorGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }
        
        if (dist->getMean()->getValue() != 0.0)
        {
            throw(RbException("GMRFUnevenGridHyperpriorGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }
        
        const TypedDagNode<double>* sd = dist->getStDev();

        // Make sure that the normal distributions have the correct stdevs for this sampler to be appropriate
        // Somewhere there is a very slight deviance being introduced, hence the addition of rounding in this comparison
        if (round(100000 * sd->getValue()) != round(100000 * zeta * global_scale->getValue() * sqrt(grid->getValue()[i]))) {
            throw(RbException("GMRFUnevenGridHyperpriorGibbsMove move only works when children are Normal(0,global_scale*zeta*sqrt(grid[i])) Distributions"));
        }

        addNode( normals[i] );
    }

    
}


/**
 * Basic destructor doing nothing.
 */
GMRFUnevenGridHyperpriorGibbsMove::~GMRFUnevenGridHyperpriorGibbsMove( void )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the GMRFUnevenGridHyperpriorGibbsMove. 
 */
GMRFUnevenGridHyperpriorGibbsMove* GMRFUnevenGridHyperpriorGibbsMove::clone( void ) const 
{
    return new GMRFUnevenGridHyperpriorGibbsMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& GMRFUnevenGridHyperpriorGibbsMove::getMoveName( void ) const 
{
    static std::string name = "GMRFUnevenGridHyperpriorGibbs";

    return name;
}


/** Perform the move */
void GMRFUnevenGridHyperpriorGibbsMove::performGibbsMove( void )
{
    
    double n = normals.size() + 1.0; // The field has one more cell than it has increments
    
    // get global scale and sample auxiliary variable
    double eta = global_scale->getValue();
    
    double xi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + 1.0/(eta*eta)); // xi_inverse ~ Gamma(1, 1 + 1/eta^2)
    
    // Get the spacing of the grid points
    RbVector<double> u = grid->getValue();
    
    // sample global scale
    double eta_squared_rate = 0;
    for (size_t i = 0; i < normals.size(); ++i) {
        eta_squared_rate += std::pow(normals[i]->getValue(),2.0)/u[i];
    }
    eta_squared_rate *= 1/(2*zeta*zeta);
    eta_squared_rate += xi_inverse;

    double eta_squared_inverse = RbStatistics::Helper::rndGamma(0.5*n, *GLOBAL_RNG) / eta_squared_rate; // eta_squared ~ InverseGamma(n/2, eta_squared_rate)

    eta = std::sqrt(1/eta_squared_inverse);

    global_scale->getValue() = eta;
    global_scale->touch();
    global_scale->keep();

}


void GMRFUnevenGridHyperpriorGibbsMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < normals.size(); ++i)
    {
        if ( normals[i] == oldN )
        {
            normals[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
    // also swap the sd parameter
    if (oldN == global_scale)
    {
        global_scale = static_cast<StochasticNode<double>* >(newN) ;
    }
    
}

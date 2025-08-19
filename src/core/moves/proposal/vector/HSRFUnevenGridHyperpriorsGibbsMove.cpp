#include <cstddef>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

#include "HalfCauchyDistribution.h"
#include "HSRFUnevenGridHyperpriorsGibbsMove.h"
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
HSRFUnevenGridHyperpriorsGibbsMove::HSRFUnevenGridHyperpriorsGibbsMove( StochasticNode<double> *g, std::vector< StochasticNode<double> *> l, std::vector< StochasticNode<double> *> n, TypedDagNode< RbVector<double> > *u, double z, double w) : AbstractGibbsMove(w),
    global_scale( g ),
    local_scales( l ),
    grid( u ),
    normals( n ),
    zeta( z )
{
    
    
    // tell the base class to add the node
    assert( not global_scale->isClamped() );
    
    HalfCauchyDistribution* dist = dynamic_cast<HalfCauchyDistribution *>( &global_scale->getDistribution() );
    
    if (dist == NULL)
    {
        throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }
    
    if (dist->getScale()->getValue() != 1.0)
    {
        throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requiresglobal scale has Half-Cauchy(0,1) Distribution"));
    }
    
    if (dist->getLocation()->getValue() != 0.0)
    {
        throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }
    
    addNode( global_scale );
    
    // tell the base class to add the node
    for (size_t i=0; i<local_scales.size(); ++i)
    {
        assert( not local_scales[i]->isClamped() );
        
        HalfCauchyDistribution* dist = dynamic_cast<HalfCauchyDistribution *>( &local_scales[i]->getDistribution() );
        
        if (dist == NULL)
        {
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }

        if (dist->getScale()->getValue() != 1.0)
        {
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }

        if (dist->getLocation()->getValue() != 0.0)
        {
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }
        
        addNode( local_scales[i] );
    }
    
    // tell the base class to add the node
    for (size_t i=0; i<normals.size(); ++i)
    {
        
        NormalDistribution* dist = dynamic_cast<NormalDistribution *>( &normals[i]->getDistribution() );
        
        if (dist == NULL)
        {
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }
        
        if (dist->getMean()->getValue() != 0.0)
        {
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }
        
        const TypedDagNode<double>* sd = dist->getStDev();
        
        // Make sure that the normal distributions have the correct stdevs for this sampler to be appropriate
        // Somewhere there is a very slight deviance being introduced, hence the addition of rounding in this comparison
        if (round(100000 * sd->getValue()) != round(100000 * zeta * global_scale->getValue() * local_scales[i]->getValue() * sqrt(grid->getValue()[i])) ) {
            std::cout << "zeta, global_scale, and local_scale[i] are " << zeta << "," << global_scale->getValue() << "," << local_scales[i]->getValue()<< std::endl;
            std::cout << "sd is " << sd->getValue() << std::endl;
            throw(RbException("HSRFUnevenGridHyperpriorsGibbsMove move only works when children are Normal(0,local_scale*global_scale*zeta*sqrt(grid[i])) Distributions"));
        }
        
        addNode( normals[i] );
    }

}


/**
 * Basic destructor doing nothing.
 */
HSRFUnevenGridHyperpriorsGibbsMove::~HSRFUnevenGridHyperpriorsGibbsMove( void )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the HSRFUnevenGridHyperpriorsGibbsMove. 
 */
HSRFUnevenGridHyperpriorsGibbsMove* HSRFUnevenGridHyperpriorsGibbsMove::clone( void ) const 
{
    return new HSRFUnevenGridHyperpriorsGibbsMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& HSRFUnevenGridHyperpriorsGibbsMove::getMoveName( void ) const 
{
    static std::string name = "HSRFUnevenGridHyperpriorsGibbs";

    return name;
}

/** Perform the move */
void HSRFUnevenGridHyperpriorsGibbsMove::performGibbsMove( void )
{
//    std::cout << "Hello from HSRFUnevenGridHyperpriorsGibbsMove::performGibbsMove( void )" << std::endl;
    
    // get global scale
    double eta_squared = std::pow(global_scale->getValue(),2.0);

    // Sample auxiliary variables and field local scale parameters, keep tabs on needed values for global scale
    size_t field_size = local_scales.size();
    double n = local_scales.size() + 1; // The field has one more cell than there are local scales
    double eta_squared_rate = 0;
    
    // sample local scales
    double two_zeta_squared = 2.0 * std::pow(zeta,2.0);
    double two_zeta_squared_eta_squared = two_zeta_squared * eta_squared;
    
    // Get the spacing of the grid points
    RbVector<double> u = grid->getValue();
    
    for (size_t i = 0; i < field_size; ++i) {
        double lambda_squared_inverse = 1.0 / std::pow(local_scales[i]->getValue(),2.0);
        
        double delta_theta_squared = std::pow(normals[i]->getValue(),2.0)/u[i];
        
        double psi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + lambda_squared_inverse); // psi_inverse ~ Gamma(1.0, 1.0 + 1.0/lambda[i]^2)
        
        double lambda_squared = 1/(RbStatistics::Helper::rndGamma(1, *GLOBAL_RNG) / (psi_inverse + delta_theta_squared/two_zeta_squared_eta_squared) ); // lambda_squared[i] ~ Gamma(0.5, psi_inverse[i] + delta_theta^2/(2*eta^2*zeta^2)
        
        local_scales[i]->getValue() = std::sqrt(lambda_squared);
        local_scales[i]->touch();
        local_scales[i]->keep();
        
        eta_squared_rate += delta_theta_squared/lambda_squared;
    }
    
    // get global scale and sample global scale auxiliary variable
    double xi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + 1.0/(eta_squared)); // xi_inverse ~ Gamma(1, 1 + 1/eta^2)
    
    // sample global scale
    eta_squared_rate *= 1/two_zeta_squared;
    eta_squared_rate += xi_inverse;
    
    double eta_squared_inverse = RbStatistics::Helper::rndGamma(0.5*n, *GLOBAL_RNG) / eta_squared_rate; // eta_squared ~ InverseGamma(n/2, eta_squared_rate)
    
    global_scale->getValue() = std::sqrt(1/eta_squared_inverse);
    global_scale->touch();
    global_scale->keep();

    
}


void HSRFUnevenGridHyperpriorsGibbsMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < local_scales.size(); ++i)
    {
        if ( local_scales[i] == oldN )
        {
            local_scales[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
    for (size_t i = 0; i < normals.size(); ++i)
    {
        if ( normals[i] == oldN )
        {
            normals[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }
    
    if (oldN == grid)
    {
        grid = static_cast<TypedDagNode< RbVector<double> >* >( newN );
    }
    
    if (oldN == global_scale)
    {
        global_scale = static_cast<StochasticNode<double>* >(newN) ;
    }

}

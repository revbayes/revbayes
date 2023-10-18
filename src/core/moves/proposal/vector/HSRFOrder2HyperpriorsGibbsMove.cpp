#include <cstddef>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

#include "HalfCauchyDistribution.h"
#include "HSRFOrder2HyperpriorsGibbsMove.h"
#include "NormalDistribution.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbStatisticsHelper.h"
#include "TypedDagNode.h"
#include "AbstractGibbsMove.h"
#include "RbException.h"
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
HSRFOrder2HyperpriorsGibbsMove::HSRFOrder2HyperpriorsGibbsMove( StochasticNode<double> *g, std::vector< StochasticNode<double> *> l, std::vector< StochasticNode<double> *> n, double p, double z, double w) : AbstractGibbsMove(w),
    global_scale( g ),
    local_scales( l ),
    normals( n ),
    prop_global_only( p ),
    zeta( z )
{


    // tell the base class to add the node
    assert( not global_scale->isClamped() );

    HalfCauchyDistribution* dist = dynamic_cast<HalfCauchyDistribution *>( &global_scale->getDistribution() );

    if (dist == NULL)
    {
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }

    if (dist->getScale()->getValue() != 1.0)
    {
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove requiresglobal scale has Half-Cauchy(0,1) Distribution"));
    }

    if (dist->getLocation()->getValue() != 0.0)
    {
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove requires global scale has Half-Cauchy(0,1) Distribution"));
    }

    addNode( global_scale );

    // tell the base class to add the node
    for (size_t i=0; i<local_scales.size(); ++i)
    {
        assert( not local_scales[i]->isClamped() );

        HalfCauchyDistribution* dist = dynamic_cast<HalfCauchyDistribution *>( &local_scales[i]->getDistribution() );

        if (dist == NULL)
        {
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }

        if (dist->getScale()->getValue() != 1.0)
        {
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }

        if (dist->getLocation()->getValue() != 0.0)
        {
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove requires all local scales have Half-Cauchy(0,1) Distribution"));
        }

        addNode( local_scales[i] );
    }

    // tell the base class to add the node
    for (size_t i=1; i<normals.size(); ++i)
    {

        NormalDistribution* dist = dynamic_cast<NormalDistribution *>( &normals[i]->getDistribution() );

        if (dist == NULL)
        {
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }

        if (dist->getMean()->getValue() != 0.0)
        {
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
        }

        const TypedDagNode<double>* sd = dist->getStDev();

        // Make sure that the normal distributions have the correct stdevs for this sampler to be appropriate
        // Somewhere there is a very slight deviance being introduced, hence the addition of rounding in this comparison
        if (round(100000 * sd->getValue()) != round(100000 * zeta * global_scale->getValue() * local_scales[i]->getValue()) ) {
            std::cout << "zeta, global_scale, and local_scale[i] are " << zeta << "," << global_scale->getValue() << "," << local_scales[i]->getValue()<< std::endl;
            std::cout << "sd is " << sd->getValue() << std::endl;
            throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children 2:n are Normal(0,local_scale*global_scale*zeta) Distributions"));
        }

        addNode( normals[i] );
    }

    NormalDistribution* dist0 = dynamic_cast<NormalDistribution *>( &normals[0]->getDistribution() );

    if (dist0 == NULL)
    {
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
    }

    if (dist0->getMean()->getValue() != 0.0)
    {
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children are Normal(0,sigma) Distributions"));
    }

    const TypedDagNode<double>* sd = dist0->getStDev();

    // Make sure that the normal distributions have the correct stdevs for this sampler to be appropriate
    // Somewhere there is a very slight deviance being introduced, hence the addition of rounding in this comparison
    if (round(100000 * sd->getValue()) != round(100000 * RbConstants::SQRT1_2 * zeta * global_scale->getValue() * local_scales[0]->getValue()) ) {
        std::cout << "zeta, global_scale, and local_scale[0] are " << zeta << "," << global_scale->getValue() << "," << local_scales[0]->getValue()<< std::endl;
        std::cout << "sd is " << sd->getValue() << std::endl;
        throw(RbException("HSRFOrder2HyperpriorsGibbsMove move only works when children[0] is a Normal(0,sqrt(0.5)*local_scale*global_scale*zeta) Distribution"));
    }

    addNode( normals[0] );

}


/**
 * Basic destructor doing nothing.
 */
HSRFOrder2HyperpriorsGibbsMove::~HSRFOrder2HyperpriorsGibbsMove( void )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the HSRFOrder2HyperpriorsGibbsMove.
 */
HSRFOrder2HyperpriorsGibbsMove* HSRFOrder2HyperpriorsGibbsMove::clone( void ) const
{
    return new HSRFOrder2HyperpriorsGibbsMove( *this );
}


/**
 * Get moves' name of object
 *
 * \return The moves' name.
 */
const std::string& HSRFOrder2HyperpriorsGibbsMove::getMoveName( void ) const
{
    static std::string name = "HSRFOrder2HyperpriorsGibbs";

    return name;
}

/**
 * The parameter this Gibbs move operates on is always in the prior.
 * If it is, then the conditional distribution is safe if the likelihood is heated.
 */
bool HSRFOrder2HyperpriorsGibbsMove::heatsAreAllowable(double prHeat, double lHeat, double pHeat)
{
  if ( prHeat == 1.0 && lHeat == 1.0 && pHeat == 1.0 )
  {
    return true;
  }

  if ( prHeat != 1.0 || pHeat != 1.0 ) {
    return false;
  }

  if ( lHeat != 1.0 )
  {
    for (size_t i=0; i<normals.size(); ++i)
    {
      if ( normals[i]->isClamped() )
      {
        return false;
      }
    }
    return true;
  }

  return false;
}

/** Perform the move */
void HSRFOrder2HyperpriorsGibbsMove::performGibbsMove( void )
{
//    std::cout << "Hello from HSRFOrder2HyperpriorsGibbsMove::performGibbsMove( void )" << std::endl;

    // get global scale
    double eta_squared = std::pow(global_scale->getValue(),2.0);

    // Sample auxiliary variables and field local scale parameters, keep tabs on needed values for global scale
    size_t field_size = local_scales.size();
    double n = local_scales.size() + 1; // The field has one more cell than there are local scales
    double eta_squared_rate = 0;
    
    // sample local scales
    double two_zeta_squared = 2.0 * std::pow(zeta,2.0);
    double two_zeta_squared_eta_squared = two_zeta_squared * eta_squared;


    // prop_global_only >= DBL_EPSILON keeps RNG output unchanged with analyses run before prop_global_only was added
    if ( prop_global_only >= DBL_EPSILON && (GLOBAL_RNG)->uniform01() < prop_global_only )
    {
        // First local scale
        double lambda_squared = std::pow(local_scales[0]->getValue(),2.0);
        double delta_theta_squared = std::pow(normals[0]->getValue(),2.0);
        eta_squared_rate += 2.0 * delta_theta_squared/lambda_squared;

        // All others
        for (size_t i = 1; i < field_size; ++i) {
            lambda_squared = std::pow(local_scales[i]->getValue(),2.0);

            delta_theta_squared = std::pow(normals[i]->getValue(),2.0);

            eta_squared_rate += delta_theta_squared/lambda_squared;
        }
    }
    else
    {
        // First local scale (different when drawing lambda, same for psi)
        double lambda_squared_inverse = 1.0 / std::pow(local_scales[0]->getValue(),2.0);

        double delta_theta_squared = std::pow(normals[0]->getValue(),2.0);

        double psi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + lambda_squared_inverse); // psi_inverse ~ Gamma(1.0, 1.0 + 1.0/lambda[i]^2)

        double lambda_squared = 1/(RbStatistics::Helper::rndGamma(1, *GLOBAL_RNG) / (psi_inverse + (2.0*delta_theta_squared)/two_zeta_squared_eta_squared) ); // lambda_squared[i] ~ Gamma(0.5, psi_inverse[i] + delta_theta^2/(2*eta^2*zeta^2)

        local_scales[0]->getValue() = std::sqrt(lambda_squared);
        local_scales[0]->touch();

        eta_squared_rate += 2.0 * delta_theta_squared/lambda_squared;

        // All others
        for (size_t i = 1; i < field_size; ++i) {
            lambda_squared_inverse = 1.0 / std::pow(local_scales[i]->getValue(),2.0);

            delta_theta_squared = std::pow(normals[i]->getValue(),2.0);

            psi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + lambda_squared_inverse); // psi_inverse ~ Gamma(1.0, 1.0 + 1.0/lambda[i]^2)

            lambda_squared = 1/(RbStatistics::Helper::rndGamma(1, *GLOBAL_RNG) / (psi_inverse + delta_theta_squared/two_zeta_squared_eta_squared) ); // lambda_squared[i] ~ Gamma(0.5, psi_inverse[i] + delta_theta^2/(2*eta^2*zeta^2)

            local_scales[i]->getValue() = std::sqrt(lambda_squared);
            local_scales[i]->touch();

            eta_squared_rate += delta_theta_squared/lambda_squared;
        }
    }

    // get global scale and sample global scale auxiliary variable
    double xi_inverse = RbStatistics::Helper::rndGamma(1.0, *GLOBAL_RNG) / (1.0 + 1.0/(eta_squared)); // xi_inverse ~ Gamma(1, 1 + 1/eta^2)

    // sample global scale
    eta_squared_rate *= 1/two_zeta_squared;
    eta_squared_rate += xi_inverse;

    double eta_squared_inverse = RbStatistics::Helper::rndGamma(0.5*n, *GLOBAL_RNG) / eta_squared_rate; // eta_squared ~ InverseGamma(n/2, eta_squared_rate)

    global_scale->getValue() = std::sqrt(1/eta_squared_inverse);
    global_scale->touch();

    // keep the variables
    for (size_t i = 0; i < field_size; ++i)
    {
      local_scales[i]->keep();
    }

}


void HSRFOrder2HyperpriorsGibbsMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
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

    if (oldN == global_scale)
    {
        global_scale = static_cast<StochasticNode<double>* >(newN) ;
    }

}

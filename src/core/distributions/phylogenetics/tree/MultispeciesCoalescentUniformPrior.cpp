#include "Clade.h"
#include "ConstantNode.h"
#include "MultispeciesCoalescentUniformPrior.h"
#include "DistributionExponential.h"
#include "DistributionUniform.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathFunctions.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

MultispeciesCoalescentUniformPrior::MultispeciesCoalescentUniformPrior(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : AbstractMultispeciesCoalescent(sp, t)
{
    
}


MultispeciesCoalescentUniformPrior::~MultispeciesCoalescentUniformPrior()
{
    
}





MultispeciesCoalescentUniformPrior* MultispeciesCoalescentUniformPrior::clone( void ) const
{
    
    return new MultispeciesCoalescentUniformPrior( *this );
}


double MultispeciesCoalescentUniformPrior::computeLnCoalescentProbability(size_t k, const std::vector<double> &times, double begin_age, double end_age, size_t index, bool add_final_interval)
{
    
    if ( k == 1 ) return 0.0;
    
    double theta_max = max_theta->getValue();
    
    double current_time = begin_age;
    double fn = 0.0;
    size_t n = times.size();
    for (size_t i=0; i<n; ++i)
    {
        
        // now we do the computation
        //a is the time between the previous and the current coalescences
        double a = times[i] - current_time;
        
        // get the number j of individuals we had before the current coalescence
        size_t j = k - i;
        double n_pairs = j * (j-1.0) / 2.0;

        fn += a * n_pairs;
    }
    
    // compute the probability of no coalescent event in the final part of the branch
    // only do this if the branch is not the root branch
    if ( add_final_interval == true )
    {
        double final_interval = end_age - current_time;
        size_t j = k - times.size();
        double n_pairs = j * (j-1.0) / 2.0;
        fn += final_interval * n_pairs;
        
    }
    
    double ln_prob_coal = RbConstants::LN2 - log( fn ) * (n-1) - log( theta_max );
    
//    shape, rate/x
//    double lowerIncompleteGamma = RbMath::incompleteGamma( rate/x, shape, RbMath::lnGamma(shape) );
//    double gamma = RbMath::gamma(shape);
    
//    Gamma(n-2,2*fn/theta_max)
    double lower_incomplete_gamma = RbMath::incompleteGamma( 2*fn/theta_max, n-1, RbMath::lnGamma(n-1) );
    
    double gamma = RbMath::lnGamma(n-1);
    gamma = 0.0;

    ln_prob_coal -= log( lower_incomplete_gamma ) - gamma;
    
    return ln_prob_coal;
}


double MultispeciesCoalescentUniformPrior::drawNe( size_t index )
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double u = RbStatistics::Uniform::rv( 0, max_theta->getValue(), *rng);
    
    return u;
}



void MultispeciesCoalescentUniformPrior::setMaxTheta(TypedDagNode<double>* m)
{
    
    removeParameter( max_theta );
    
    max_theta = m;
    
    addParameter( max_theta );
}


/** Swap a parameter of the distribution */
void MultispeciesCoalescentUniformPrior::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == max_theta )
    {
        max_theta = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    AbstractMultispeciesCoalescent::swapParameterInternal(oldP, newP);
    
}

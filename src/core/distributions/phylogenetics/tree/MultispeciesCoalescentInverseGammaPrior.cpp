#include <stddef.h>
#include <cmath>
#include <vector>

#include "MultispeciesCoalescentInverseGammaPrior.h"
#include "DistributionInverseGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "AbstractMultispeciesCoalescent.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

MultispeciesCoalescentInverseGammaPrior::MultispeciesCoalescentInverseGammaPrior(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : AbstractMultispeciesCoalescent(sp, t)
{

}


MultispeciesCoalescentInverseGammaPrior::~MultispeciesCoalescentInverseGammaPrior()
{

}





MultispeciesCoalescentInverseGammaPrior* MultispeciesCoalescentInverseGammaPrior::clone( void ) const
{

    return new MultispeciesCoalescentInverseGammaPrior( *this );
}


double MultispeciesCoalescentInverseGammaPrior::computeLnCoalescentProbability(size_t k, const std::vector<double> &times, double begin_age, double end_age, size_t index, bool add_final_interval)
{
    if ( k == 1 ) return 0.0;

    double alpha = shape->getValue();
    double beta = rate->getValue();

    double current_time = begin_age;

    size_t n = times.size();
    double a = n;

    double ln_prob_coal = 0.0;

    // If the number of gene copies is 1, then there can be no coalescence event
    // and the probability is equal to 1.0 for the only possible event (no coalescence)
    if ( n <= 1 )
    {
        ln_prob_coal = 0.0;
    }
    else
    {
        double b = 0.0;
        for (size_t i=0; i<n; ++i)
        {
            // now we do the computation
            // t is the time between the previous and the current coalescences
            double t = times[i] - current_time;
            current_time = times[i];

            // get the number j of individuals we had before the current coalescence
            size_t j = k - i;
            double n_pairs = j * (j-1.0) / 2.0;

            b += t * n_pairs;
        }

        // compute the probability of no coalescent event in the final part of the branch
        // only do this if the branch is not the root branch
        if ( add_final_interval == true )
        {
            double final_interval = end_age - current_time;
            size_t j = k - times.size();
            double n_pairs = j * (j-1.0) / 2.0;
            b += final_interval * n_pairs;
        }

        b *= 2.0;

        ln_prob_coal = RbConstants::LN2 * a + log(beta) * alpha + RbMath::lnGamma(a+alpha) - RbMath::lnGamma(alpha) - log(b+beta)*(a+alpha);
    }

    return ln_prob_coal;
}


double MultispeciesCoalescentInverseGammaPrior::drawNe( size_t index )
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    double u = RbStatistics::InverseGamma::rv(shape->getValue(), rate->getValue(), *rng);

    return u;
}



void MultispeciesCoalescentInverseGammaPrior::setShape(TypedDagNode<double>* s)
{

    removeParameter( shape );

    shape = s;

    addParameter( shape );
}



void MultispeciesCoalescentInverseGammaPrior::setRate(TypedDagNode<double>* r)
{

    removeParameter( rate );

    rate = r;

    addParameter( rate );
}


/** Swap a parameter of the distribution */
void MultispeciesCoalescentInverseGammaPrior::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == rate )
    {
        rate = static_cast<const TypedDagNode< double >* >( newP );
    }

    if ( oldP == shape )
    {
        shape = static_cast<const TypedDagNode< double >* >( newP );
    }
    AbstractMultispeciesCoalescent::swapParameterInternal(oldP, newP);

}

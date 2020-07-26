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
    // k is the number of entering lineages, so the log like is 0 if
    // there is only one lineages (as the probability of no coalescence
    // is equal to 1.0 in this case, as it is the only possible outcome)
    if ( k == 1 ) return 0.0;

    double alpha = shape->getValue();
    double beta = rate->getValue();

    double ln_prob_coal = 0.0;
    double current_time = begin_age;

    // Get the number of coalescences
    size_t n = times.size();

    // Get the rb term from Jones (2017)
    // We assume autosomal nuclear genes, so ploidy = 2
    double a = n;
    double r = -a * log(2.0);

    // We need to get the branch gamma term (gamma_b in Jones 2017)
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
        size_t j = k - n;
        double n_pairs = j * (j-1.0) / 2.0;
        b += final_interval * n_pairs;
    }

    // Divide by ploidy
    b /= 2.0;

    // Calculate the log gamma ratio
    double log_gamma_ratio = 0.0;
    for (size_t i=0; i<n; ++i)
    {
        log_gamma_ratio += log(alpha + i);
    }

    ln_prob_coal += r + (alpha * log(beta)) + log_gamma_ratio - ((alpha + a) * log(beta + b));

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

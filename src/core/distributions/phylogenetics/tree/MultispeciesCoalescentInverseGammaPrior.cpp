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
    a = 0.0;
    b = 0.0;
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
    a += n;

    for (size_t i=0; i<n; ++i)
    {
        // now we do the computation
        // t is the time between the previous and the current coalescences
        double t = times[i] - current_time;
        current_time = times[i];

        // get the number j of individuals we had before the current coalescence
        size_t j = k - i;
        double n_pairs = j * (j-1.0);

        b += t * n_pairs;
    }

    // compute the probability of no coalescent event in the final part of the branch
    // only do this if the branch is not the root branch
    if ( add_final_interval == true )
    {
        double final_interval = end_age - current_time;
        size_t j = k - times.size();
        double n_pairs = j * (j-1.0);
        b += final_interval * n_pairs;
    }


    // If we've gotten to the last node of the tree, then we can calculate the likelihood
    // for the entire gene tree given the species tree using the total number of gene
    // copies and the total coalescent rate over the entire genealogy
    double ln_prob_coal = 0.0;

    double num_tips = getNumberOfSpeciesTreeTips();

    if ( index == 2*(num_tips-1) )
    {
        ln_prob_coal += RbConstants::LN2 * a + log(beta) * alpha + RbMath::lnGamma(a+alpha) - RbMath::lnGamma(alpha) - log(b+beta)*(a+alpha);

        // Remember to reset the total coalescent rate and number of coalescent times so that
        // we don't just keep adding to them; we're done with them for this particular gene tree
        resetAB();
    }
    else
    {
        ln_prob_coal += 0.0;
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


double MultispeciesCoalescentInverseGammaPrior::getNumberOfSpeciesTreeTips( void )
{
    double num_species_tree_tips = species_tree->getValue().getNumberOfTips();

    return num_species_tree_tips;
}


void MultispeciesCoalescentInverseGammaPrior::resetAB( void )
{

    a = 0.0;
    b = 0.0;

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

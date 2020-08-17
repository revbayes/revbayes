#include <stddef.h>
#include <cmath>
#include <vector>

#include "ModelVector.h"
#include "MultispeciesCoalescentInverseGammaPrior.h"
#include "DistributionInverseGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "AbstractMultispeciesCoalescentGenewise.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

MultispeciesCoalescentInverseGammaPrior::MultispeciesCoalescentInverseGammaPrior(const TypedDagNode<Tree> *sp, RbVector< RbVector<Taxon> > t, size_t ngt) : AbstractMultispeciesCoalescentGenewise(sp, t, ngt)
{

}


MultispeciesCoalescentInverseGammaPrior::~MultispeciesCoalescentInverseGammaPrior()
{

}





MultispeciesCoalescentInverseGammaPrior* MultispeciesCoalescentInverseGammaPrior::clone( void ) const
{

    return new MultispeciesCoalescentInverseGammaPrior( *this );
}


double MultispeciesCoalescentInverseGammaPrior::computeLnCoalescentProbability(std::vector<size_t> k, const std::vector< std::vector<double> > &times, double begin_age, double end_age, size_t index, bool add_final_interval)
{
    // Index is the index of the species node

    // k is a vector holding the number of entering lineages per gene
    // So the log like is 0 for a particular gene i in this branch of the
    // species tree if k[i] = 1, as there is only one lineage and the
    // probability of no coalescence is equal to 1.0 in this case (as it
    // is the only possible outcome)

    double alpha = shape->getValue();
    double beta = rate->getValue();

    double current_time = begin_age;

    double ln_prob_coal = 0.0;

    // Initialize terms that are summed over all genes
    double a = 0; // q_b term in Jones (2017)
    double b = 0; // gamma_b term in Jones (2017)
    double log_r = 0; // log(r_b) term in Jones (2017)

    for (size_t i=0; i<num_gene_trees; i++)
    {
        // We only need to calculate terms if k > 1
        if ( k[i] > 1 )
        {
            // Get the number of coalescences
            size_t n = times[i].size();

            // Branch ploidy term (log)
            // We assume autosomal nuclear genes, so ploidy = 2
            log_r += -n * log(2.0);

            // Branch event term
            a += n;

            // Branch gamma term
            for (size_t m=0; m<n; ++m)
            {
                // Get the time t between the previous and the current coalescences
                double t = times[i][m] - current_time;
                current_time = times[i][m];

                // Get the number j of individuals we had before the current coalescence
                size_t j = k[i] - m;
                double n_pairs = j * (j-1.0) / 2.0;

                b += t * n_pairs;
            }

            // compute the probability of no coalescent event in the final part of the branch
            // only do this if the branch is not the root branch
            if ( add_final_interval == true )
            {
                double final_interval = end_age - current_time;
                size_t j = k[i] - n;
                double n_pairs = j * (j-1.0) / 2.0;
                b += final_interval * n_pairs;
            }
        }
    }

    // Get final branch gamma term by dividing sum by ploidy
    b /= 2.0;

    // Calculate the log gamma ratio
    double log_gamma_ratio = 0.0;
    for (size_t i=0; i<a; ++i)
    {
        log_gamma_ratio += log(alpha + i);
    }

    // Finally calculate the total log probability over all gene trees for this branch of the species tree
    ln_prob_coal += log_r + (alpha * log(beta)) - ((alpha + a) * log(beta + b)) + log_gamma_ratio;

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

    AbstractMultispeciesCoalescentGenewise::swapParameterInternal(oldP, newP);

}

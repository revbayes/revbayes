#include <cstddef>
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

MultispeciesCoalescentInverseGammaPrior::MultispeciesCoalescentInverseGammaPrior(const TypedDagNode<Tree> *sp, TypedDagNode<double>* sh, TypedDagNode<double>* sc, RbVector< RbVector<Taxon> > t, size_t ngt) : AbstractMultispeciesCoalescentGenewise(sp, t, ngt),
    shape( sh ),
    scale( sc )
{
    addParameter( shape );
    addParameter( scale );
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

    // k is a vector holding the number of entering lineages per gene.
    // So the log like is 0 for a particular gene i in this branch of the
    // species tree if k[i] = 1, as there is only one lineage and the
    // probability of no coalescence is equal to 1.0 in this case (as it
    // is the only possible outcome)

    double alpha = shape->getValue();
    double beta = scale->getValue();

    // Initialize terms that are summed over all genes
    int a = 0; // q_b term in Jones (2017); branchQ term in *BEAST
    double b = 0.0; // gamma_b term in Jones (2017); branch_gamma term in *BEAST
    //double log_r = 0.0; // log(r_b) term in Jones (2017); branchLogR term in *BEAST

    // std::cout << "-------------------" << std::endl;
    //
    // std::cout << "beta: " << beta << std::endl;
    //
    // std::cout << "start: " << begin_age << std::endl;
    // std::cout << "end: " << end_age << std::endl;

    for (size_t i=0; i<num_gene_trees; i++)
    {
        // std::cout << "k: " << k[i] << std::endl;

        // We only need to calculate terms if k > 1
        if ( k[i] > 1 )
        {
            double current_time = begin_age;
            double gene_b = 0.0;

            // Get the number of coalescences
            size_t n = times[i].size();
            // double nc = n;

            // std::cout << "n: " << n <<std::endl;

            // Branch ploidy term (log)
            // We assume autosomal nuclear genes, so ploidy = 2
            //log_r -= nc * RbConstants::LN2;

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

                gene_b += t * n_pairs;

                // std::cout << "t: " << t << std::endl;
                // std::cout << "n_pairs: " << n_pairs << std::endl;

            }

            // compute the probability of no coalescent event in the final part of the branch
            // only do this if the branch is not the root branch
            if ( add_final_interval == true )
            {
                double final_interval = end_age - current_time;
                size_t j = k[i] - n;
                double n_pairs = j * (j-1.0) / 2.0;
                gene_b += final_interval * n_pairs;

                // std::cout << "t: " << final_interval << std::endl;
                // std::cout << "n_pairs: " << n_pairs << std::endl;

            }

        b += gene_b * 2.0;
        // std::cout << "b for gene[" << i << "]: " << gene_b << std::endl;
        // std::cout << "current b: " << b << std::endl;
        }
    }

    // std::cout << "final a: " << a << std::endl;
    // std::cout << "final b: " << b << std::endl;

    // Calculate the log gamma ratio
    double log_gamma_ratio = 0.0;
    for (size_t i=0; i<a; ++i)
    {
        log_gamma_ratio += log(alpha + i);
    }

    // Finally calculate the total log probability over all gene trees for this branch of the species tree
    //double ln_prob_coal = log_r + (alpha * log(beta)) - ((alpha + a) * log(beta + b)) + log_gamma_ratio;

    if ((a == 0) && (b == 0.0))
    {
        return 0.0;
    }
    else
    {
        double ln_prob_coal = (a * RbConstants::LN2) + (alpha * log(beta)) - ((alpha + a) * log(beta + b)) + log_gamma_ratio;

        // std::cout << "ln prob coal: " << ln_prob_coal << "\n" << std::endl;

        return ln_prob_coal;
    }


}


double MultispeciesCoalescentInverseGammaPrior::drawNe( size_t index )
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    double u = RbStatistics::InverseGamma::rv(shape->getValue(), scale->getValue(), *rng);

    return u;
}


// double  MultispeciesCoalescentInverseGammaPrior::getShape(size_t index) const
// {

//     if ( shape != NULL )
//     {
//         return shape->getValue();
//     }
//     else
//     {
//         std::cerr << "Error: Null Pointers for shape." << std::endl;
//         exit(-1);
//     }
// }


// double  MultispeciesCoalescentInverseGammaPrior::getScale(size_t index) const
// {

//     if ( scale != NULL )
//     {
//         return scale->getValue();
//     }
//     else
//     {
//         std::cerr << "Error: Null Pointers for scale." << std::endl;
//         exit(-1);
//     }
// }


// void MultispeciesCoalescentInverseGammaPrior::setShape(TypedDagNode<double>* s)
// {

//     removeParameter( shape );

//     shape = s;

//     addParameter( shape );
// }


// void MultispeciesCoalescentInverseGammaPrior::setScale(TypedDagNode<double>* s)
// {

//     removeParameter( scale );

//     scale = s;

//     addParameter( scale );
// }


/** Swap a parameter of the distribution */
void MultispeciesCoalescentInverseGammaPrior::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == scale )
    {
        scale = static_cast<const TypedDagNode< double >* >( newP );
    }

    if ( oldP == shape )
    {
        shape = static_cast<const TypedDagNode< double >* >( newP );
    }

    AbstractMultispeciesCoalescentGenewise::swapParameterInternal(oldP, newP);

}

#include <cstddef>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/expint.hpp>
#include <iostream>

#include "MultispeciesCoalescentUniformPrior.h"
#include "DistributionUniform.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "AbstractMultispeciesCoalescent.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

MultispeciesCoalescentUniformPrior::MultispeciesCoalescentUniformPrior(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : AbstractMultispeciesCoalescent(sp, t)
{
    fn = 0.0;
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

    size_t n = times.size();
    double nc = n; // This is the number of coalescences, i.e., (m-n)

    for (size_t i=0; i<n; ++i)
    {
        // now we do the computation
        // a is the time between the previous and the current coalescences
        double a = times[i] - current_time;
        current_time = times[i];

        // get the number j of individuals we had before the current coalescence
        size_t j = k - i;
        double n_pairs = j * (j-1.0);

        fn += a * n_pairs;
    }

    // compute the probability of no coalescent event in the final part of the branch
    // only do this if the branch is not the root branch
    if ( add_final_interval == true )
    {
        double final_interval = end_age - current_time;
        size_t j = k - times.size();
        double n_pairs = j * (j-1.0);
        fn += final_interval * n_pairs;
    }

    // Now calculate the likelihood
    double ln_prob_coal = 0.0;

    // double ngc = getNumberOfGeneCopies();
    double integral_limit = 2.0 * fn / theta_max;
    double upper_incomplete_gamma = 0.0;

    if ( nc == 2 )
    {
        upper_incomplete_gamma = recursiveIncompleteGamma( nc-2.0, integral_limit );
    }
    else
    {
        double lower_incomplete_gamma = RbMath::incompleteGamma( integral_limit, nc-2.0, RbMath::lnGamma( nc-2.0 ) ) * RbMath::gamma( nc-2.0 );
        upper_incomplete_gamma = RbMath::gamma( nc-2.0 ) - lower_incomplete_gamma;
    }

    ln_prob_coal += RbConstants::LN2 + (( -nc+2 ) * log( fn )) + log( upper_incomplete_gamma ) - log( theta_max );


    // If we've gotten to the last node of the tree, then we can calculate the likelihood
    // for the entire gene tree given the species tree using the total number of gene
    // copies and the total coalescent rate over the entire genealogy

    // double num_tips = getNumberOfSpeciesTreeTips();

    // if ( index == 2*(num_tips-1) )
    // {
    //     double ngc = getNumberOfGeneCopies();
    //     double integral_limit = 2 * fn / theta_max;
    //
    //     // std::cout << "fn total: " << fn << std::endl;
    //     // std::cout << "ngc: " << ngc << std::endl;
    //     // std::cout << "integral limit: " << integral_limit << std::endl;
    //
    //     double upper_incomplete_gamma = 0.0;
    //     if ( ngc <= 2 )
    //     {
    //         upper_incomplete_gamma = recursiveIncompleteGamma( ngc-2.0, integral_limit );
    //     }
    //     else
    //     {
    //         double lower_incomplete_gamma = RbMath::incompleteGamma( integral_limit, ngc-2.0, RbMath::lnGamma( ngc-2.0 ) ) * RbMath::gamma( ngc-2.0 );
    //         upper_incomplete_gamma = RbMath::gamma( ngc-2.0 ) - lower_incomplete_gamma;
    //     }
    //
    //     ln_prob_coal += RbConstants::LN2 + (( -ngc+2 ) * log( fn )) + log( upper_incomplete_gamma ) - log( theta_max );
    //
    //     // Remember to reset the total coalescent rate so that we don't just keep adding to it
    //     // because we're now done with it for this particular gene tree
    //     resetFn();
    // }
    // // Otherwise we don't change the likelihood because we haven't finished adding up the
    // // total coalescent rate over the genealogy
    // else
    // {
    //     ln_prob_coal += 0.0;
    // }


    return ln_prob_coal;
}


/** Calculates the incomplete gamma recursively for shape parameter <= 0 **/
/** a is the shape parameter and x is the integral limit **/
double MultispeciesCoalescentUniformPrior::recursiveIncompleteGamma( double a, double x )
{

    if ( a == 0 )
    {
        // Base case
        // We need to get the exponential integral G(0,x) where x = the integral limit
        double ei = boost::math::expint( -x );
        double e1 = -ei;

        return e1;
    }
    else
    {
        return (-pow(x,a) * exp(-x) / a) + (1.0/a) * recursiveIncompleteGamma( a+1.0, x );
    }

}


double MultispeciesCoalescentUniformPrior::drawNe( size_t index )
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    double u = RbStatistics::Uniform::rv( 0, max_theta->getValue(), *rng);

    return u;

}


/** Get number of gene copies (i.e., number of tips in gene tree) **/
double MultispeciesCoalescentUniformPrior::getNumberOfGeneCopies( void )
{
    double nt = num_taxa;

    return nt;
}


double MultispeciesCoalescentUniformPrior::getNumberOfSpeciesTreeTips( void )
{
    double num_species_tree_tips = species_tree->getValue().getNumberOfTips();

    return num_species_tree_tips;
}


void MultispeciesCoalescentUniformPrior::setMaxTheta(TypedDagNode<double>* m)
{

    removeParameter( max_theta );

    max_theta = m;

    addParameter( max_theta );

}


void MultispeciesCoalescentUniformPrior::resetFn( void )
{

    fn = 0.0;

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

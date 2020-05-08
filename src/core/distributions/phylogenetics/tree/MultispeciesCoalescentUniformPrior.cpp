#include <stddef.h>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/expint.hpp>

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
    double nt = n;

    for (size_t i=0; i<n; ++i)
    {
        // now we do the computation
        //a is the time between the previous and the current coalescences
        double a = times[i] - current_time;
        current_time = times[i];

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

    double ln_prob_coal = RbConstants::LN2 - log( fn ) * (nt-2) - log( theta_max );

//    shape, rate/x
//    double lowerIncompleteGamma = RbMath::incompleteGamma( rate/x, shape, RbMath::lnGamma(shape) );
//    double gamma = RbMath::gamma(shape);

//    Gamma(n-2,2*fn/theta_max)

    // Now we need to deal with the incomplete gamma term

    // If the number of gene copies is 1, then there can be no coalescence event
    // and the probability is equal to 1.0 for the only possible event (no coalescence)
    if (n <= 1)
    {
        ln_prob_coal = 0.0;
    }
    else
    {
        double integral_limit = 2 * fn / theta_max;

        // When the shape term is 0 (as when n == 2), then we calculate
        // the upper incomplete gamma using an exponential integral instead
        // (see )
        if (n == 2)
        {
            if (integral_limit > 0)
            {
                double ei = boost::math::expint( -integral_limit );
                double e1 = -ei;

                ln_prob_coal += log(e1);
            }
            else
            {
                std::cerr << "The integral limit in dnMultiSpeciesCoalescentUniformPrior is negative." << std::endl;
            }

        }
        // Otherwise we calculate the incomplete gamma directly
        else
        {
            double lower_incomplete_gamma = RbMath::incompleteGamma( integral_limit, nt-2, RbMath::lnGamma(nt-2) ) * gamma(nt-2);
            double upper_incomplete_gamma = gamma(nt-2) - lower_incomplete_gamma;

            ln_prob_coal += log( upper_incomplete_gamma );
        }
    }

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

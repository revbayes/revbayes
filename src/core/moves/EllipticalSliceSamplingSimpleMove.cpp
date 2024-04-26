#include <cstddef>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "DagNode.h"
#include "DistributionNormal.h"
#include "NormalDistribution.h"
#include "EllipticalSliceSamplingSimpleMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TypedDagNode.h"
#include "AbstractMove.h"
#include "RbException.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

using namespace RevBayesCore;


/**
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
EllipticalSliceSamplingSimpleMove::EllipticalSliceSamplingSimpleMove( std::vector< StochasticNode<double> *> n, double window_, double weight_, bool t , bool fa) : AbstractMove( std::vector<DagNode*>(), weight_ ,t),
    variables( n ),
    window( window_ ),
    total_movement( 0.0 ),
    force_accept( fa )
{
    // Check window is valid (in [0,2*pi])
    if (window > RbConstants::TwoPI)
    {
        throw(RbException("mvEllipticalSliceSamplingSimple requires the window to be in the interval(0,2*pi]"));
    }

    // tell the base class to add the node
    for (size_t i=0; i<variables.size(); ++i)
    {
        assert( not variables[i]->isClamped() );

        NormalDistribution* dist = dynamic_cast<NormalDistribution *>( &variables[i]->getDistribution() );

        if (dist == NULL)
        {
            throw(RbException("mvEllipticalSliceSamplingSimple only works for vector of variables with Normal(0,sigma) Distribution"));
        }

        if (dist->getMean()->getValue() != 0.0)
        {
            throw(RbException("mvEllipticalSliceSamplingSimple only works for vector of variables with Normal(0,sigma) Distribution"));
        }

        addNode( variables[i] );
    }

}


/**
 * Basic destructor doing nothing.
 */
EllipticalSliceSamplingSimpleMove::~EllipticalSliceSamplingSimpleMove( void )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the EllipticalSliceSamplingSimpleMove.
 */
EllipticalSliceSamplingSimpleMove* EllipticalSliceSamplingSimpleMove::clone( void ) const
{
    return new EllipticalSliceSamplingSimpleMove( *this );
}


/**
 * Get moves' name of object
 *
 * \return The moves' name.
 */
const std::string& EllipticalSliceSamplingSimpleMove::getMoveName( void ) const
{
    static std::string name = "EllipticalSliceSamplingSimpleMove";

    return name;
}


double EllipticalSliceSamplingSimpleMove::getMoveTuningParameter( void ) const
{
    return window;
}


double unif()
{
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    return rng->uniform01();
}

namespace  {
/// This object allow computing the probability of the current point, and also store the variables's range
class slice_function
{
    std::vector<StochasticNode<double> *> variables;
    // std::vector<DagNode*> variables_nodes; // Copy of variables, as vector of DagNodes, to avoid repeatedly constructing this vector
    double lHeat;
    RbOrderedSet<DagNode*> affectedNodes;
    int num_evals;

public:

    int get_num_evals() const {return num_evals;}

    double operator()()
    {

        double lnPrior = 0.0;
        double lnLikelihood = 0.0;


        // Under Murray's sampler, we're drawing new values from this prior, and we don't need to compute their probability
        // 2. then we recompute the probability for all the affected nodes
        for (RbOrderedSet<DagNode*>::iterator it = affectedNodes.begin(); it != affectedNodes.end(); ++it)
        {
            if ( (*it)->isClamped() )
                lnLikelihood += (*it)->getLnProbability();
            else
                lnPrior += (*it)->getLnProbability();
        }

        // 3. exponentiate with the chain heat
        double lnPosterior = lHeat * lnLikelihood + lnPrior;

        return lnPosterior;
    }


    double operator()(std::vector<double> x)
    {

        ++num_evals;

        for (size_t i=0; i<variables.size(); ++i)
        {
            variables[i]->getValue() = x[i];

            // first we touch all the nodes
            // that will set the flags for recomputation
            variables[i]->touch();
        }


        double Pr_ = (*this)();

        // call accept for each node  --  automatically includes affected nodes
        for (size_t i=0; i<variables.size(); ++i)
        {
            variables[i]->keep();
        }

        // // keep the nodes, vectorized version avoids redundant calls
        // (*variables_nodes.begin())->keepVector(variables_nodes);

        return Pr_;
    }

    std::vector<double> current_value() const
    {
        std::vector<double> vals;

        for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
        {
            vals.push_back( (*it)->getValue() );
        }

        return vals;
    }

    slice_function(std::vector<StochasticNode<double> *> n, double l)
        :variables(n),
          lHeat(l),
          num_evals(0)
        // variables_nodes()
    {
        for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
        {
            (*it)->initiateGetAffectedNodes( affectedNodes );
        }
        // // Assemble std::vector<DagNode*> of the variables to enable vectorized calls to keep and initiateGetAffected
        // for (std::vector< StochasticNode<double> *>::const_iterator it = variables.begin(); it != variables.end(); it++)
        // {
        //     variables_nodes.push_back(*it);
        // }
        //
        // // Get all affected nodes of entire vector, vectorized call avoids redundant calls
        // (*variables_nodes.begin())->initiateGetAffectedNodesVector(affectedNodes, variables_nodes);
    }
};

}


void EllipticalSliceSamplingSimpleMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    if ( prHeat != 1 || pHeat != 1)
    {
        throw RbException("Elliptical slice sampling moves are invalid with heated priors or posteriors.");
    }

    slice_function lnL(variables, lHeat);

    // Current values
    std::vector<double> f = lnL.current_value();

    // Choose ellipse (Murray step 1)
    std::vector<double> nu = std::vector<double>(variables.size(),0.0);
    for (size_t i=0; i<variables.size(); ++i)
    {
        NormalDistribution* normal = dynamic_cast<NormalDistribution *>( &variables[i]->getDistribution() );

        const TypedDagNode<double>* sd = normal->getStDev();

        nu[i] = sd->getValue() * RbStatistics::Normal::rv(*GLOBAL_RNG);

    }

    // Log-likelihood threshold (Murray step 2)
    double L_f = lnL();

    assert( fabs(L_f - lnL(f)) < 1E-8 );

    double logy = L_f + log(unif());

    // Find initial boundaries and draw first value (Murray step 3)
    double theta = unif()*window;
    double L = theta - window;
    double R = theta;

    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);

    // Run slice loop (Murray steps 4-10, with automatic termination after 200 steps)
    std::vector<double> f_prime = std::vector<double>(variables.size(),0.0);
    size_t loop_iterations = 0;
    for (size_t i=0;i<1000;++i)
    {
        ++loop_iterations;
        // step 4
        for (size_t j = 0; j < f.size(); ++j)
        {
            f_prime[j] = f[j] * cos_theta + nu[j] * sin_theta;
        }

        // step 5 (step 6 is return value)
        double L_f_prime = lnL(f_prime);

        if (L_f_prime > logy)
        {
            break;
        }

        // steps 7-9 (10 being loop again)
        if (theta > 0.0)
        {
            R = theta;
        }
        else {
            L = theta;
        }
        // Draw next theta, compute sine/cosine
        theta = L + unif()*(R-L);
        sin_theta = std::sin(theta);
        cos_theta = std::cos(theta);
    }

    if ( loop_iterations == 1000 ) {
      if ( force_accept )
      {
        // After 1000 steps we should have hit theta = 0, so it should be possible to set f_prime = f and accept
        cos_theta = 1.0;

        for (size_t j = 0; j < f.size(); ++j)
        {
            f_prime[j] = f[j];
        }

        // step 5 (step 6 is return value)
        double L_f_prime = lnL(f_prime);

      }
      else
      {
        throw(RbException("mvEllipticalSliceSamplingSimple has iterated 1000 times without success. This should not be possible, something may be wrong."));
      }
    }

    total_movement += cos_theta;

    numPr += lnL.get_num_evals();

}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void EllipticalSliceSamplingSimpleMove::printSummary(std::ostream &o, bool current_period) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o << std::fixed;
    o << std::setprecision(4);

    // print the name
    const std::string &n = getMoveName();
    size_t spaces = 40 - (n.length() > 40 ? 40 : n.length());
    o << n;
    for (size_t i = 0; i < spaces; ++i) {
        o << " ";
    }
    o << " ";

    // print the DagNode name
    const std::string &dn_name = (*nodes.begin())->getName();
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i) {
        o << " ";
    }
    o << " ";

    // print the weight
    int w_length = 4;
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i) {
        o << " ";
    }
    o << weight;
    o << " ";

    // print the number of tries
    int t_length = 9;
    if (num_tried_total > 0) t_length -= (int)log10(num_tried_total);
    for (int i = 0; i < t_length; ++i) {
        o << " ";
    }
    o << num_tried_total;
    o << " ";

    // print the average distance moved
    o<<"\n";
    if (num_tried_current_period > 0)
    {
      o<<"  Ave. cos(angle(x1,x2)) = "<<total_movement/num_tried_current_period<<std::endl;
    }

    // print the average distance moved
    if (num_tried_current_period > 0)
    {
      o<<"  Ave. # of Pr evals = "<<double(numPr)/num_tried_current_period<<std::endl;
    }

    //    proposal->printParameterSummary( o );
    o<<"  window = "<<window<<std::endl;

    o << std::endl;

    o.setf(previousFlags);
    o.precision(previousPrecision);


}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void EllipticalSliceSamplingSimpleMove::resetMoveCounters( void )
{
    total_movement = 0.0;
    num_tried_current_period = 0;
    numPr = 0;
}


/**
 * Swap the current variables for a new one.
 *
 * \param[in]     oldN     The old variables that needs to be replaced.
 * \param[in]     newN     The new Revvariables.
 */
void EllipticalSliceSamplingSimpleMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<double> *>(newN);
            assert( not variables[i]->isClamped() );
        }
    }

}


void EllipticalSliceSamplingSimpleMove::setMoveTuningParameter(double tp)
{
    window = tp;
}


/**
 * Tune the move to accept the desired acceptance ratio.
 * We only compute the acceptance ratio here and delegate the call to the proposal.
 */
void EllipticalSliceSamplingSimpleMove::tune( void )
{

    if ( num_tried_current_period > 2 )
    {
        if (numPr/num_tried_current_period > 9)
        {
            window *= 0.9;
        }
    }

}

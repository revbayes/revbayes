#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

#include "DistributionNormal.h"
#include "UpDownSlideBactrianProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Proposal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
UpDownSlideBactrianProposal::UpDownSlideBactrianProposal( double l ) : Proposal(),
    up_scalar_variables(  ),
    up_vector_variables(  ),
    down_scalar_variables(  ),
    down_vector_variables(  ),
    lambda( l )
{

}


/**
 * Add a variable.
 */
void UpDownSlideBactrianProposal::addVariable( StochasticNode<double> *v, bool up)
{


    // add it to the nodes vector
    addNode( v );

    // add it to my internal vector
    if ( up == true )
    {
        up_scalar_variables.push_back( v );
    }
    else
    {
        down_scalar_variables.push_back( v );
    }

}


/**
 * Add a variable.
 */
void UpDownSlideBactrianProposal::addVariable( StochasticNode<RbVector<double> > *v, bool up)
{


    // add it to the nodes vector
    addNode( v );

    // add it to my internal vector
    if ( up == true )
    {
        up_vector_variables.push_back( v );
    }
    else
    {
        down_vector_variables.push_back( v );
    }

}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void UpDownSlideBactrianProposal::cleanProposal( void )
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
UpDownSlideBactrianProposal* UpDownSlideBactrianProposal::clone( void ) const
{

    return new UpDownSlideBactrianProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& UpDownSlideBactrianProposal::getProposalName( void ) const
{
    static std::string name = "upDownSlideBactrian";

    return name;
}


double UpDownSlideBactrianProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and Slides the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double UpDownSlideBactrianProposal::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    double u = rng->uniform01();
    double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * RbConstants::BACT_SD;

    // See Yang and Rodriguez (2013) SI eqns 19 and 20
    // Currently hard-coding m = 0.95
    if ( u < 0.5 )
    {
        delta += 0.95;
    }
    else
    {
        delta -= 0.95;
    }

    delta = ( lambda * delta );

    stored_up_scalar_values = std::vector<double>(up_scalar_variables.size(), 0.0);
    // Slide all the single variable up
    for (size_t i=0; i<up_scalar_variables.size(); ++i)
    {
        stored_up_scalar_values[i] = up_scalar_variables[i]->getValue();
        up_scalar_variables[i]->getValue() += delta;
    }

    stored_up_vector_values = std::vector<std::vector<double> >(up_vector_variables.size(), std::vector<double>());
    // Slide all the vector variables up
    for (size_t i=0; i<up_vector_variables.size(); ++i)
    {
        RbVector<double> &v = up_vector_variables[i]->getValue();

        stored_up_vector_values[i] = std::vector<double>(v.size(), 0.0);
        for (size_t j=0; j<v.size(); ++j)
        {
            stored_up_vector_values[i][j] = v[j];
            v[j] += delta;
        }
    }

    stored_down_scalar_values = std::vector<double>(down_scalar_variables.size(), 0.0);
    // Slide all the single variable down
    for (size_t i=0; i<down_scalar_variables.size(); ++i)
    {
        stored_down_scalar_values[i] = down_scalar_variables[i]->getValue();
        down_scalar_variables[i]->getValue() -= delta;
    }

    stored_down_vector_values = std::vector<std::vector<double> >(down_vector_variables.size(), std::vector<double>());
    // Slide all the vector variables down
    for (size_t i=0; i<down_vector_variables.size(); ++i)
    {
        RbVector<double> &v = down_vector_variables[i]->getValue();

        stored_down_vector_values[i] = std::vector<double>(v.size(), 0.0);
        for (size_t j=0; j<v.size(); ++j)
        {
            stored_down_vector_values[i][j] = v[j];
            v[j] -= delta;
        }
    }

    // compute the Hastings ratio (symmetric move)
    double lnHastingsratio = 0.0;

    return lnHastingsratio;
}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void UpDownSlideBactrianProposal::prepareProposal( void )
{

}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void UpDownSlideBactrianProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }

}


/**
 * Remove a variable.
 */
void UpDownSlideBactrianProposal::removeVariable( StochasticNode<double> *v, bool up)
{


    // add it to the nodes vector
    removeNode( v );

    // add it to my internal vector
    if ( up == true )
    {
        up_scalar_variables.erase(std::remove(up_scalar_variables.begin(), up_scalar_variables.end(), v), up_scalar_variables.end());
    }
    else
    {
        down_scalar_variables.erase(std::remove(down_scalar_variables.begin(), down_scalar_variables.end(), v), down_scalar_variables.end());
    }

}


/**
 * Remove a variable.
 */
void UpDownSlideBactrianProposal::removeVariable( StochasticNode<RbVector<double> > *v, bool up)
{


    // add it to the nodes vector
    removeNode( v );

    // add it to my internal vector
    if ( up == true )
    {
        up_vector_variables.erase(std::remove(up_vector_variables.begin(), up_vector_variables.end(), v), up_vector_variables.end());
    }
    else
    {
        down_vector_variables.erase(std::remove(down_vector_variables.begin(), down_vector_variables.end(), v), down_vector_variables.end());
    }

}

/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void UpDownSlideBactrianProposal::undoProposal( void )
{

    // Slide all the single variable up
    for (size_t i=0; i<up_scalar_variables.size(); ++i)
    {
        up_scalar_variables[i]->getValue() = stored_up_scalar_values[i];
    }

    // Slide all the vector variables up
    for (size_t i=0; i<up_vector_variables.size(); ++i)
    {
        RbVector<double> &v = up_vector_variables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] = stored_up_vector_values[i][j];
        }
    }

    // Slide all the single variable down
    for (size_t i=0; i<down_scalar_variables.size(); ++i)
    {
        down_scalar_variables[i]->getValue() = stored_down_scalar_values[i];
    }

    // Slide all the vector variables down
    for (size_t i=0; i<down_vector_variables.size(); ++i)
    {
        RbVector<double> &v = down_vector_variables[i]->getValue();
        for (size_t j=0; j<v.size(); ++j)
        {
            v[j] = stored_down_vector_values[i][j];
        }
    }

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void UpDownSlideBactrianProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    for (size_t i = 0; i < up_scalar_variables.size(); ++i)
    {
        if ( up_scalar_variables[i] == oldN )
        {
            up_scalar_variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

    for (size_t i = 0; i < up_vector_variables.size(); ++i)
    {
        if ( up_vector_variables[i] == oldN )
        {
            up_vector_variables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }

    for (size_t i = 0; i < down_scalar_variables.size(); ++i)
    {
        if ( down_scalar_variables[i] == oldN )
        {
            down_scalar_variables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

    for (size_t i = 0; i < down_vector_variables.size(); ++i)
    {
        if ( down_vector_variables[i] == oldN )
        {
            down_vector_variables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }

}


void UpDownSlideBactrianProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void UpDownSlideBactrianProposal::tune( double rate )
{

    if ( rate > 0.44 )
    {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        lambda /= (2.0 - rate/0.44 );
    }

    if ( lambda > 1 ) {
        lambda = fmin(10000, lambda);
    } else {
        lambda = fmax(1/10000, lambda);
    }

}

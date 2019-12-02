#include "DistributionNormal.h"
#include "MirrorMultiplierProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MirrorMultiplierProposal::MirrorMultiplierProposal( StochasticNode<double> *n, double l, double m0, double n0, double c0, double m, bool a, double p) : Proposal(p),
    variable( n ),
    storedValue( 0.0 ),
    lambda( l ),
    waitBeforeLearning( n0 ),
    waitBeforeUsing ( c0 ),
    maxUpdates ( m ),
    muHat ( m0 ),
    sigmaHat ( l ),
    learnedMuHat( 0.0 ),
    learnedSigmaHat ( l ),
    updates ( 0 ),
    nTried ( 0 ),
    adaptOnly ( a )
{
    // tell the base class to add the node
    addNode( variable );
// std::cout << "waitBeforeUsing " << waitBeforeUsing << "; waitBeforeLearning " << waitBeforeLearning << std::endl;
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MirrorMultiplierProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MirrorMultiplierProposal* MirrorMultiplierProposal::clone( void ) const
{

    return new MirrorMultiplierProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MirrorMultiplierProposal::getProposalName( void ) const
{
    static std::string name = "MirrorMultiplier";

    return name;
}


double MirrorMultiplierProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and slides the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MirrorMultiplierProposal::doProposal( void )
{
    ++nTried;

    double &val = variable->getValue();

    // copy value
    storedValue = val;

    // Perform move on log-scale
    val = log(val);

    if ( nTried > waitBeforeLearning && updates <= maxUpdates )
    {
      ++updates;
      double learnedMuHatOld = learnedMuHat;

      // Update average
      double n = double(updates);
      learnedMuHat = 1/n * val + (n - 1)/n * learnedMuHat;

      // Update sd
      double ssq = learnedSigmaHat*learnedSigmaHat;
      learnedSigmaHat = sqrt(ssq + ( ((val - learnedMuHatOld) * (val - learnedMuHat)) - ssq)/n);

      // std::cout << "learnedMuHat " << learnedMuHat << "; learnedSigmaHat " << learnedSigmaHat << std::endl;
    }

    double delta = RbStatistics::Normal::rv(*GLOBAL_RNG) * lambda * sigmaHat;
    double fwdMean = 2.0 * muHat - val;
    double newVal = fwdMean + delta;

    double lnHastingsratio = -RbStatistics::Normal::lnPdf(fwdMean, lambda * sigmaHat, newVal);
    double revMean = 2.0 * muHat - newVal;
    lnHastingsratio += RbStatistics::Normal::lnPdf(revMean, lambda * sigmaHat, val);
    lnHastingsratio += -val + newVal;

    val = exp(newVal);

    return lnHastingsratio;
}


/**
 *
 */
void MirrorMultiplierProposal::prepareProposal( void )
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
void MirrorMultiplierProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "sigma = ";
    if (name_only == false)
    {
        o << lambda;
    }

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MirrorMultiplierProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( new double(storedValue) );

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MirrorMultiplierProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    variable = static_cast< ContinuousStochasticNode* >(newN) ;

}


void MirrorMultiplierProposal::setProposalTuningParameter(double tp)
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
void MirrorMultiplierProposal::tune( double rate )
{
  // Early moves will be difficult until mean/sd are in ballpark, but until then
  // we should avoid responding to early rejections by shrinking proposal covariance
  // We also want to ensure that the move only updates when tuning is allowed,
  // Hence we only update the usable mean/variance here
  if ( nTried > waitBeforeUsing )
  {
    muHat = learnedMuHat;
    sigmaHat = learnedSigmaHat;

    // In case user wants to learn mean/variance but not adjust proposal variance further
    if ( !adaptOnly )
    {
      double p = this->targetAcceptanceRate;
      if ( rate > p )
      {
          lambda *= 1.1;
      }
      else
      {
          lambda /= 1.1;
      }

      lambda = fmin(10000, lambda);
    }
  }
    // std::cout << "Updating parameters, muHat " << muHat << "; sigmaHat " << sigmaHat << "; sigma " << lambda << std::endl;
}

#include "MPQRateMatrixProposal.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RateMatrix_MPQ.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"


#define MIN_FREQ    10e-4
#define A           0
#define C           1
#define G           2
#define T           3

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MPQRateMatrixProposal::MPQRateMatrixProposal( StochasticNode<RateGenerator> *n ) : Proposal(),
variable( n )
{
    
    rev_alpha_er = 1.0;
    rev_alpha_pi = 1.0;
    non_rev_alpha = 1.0;
    rj_alpha = 1.0;
    
    W.resize( 6 );

    // tell the base class to add the node
    addNode( variable );
}


MPQRateMatrixProposal::MPQRateMatrixProposal( const MPQRateMatrixProposal& p ) : Proposal( p ),
variable( p.variable ),
rev_alpha_pi( p.rev_alpha_pi ),
rev_alpha_er( p.rev_alpha_er),
non_rev_alpha( p.non_rev_alpha ),
rj_alpha( p.rj_alpha ),
stored_Q( p.stored_Q ),
W( p.W )
{
        
    // tell the base class to add the node
    addNode( variable );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MPQRateMatrixProposal::cleanProposal( void ) {

    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MPQRateMatrixProposal* MPQRateMatrixProposal::clone( void ) const {
    
    return new MPQRateMatrixProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MPQRateMatrixProposal::getProposalName( void ) const {

    static std::string name = "MPQRateMatrixProposal";
    return name;
}


double MPQRateMatrixProposal::getProposalTuningParameter( void ) const
{
    return rev_alpha_pi;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MPQRateMatrixProposal::doProposal( void ) {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    double lnProb = 0.0;
    if (Q.getIsReversible() == true)
        {
        // update GTR model
        double u = rng->uniform01();
        if (u < 0.55)
            {
            double u2 = rng->uniform01();
            bool all = u2 < 0.25;
            lnProb = updateStationaryFrequencies(all);
            }
        else if (u >= 0.55 && u < 0.90)
            {
            double u2 = rng->uniform01();
            bool all = u2 < 0.166;
            lnProb = updateExchangabilityRates(all);
            }
        else
            {
            lnProb = updateToNonReversible();
            }
        }
    else
        {
        // update non-reversible model
        double u = rng->uniform01();
        if (u < 0.90)
            {
            double u2 = rng->uniform01();
            bool all = u2 < 0.1;
            lnProb = updateNonreversibleRates(all);
            }
        else
            {
            lnProb = updateToReversible();
            }
        }
        
    
    return lnProb;
    
}


/**
 *
 */
void MPQRateMatrixProposal::prepareProposal( void ) {
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void MPQRateMatrixProposal::printParameterSummary(std::ostream &o, bool name_only) const {
    
    o << "lambda = ";
    if (name_only == false)
    {
        o << rev_alpha_pi;
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MPQRateMatrixProposal::undoProposal( void ) {
    
    RateMatrix_MPQ& v = static_cast<RateMatrix_MPQ&>( variable->getValue() );

    v = stored_Q;
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MPQRateMatrixProposal::swapNodeInternal(DagNode *oldN, DagNode *newN) {
    
    variable = static_cast< StochasticNode<RateGenerator>* >(newN) ;
}


void MPQRateMatrixProposal::setProposalTuningParameter(double tp)
{
    rev_alpha_pi = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void MPQRateMatrixProposal::tune( double rate ) {
    
//    if ( rate > 0.44 )
//        {
//        rev_alpha_pi *= (1.0 + ((rate-0.44)/0.56) );
//        }
//    else
//        {
//        rev_alpha_pi /= (2.0 - rate/0.44 );
//        }
//    rev_alpha_pi = fmin(10000, rev_alpha_pi);
}


double MPQRateMatrixProposal::updateExchangabilityRates(bool all) {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());

    if (Q.getIsReversible() == false)
        throw RbException("Can only update the exchangability rates for time reversible models");
            
    // store old values
    stored_Q = Q;

    if ( all = true )
    {
        // update and return log probability of Hastings ratio
        return Q.updateExchangeabilityRates(rng, 100.0, 0.5);
    }
    else
    {
        // update and return log probability of Hastings ratio
        return Q.updateExchangeabilityRatesSingle(rng, 50);
    }
}



double MPQRateMatrixProposal::updateNonreversibleRates(bool all) {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == true)
        throw RbException("Can only update the 12 rates (directly) for non-reversible models");

    // store old values
    stored_Q = Q;

    if ( all == true )
    {
        // update and return log probability of Hastings ratio
        return Q.updateNonReversibleRates(rng, 500.0, 0.5);
    }
    else
    {
        // update and return log probability of Hastings ratio
        return Q.updateNonReversibleRatesSingle(rng, 50.0);
    }
}

double MPQRateMatrixProposal::updateToNonReversible(void) {
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == false)
        throw RbException("Can only move to a non-reversible model from a reversible one");

    // store old values
    stored_Q = Q;

    // calculate the Jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double wCG = piC * stored_Q(C,G).get_d();
    double wCT = piC * stored_Q(C,T).get_d();
    double wGT = piG * stored_Q(G,T).get_d();
    double lnJacobian = (log(64.0) + log(wCG) + log(wCT) + log(wGT));
    lnJacobian -= log(piC) + 2.0 * log(piG);

    // polyhedron density for forward move
    Q.calculateWeights(W);
    Vector randomPoint;
    double lnRv = -poly.lnProbabilityForward(W, randomPoint);
    
    // update rates (off diagonal components)
    mpq_class u1 = randomPoint.getX();
    mpq_class u2 = randomPoint.getY();
    mpq_class u3 = randomPoint.getZ();
#   if 1
    Q.nonreversibilize(u1, u2, u3);
#   else
    Q.setIsReversible(false);
    Q(A,C) = stored_Q(A,C) + (pi[C]/pi[A]) * stored_Q(C,G) * (2 * u1 - 1) + (pi[C]/pi[A]) * stored_Q(C,T) * (2 * u2 - 1);
    Q(A,G) = stored_Q(A,G) - (pi[C]/pi[A]) * stored_Q(C,G) * (2 * u1 - 1) + (pi[G]/pi[A]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(A,T) = stored_Q(A,T) - (pi[C]/pi[A]) * stored_Q(C,T) * (2 * u2 - 1) - (pi[G]/pi[A]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(C,G) = 2 * stored_Q(C,G) * u1;
    Q(C,T) = 2 * stored_Q(C,T) * u2;
    Q(G,T) = 2 * stored_Q(G,T) * u3;
    Q(C,A) = (pi[A]/pi[C]) * stored_Q(A,C) - stored_Q(C,G) * (2 * u1 - 1) - stored_Q(C,T) * (2 * u2 - 1);
    Q(G,A) = (pi[A]/pi[G]) * stored_Q(A,G) + (pi[C]/pi[G]) * stored_Q(C,G) * (2 * u1 - 1) - stored_Q(G,T) * (2 * u3 - 1);
    Q(T,A) = (pi[A]/pi[T]) * stored_Q(A,T) + (pi[C]/pi[T]) * stored_Q(C,T) * (2 * u2 - 1) + (pi[G]/pi[T]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(G,C) = 2 * (pi[C]/pi[G]) * stored_Q(C,G) * (1 - u1);
    Q(T,C) = 2 * (pi[C]/pi[T]) * stored_Q(C,T) * (1 - u2);
    Q(T,G) = 2 * (pi[G]/pi[T]) * stored_Q(G,T) * (1 - u3);
    
    // update the diagonal components of the rate matrix
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
#   endif
    
    return lnRv + lnJacobian;
}

double MPQRateMatrixProposal::updateToReversible(void) {
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == true)
        throw RbException("Can only move to a time reversible model from a non-reversible one");

    // store old values
    stored_Q = Q;

    // average rates with same stationary frequencies
#   if 1
    Q.reversibilize();
#   else
    Q.setIsReversible(true);
    std::vector<mpq_class>& pi = Q.getPi();
    mpq_class averageRate;
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                Q(i,j) = (pi[i] * stored_Q(i,j) + pi[j] * stored_Q(j,i)) / (2 * pi[i]);
                sum += Q(i,j);
                averageRate += pi[i] * Q(i,j);
                }
            }
        Q(i,i) = -sum;
        }
        
    Q.setExchangeabilityRates();
    
    // make certain average rate is one
    if (averageRate != 1)
        {
        mpq_class factor = 1 / averageRate;
        std::cout << "Warning: the average rate should be one in updateToReversible" << std::endl;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                Q(i,j) *= factor;
        }
#   endif

    // jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double piT = pi[T].get_d();
    double wCG = piC * stored_Q(C,G).get_d();
    double wGC = piG * stored_Q(G,C).get_d();
    double wCT = piC * stored_Q(C,T).get_d();
    double wTC = piT * stored_Q(T,C).get_d();
    double wGT = piG * stored_Q(G,T).get_d();
    double wTG = piT * stored_Q(T,G).get_d();
    //double lnJacobian = (piC * piG * piG) / (8.0 * (wCG + wGC) * (wCT + wTC) * (wGT + wTG));
    double lnJacobian = log(piC) + 2.0 * log(piG);
    lnJacobian -= (log(8.0) + log(wCG + wGC) + log(wCT + wTC) + log(wGT + wTG));
    
    // polyhedron parameters
    Q.calculateWeights(W);
    
    // figure out where in the (u1,u2,u3) space the non-reversible model lives
    mpq_class u1 = stored_Q(C,G) / (2 * Q(C,G));
    mpq_class u2 = stored_Q(C,T) / (2 * Q(C,T));
    mpq_class u3 = stored_Q(G,T) / (2 * Q(G,T));
    Vector pt(u1, u2, u3);
    double lnRv = poly.lnProbabilityReverse(W, pt);

    return lnRv + lnJacobian;
}

double MPQRateMatrixProposal::updateStationaryFrequencies(bool all) {

    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == false)
        throw RbException("Can only update the stationary frequencies (directly) for time reversible models");

    // store old values
    stored_Q = Q;
    
    if ( all == true )
    {
        // update and return log probability of Hastings ratio
        return Q.updateStationaryFrequencies(rng, 100.0, 0.5);
    }
    else
    {
        // update and return log probability of Hastings ratio
        return Q.updateStationaryFrequenciesSingle(rng, 50.0);
    }

}


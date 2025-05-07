#include "CorrelationMatrixProposal.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
CorrelationMatrixProposal::CorrelationMatrixProposal( StochasticNode<MatrixReal> *n, RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s, double l ) : Proposal(),
variable( n ),
sigma( s ),
lambda( l ),
storedValue( 0.0 )
{

    // tell the base class to add the node
    addNode( variable );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void CorrelationMatrixProposal::cleanProposal( void ) {

    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
CorrelationMatrixProposal* CorrelationMatrixProposal::clone( void ) const {
    
    return new CorrelationMatrixProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& CorrelationMatrixProposal::getProposalName( void ) const {

    static std::string name = "CorrelationMatrixProposal";
    return name;
}


double CorrelationMatrixProposal::getProposalTuningParameter( void ) const
{
    return lambda;
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
double CorrelationMatrixProposal::doProposal( void ) {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    MatrixReal& v = variable->getValue();
    //const RbVector<double>& s = sigma->getValue();
    
    // choose an entry (a,b) in the matrix to change
    idxA = size_t( rng->uniform01() * v.getNumberOfRows() );
    idxB = size_t( rng->uniform01() * v.getNumberOfColumns() );
    
    // copy value
    storedValue = v[idxA][idxB];
    
    // find the range of values over which r_ab is positive definite
    MatrixReal m1 = MatrixReal(v);
    MatrixReal m2 = MatrixReal(v);
    MatrixReal m3 = MatrixReal(v);
    m1[idxA][idxB] = 1.0;
    m1[idxB][idxA] = 1.0;
    m2[idxA][idxB] = -1.0;
    m2[idxB][idxA] = -1.0;
    m3[idxA][idxB] = 0.0;
    m3[idxB][idxA] = 0.0;
    double det1 = m1.getDet();
    double det2 = m2.getDet();
    double det3 = m3.getDet();
    double a = 0.5 * (det1 + det2 - 2.0 * det3);
    double b = 0.5 * (det1 - det2);
    double c = det3;
    double root1 = (-b + sqrt(b*b-4.0*a*c)) / (2.0 * a);
    double root2 = (-b - sqrt(b*b-4.0*a*c)) / (2.0 * a);
    double lLim = -1.0, uLim = 1.0;
    if (root1 < root2)
        {
        lLim = root1;
        uLim = root2;
        }
    else
        {
        lLim = root2;
        uLim = root1;
        }
    if (lLim < -1.0 || uLim > 1.0)
        {
        // better bail out with an error here
        }
    
    
    std::cout << det1 << " " << det2 << " " << det3 << std::endl;
    std::cout << lLim << "-" << uLim << std::endl;
    
    
    
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double newVal = lLim + (rng->uniform01() * (uLim-lLim));
    std::cout << "newValue = " << newVal << std::endl;
    v[idxA][idxB] = newVal;
    v[idxB][idxA] = newVal;

    // check that the resulting variance-covariance matrix is positive definite

    std::cout << v.isPositiveDefinite() << std::endl;
    std::cout << v << std::endl;
    exit(1);

    // this code is commented out because it the exit(1) command ensures it is never called
    // variable->addTouchedElementIndex(idxA);
    // variable->addTouchedElementIndex(idxB);
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
    
#if 0
    std::vector<double> weights = std::vector<double>(new_brothers.size(), 0.0);
    double sumOfWeights = 0.0;
    for (size_t i = 0; i<new_brothers.size(); ++i)
    {
        // get the new brother
        TopologyNode* newBro = new_brothers[i];
        
        // do the proposal
        TopologyNode *newGrandparent = pruneAndRegraft(&brother, newBro, parent, grandparent);
        
        // flag for likelihood recomputation
        variable->touch();
        
        // compute the likelihood of the new value
        double priorRatio = variable->getLnProbability();
        double likelihoodRatio = 0.0;
        for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
        {
            likelihoodRatio += (*it)->getLnProbability();
        }
        weights[i] = exp(priorRatio + likelihoodRatio + offset);
        sumOfWeights += weights[i];
        
        // undo proposal
        pruneAndRegraft(newBro, &brother, parent, *newGrandparent);
        
        // restore the previous likelihoods;
        variable->restore();
    }
    

#endif
}


/**
 *
 */
void CorrelationMatrixProposal::prepareProposal( void ) {
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void CorrelationMatrixProposal::printParameterSummary(std::ostream &o, bool name_only) const {
    
    o << "lambda = ";
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
void CorrelationMatrixProposal::undoProposal( void ) {
    
    MatrixReal& v = variable->getValue();
    v[idxA][idxB] = storedValue;
    v[idxB][idxA] = storedValue;
    variable->clearTouchedElementIndices();
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void CorrelationMatrixProposal::swapNodeInternal(DagNode *oldN, DagNode *newN) {
    
    variable = static_cast< StochasticNode<MatrixReal>* >(newN) ;
}


void CorrelationMatrixProposal::setProposalTuningParameter(double tp)
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
void CorrelationMatrixProposal::tune( double rate ) {
    
    if ( rate > 0.44 )
        {
        lambda *= (1.0 + ((rate-0.44)/0.56) );
        }
    else
        {
        lambda /= (2.0 - rate/0.44 );
        }
    lambda = fmin(10000, lambda);
}


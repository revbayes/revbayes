#include "DirichletSimplexProposal.h"

#include <cmath>
#include <iostream>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "RbStatisticsHelper.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
DirichletSimplexProposal::DirichletSimplexProposal( StochasticNode<Simplex> *n, double a, size_t nc, double o, double k /*=0.0*/, double p) : Proposal(p),
    variable( n ),
    storedValue( 0.0 ),
    alpha( a ),
    nCategories( nc ),
    offset( o ),
    kappa( k )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void DirichletSimplexProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
DirichletSimplexProposal* DirichletSimplexProposal::clone( void ) const
{
    
    return new DirichletSimplexProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& DirichletSimplexProposal::getProposalName( void ) const
{
    static std::string name = "DirichletSimplex";
    
    return name;
}


double DirichletSimplexProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A Dirichlet-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Dirichlet(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double DirichletSimplexProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // store the value
    storedValue = variable->getValue();
    
    const RbVector<double>& curVal = variable->getValue();
    RbVector<double> newVal = curVal;
    size_t           n      = curVal.size();
    
    /* We update the simplex values by proposing new values from a Dirichlet centered
     on the current values. The i-th parameter of the Dirichlet is the i-th value
     in the simplex multiplied by a parameter (alpha0, AKA tuning) that controls the
     variance of the Dirichlet. We implement two cases of this general move. In one
     case, all of the elements of the simplex are targeted for update (n == k). In the
     other, more complicated, case a subset of the elements of the simplex are updated
     (k < n). Here, we construct a smaller simplex with k+1 elements. The first k of the
     elements are the values from the full simplex that were targeted for update. The last
     element of the smaller simplex accumulates the probabilities of all of the simplex
     values in the full simplex that were not targeted for update. We then update the
     small simplex by centering a Dirichlet on the small simplex. The values for those elements
     in the full simplex that were not targeted for update are all changed proportionally.
     This means that we need to calculate the Jacobian for the Hastings ratio in this case. */
    double lnProposalRatio = 0.0;
    if ( nCategories > n )
    {
        // we can't update more elements than there are elements in the simplex
        throw RbException( "Attempting to update too many simplex variables" );
    }
    else if ( nCategories < 1 )
    {
        // we also can't update 0 or a negative number of elements
        throw RbException( "Attempting to update too few simplex variables" );
    }
    else if ( nCategories < n )
    {
        // we update a subset of the elements in the full simplex
        // pick k values at random, producing a map from the index in the full vector (curVal) to
        // the index in the reduced vector (x, alpha, and z)
        std::vector<size_t> indicesToUpdate;
        std::vector<size_t> tmpV;
        for (size_t i=0; i<n; i++)
        {
            tmpV.push_back(i);
        }
        
        RbStatistics::Helper::randomlySelectFromVectorWithoutReplacement<size_t>(tmpV, indicesToUpdate, nCategories, *rng);
        std::map<size_t,size_t> mapper;
        for (size_t i=0; i<indicesToUpdate.size(); i++)
        {
            mapper.insert( std::make_pair(indicesToUpdate[i], i) );
        }
        
        // set up the vectors
        std::vector<double> x(indicesToUpdate.size()+1, 0.0);
        std::vector<double> kappaV(indicesToUpdate.size()+1, 0.0);
        std::vector<double> alphaForward(indicesToUpdate.size()+1, 0.0);
        std::vector<double> alphaReverse(indicesToUpdate.size()+1, 0.0);
        std::vector<double> z(indicesToUpdate.size()+1, 0.0);
        for (size_t i=0; i<n; i++)
        {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
            {
                x[it->second] += curVal[it->first];
                kappaV[it->second] += kappa;
            }
            else
            {
                x[x.size()-1] += curVal[i];
                kappaV[kappaV.size()-1] += kappa;
            }
        }
        
        for (size_t i=0; i<x.size(); i++)
        {
            alphaForward[i] = (x[i]+offset) * alpha + kappaV[i];
        }
        
        // draw a new value for the reduced vector
        z = RbStatistics::Dirichlet::rv( alphaForward, *rng );
        
        // fill in the Dirichlet parameters for the reverse probability calculations
        for (size_t i=0; i<z.size(); i++)
        {
            alphaReverse[i] = (z[i]+offset) * alpha + kappaV[i];
        }
        
        // fill in the full vector
        double factor = z[z.size()-1] / x[x.size()-1];
        for (size_t i=0; i<n; i++)
        {
            std::map<size_t,size_t>::iterator it = mapper.find(i);
            if (it != mapper.end())
            {
                newVal[i] = z[it->second];
            }
            else
            {
                newVal[i] = curVal[i] * factor;
            }
            
            // test for 0-values
            if ( newVal[i] < 1E-100)
            {
                return RbConstants::Double::neginf;
            }
        }
        
        try
        {
            // Hastings ratio
            lnProposalRatio  = RbStatistics::Dirichlet::lnPdf(alphaReverse, x) - RbStatistics::Dirichlet::lnPdf(alphaForward, z); // Hastings Ratio
            lnProposalRatio += (n - nCategories - 1) * log(factor); // Jacobian
        }
		catch (RbException &e)
        {
            lnProposalRatio = RbConstants::Double::neginf;
        }
        
        
    }
    else
    {
        // we update all of the elements in the vector
        // first, we get the parameters of the Dirichlet for the forward move
        std::vector<double> alphaForward(curVal.size());
        for (size_t i=0; i<curVal.size(); i++)
        {
            alphaForward[i] = (curVal[i]+offset) * alpha + kappa;
            // we need to check for 0 values
            if (alphaForward[i] < 1E-100)
            {
                // very low proposal probability which will hopefully result into a rejected proposal
                return RbConstants::Double::neginf;
            }
        }
        
        // then, we propose new values
        newVal = RbStatistics::Dirichlet::rv( alphaForward, *rng );
        
        // and calculate the Dirichlet parameters for the (imagined) reverse move
        std::vector<double> alphaReverse(newVal.size());
        for (size_t i=0; i<curVal.size(); i++)
        {
            alphaReverse[i] = (newVal[i]+offset) * alpha + kappa;
            // we need to check for 0 values
            if (alphaReverse[i] < 1E-100)
            {
                // very low proposal probability which will hopefully result into a rejected proposal
                return RbConstants::Double::neginf;
            }
        }
        
                
        try
        {
            // finally, we calculate the log of the Hastings ratio
            lnProposalRatio = RbStatistics::Dirichlet::lnPdf(alphaReverse, curVal) - RbStatistics::Dirichlet::lnPdf(alphaForward, newVal);
        }
		catch (RbException &e)
        {
            lnProposalRatio = RbConstants::Double::neginf;
        }
    }
    
    variable->setValue( new Simplex(newVal), false );
    
    return lnProposalRatio;
}


/**
 *
 */
void DirichletSimplexProposal::prepareProposal( void )
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
void DirichletSimplexProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void DirichletSimplexProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( storedValue.clone() );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void DirichletSimplexProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Simplex>* >(newN) ;
    
}


void DirichletSimplexProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void DirichletSimplexProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate-p)/(1.0 - p)));
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }
    
    // set a hard maximum of 100
    alpha = fmin(100, alpha);
}


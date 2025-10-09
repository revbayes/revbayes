#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "DistributionNormal.h"
#include "NodeTimeSlidePathTruncatedNormalProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Proposal.h"
#include "RbConstants.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NodeTimeSlidePathTruncatedNormalProposal::NodeTimeSlidePathTruncatedNormalProposal( StochasticNode<Tree> *n, double s, bool a ) : Proposal(),
variable( n ),
sigma( s ),
scaleByAge( a )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NodeTimeSlidePathTruncatedNormalProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NodeTimeSlidePathTruncatedNormalProposal* NodeTimeSlidePathTruncatedNormalProposal::clone( void ) const
{
    
    return new NodeTimeSlidePathTruncatedNormalProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NodeTimeSlidePathTruncatedNormalProposal::getProposalName( void ) const
{
    static std::string name = "NodeTimeSlidePathTruncatedNormal";
    
    return name;
}


double NodeTimeSlidePathTruncatedNormalProposal::getProposalTuningParameter( void ) const
{
    return sigma;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NodeTimeSlidePathTruncatedNormalProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    if (tau.getNumberOfTips() <= 2)
    {
        return 0.0;
    }
    
    storedNodes.clear();
    storedAges.clear();
    TopologyNode* node;
    TopologyNode* tmpNode;
    double prop_sd;
    double current_age;
    double my_new_age;
    double UB;
    double LB;
    double lnHastingsRatio = 0.0;
    
    // For scaling proposal variances to match projections from BDP
    double ra = tau.getRoot().getAge();
    double r = log(tau.getNumberOfTips() / 2.0) / ra;
    
    // Move even or odd numbered nodes?
    size_t odds_evens = rng->uniform01() < 0.5 ? 0 : 1;
    
    // Pick tip node to start path
    size_t index = size_t( std::floor(tau.getNumberOfTips() * rng->uniform01()) );
    node = &tau.getTipNode(index);
    
    node = &node->getParent();
    
    size_t path_index = 1;
    bool is_even;
    // Trace path of all nodes between this tip and the root
    while ( !(node->isRoot()) )
    {
        is_even = path_index % 2 == 0 ? false : true;
        if ( (is_even && odds_evens == 0) || ((!is_even) && odds_evens == 1) )
        {
            
            current_age = node->getAge();
            
            tmpNode = &node->getParent();
            UB = tmpNode->getAge();
            
            LB = 0.0;
            for (size_t i=0; i < node->getNumberOfChildren(); ++i)
            {
                tmpNode = &node->getChild(i);
                if ( tmpNode->getAge() > LB )
                {
                    LB = tmpNode->getAge();
                }
            }
            
            prop_sd = scaleByAge == true ? sigma / 2 / exp(r * ((ra - (LB+UB)/2.0))) : sigma;
            
            my_new_age = RbStatistics::Normal::rv(current_age, prop_sd, LB, UB, *rng);
            
            tau.getNode(node->getIndex()).setAge( my_new_age );
            
            // q(x | x')
//            double new_prop_sd = scaleByAge == true ? sigma / 2 / exp(r * (ra - my_new_age)) : sigma;
            lnHastingsRatio += RbStatistics::Normal::lnPdf(my_new_age, prop_sd, current_age, LB, UB);
            // q(x' | x)
            lnHastingsRatio -= RbStatistics::Normal::lnPdf(current_age, prop_sd, my_new_age, LB, UB);
            
            storedNodes.push_back(node);
            storedAges.push_back(current_age);
            
            if (my_new_age == LB || my_new_age == UB)
            {
                // Can't have node ages equal, reject move
                return RbConstants::Double::neginf;
            }

        }
        node = &node->getParent();
        ++path_index;
    }
    
    return lnHastingsRatio;
    
}


/**
 *
 */
void NodeTimeSlidePathTruncatedNormalProposal::prepareProposal( void )
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
void NodeTimeSlidePathTruncatedNormalProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "sigma = ";
    if (name_only == false)
    {
        o << sigma;
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NodeTimeSlidePathTruncatedNormalProposal::undoProposal( void )
{
    
    // undo the proposal
    for (size_t i=0; i < storedNodes.size(); ++i)
    {
        storedNodes[i]->setAge( storedAges[i] );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NodeTimeSlidePathTruncatedNormalProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NodeTimeSlidePathTruncatedNormalProposal::setProposalTuningParameter(double tp)
{
    sigma = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void NodeTimeSlidePathTruncatedNormalProposal::tune( double rate )
{
    
    // Try tuning
    if ( rate > 0.44 )
    {
        sigma *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        sigma /= (2.0 - rate/0.44 );
    }
    
    // Get scale of tree
    double ra = variable->getValue().getRoot().getAge();
    
    // Avoid crazy operator values
    if ( sigma > ra ) {
        sigma = ra;
    } else if ( sigma < ra/10000 ) {
        sigma = sigma < ra/10000;
    }

    
}


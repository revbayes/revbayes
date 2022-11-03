#ifndef GibbsMixtureAllocationProposal_H
#define GibbsMixtureAllocationProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The allocation proposal between mixture elements.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    template <class mixtureType>
    class GibbsMixtureAllocationProposal : public Proposal {
        
    public:
        GibbsMixtureAllocationProposal( StochasticNode<mixtureType> *n, size_t d=0 );                                                                    //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        GibbsMixtureAllocationProposal*         clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<mixtureType>*            variable;                                                                           //!< The variable the Proposal is working on
        size_t                                  new_category;
        size_t                                  old_category;
        
    };
    
}


#include "MixtureDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "ReversibleJumpMixtureConstantDistribution.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template <class mixtureType>
RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::GibbsMixtureAllocationProposal( StochasticNode<mixtureType> *n, size_t d ) : Proposal(),
    variable( n ),
    new_category( 0 ),
    old_category( 0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template <class mixtureType>
void RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class mixtureType>
RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>* RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::clone( void ) const
{
    
    return new GibbsMixtureAllocationProposal<mixtureType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class mixtureType>
const std::string& RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::getProposalName( void ) const
{
    static std::string name = "Gibbs-Mixture-Allocation";
    
    return name;
}


/**
 * Perform the proposal.
 *
 * The reversible jump proposal switches the current "dimension".
 *
 * \return The hastings ratio.
 */
template <class mixtureType>
double RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::doProposal( void )
{
    
    // potential affected nodes for likelihood computation
    RbOrderedSet<DagNode *> affected;
    variable->initiateGetAffectedNodes( affected );
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    MixtureDistribution<mixtureType>& dist = static_cast<MixtureDistribution<mixtureType> &>( variable->getDistribution() );
    
    // get the number of categories
    size_t n = dist.getNumberOfMixtureElements();
    
    // create a vector for the weights
    std::vector<double> weights = std::vector<double>(n,0);
    double sum_of_weights = 0.0;
    double max_weight = RbConstants::Double::neginf;
    
    // get the current index
    old_category = dist.getCurrentIndex();
    
    for (size_t i=0; i<n; ++i)
    {
        // set our new value
        dist.setCurrentIndex( i );
        
        // flag for likelihood recomputation
        variable->touch();
        
        // compute the likelihood of the new value
        double prior_ratio = variable->getLnProbability();
        double likelihood_ratio = 0.0;
        for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
        {
            likelihood_ratio += (*it)->getLnProbability();
        }
        weights[i] = prior_ratio + likelihood_ratio;
        
        if (max_weight < weights[i])
        {
            max_weight = weights[i];
        }
        
    }
    
    // normalize weights
    for (size_t i=0; i<n; ++i)
    {
        weights[i] = exp(weights[i] - max_weight);
        sum_of_weights += weights[i];
    }
    
    double u = rng->uniform01() * sum_of_weights;
    new_category = 0;
    do
    {
        u -= weights[new_category];
        ++new_category;
    } while ( u > 0.0 );
    --new_category;
    
    double forwards  = weights[new_category];
    double backwards = weights[old_category];
    
    // set our new value
    dist.setCurrentIndex( new_category );
    
    return log( backwards / forwards );
}


/**
 *
 */
template <class mixtureType>
void RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::prepareProposal( void )
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
template <class mixtureType>
void RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    // nothing to print
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
template <class mixtureType>
void RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::undoProposal( void )
{
    
    MixtureDistribution<mixtureType>& dist = static_cast<MixtureDistribution<mixtureType> &>( variable->getDistribution() );
    dist.setCurrentIndex( old_category );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class mixtureType>
void RevBayesCore::GibbsMixtureAllocationProposal<mixtureType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<mixtureType>* >(newN) ;
    
}

#endif


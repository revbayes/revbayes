#ifndef UPPAllocationProposal_H
#define UPPAllocationProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The allocation proposal for the UPP.
     *
     * Will Freyman 6/30/17
     *
     */
    template <class mixtureType>
    class UPPAllocationProposal : public Proposal {
        
    public:
        UPPAllocationProposal( StochasticNode< RbVector<mixtureType> >* n );                                                                    //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        UPPAllocationProposal*              clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:

        // parameters
        StochasticNode< RbVector<mixtureType> >*  variable;                                                                           //!< The variable the Proposal is working on
        size_t                                    old_partition;
        std::vector<int>                          old_value_assignments;
        
    };
    
}


#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "ReversibleJumpMixtureConstantDistribution.h"
#include "TypedDagNode.h"
#include "UniformPartitioningDistribution.h"

#include <cmath>
#include <iostream>


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template <class mixtureType>
RevBayesCore::UPPAllocationProposal<mixtureType>::UPPAllocationProposal( StochasticNode< RbVector<mixtureType> >* n ) : Proposal(),
    variable( n ),
    old_partition( 0 )
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
void RevBayesCore::UPPAllocationProposal<mixtureType>::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class mixtureType>
RevBayesCore::UPPAllocationProposal<mixtureType>* RevBayesCore::UPPAllocationProposal<mixtureType>::clone( void ) const
{
    
    return new UPPAllocationProposal<mixtureType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class mixtureType>
const std::string& RevBayesCore::UPPAllocationProposal<mixtureType>::getProposalName( void ) const
{
    static std::string name = "UPP-Allocation";
    
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
double RevBayesCore::UPPAllocationProposal<mixtureType>::doProposal( void )
{
    
    UniformPartitioningDistribution<mixtureType>& dist = static_cast<UniformPartitioningDistribution<mixtureType> &>( variable->getDistribution() );
    
    // get the current index
    old_partition = dist.getCurrentIndex();
    
    // get the value assignments
    old_value_assignments = dist.getValueAssignments();

    // draw a new random partition
    dist.redrawValue();
    
    return 0.0;
    
}


/**
 *
 */
template <class mixtureType>
void RevBayesCore::UPPAllocationProposal<mixtureType>::prepareProposal( void )
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
void RevBayesCore::UPPAllocationProposal<mixtureType>::printParameterSummary(std::ostream &o, bool name_only) const
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
void RevBayesCore::UPPAllocationProposal<mixtureType>::undoProposal( void )
{
    
    UniformPartitioningDistribution<mixtureType>& dist = static_cast<UniformPartitioningDistribution<mixtureType> &>( variable->getDistribution() );
    dist.setCurrentIndex( (int)old_partition );
    dist.setValueAssignments( old_value_assignments );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class mixtureType>
void RevBayesCore::UPPAllocationProposal<mixtureType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode< RbVector<mixtureType> >* >(newN) ;
    
}

#endif


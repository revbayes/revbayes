#ifndef MixtureAllocationProposal_H
#define MixtureAllocationProposal_H

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
    class MixtureAllocationProposal : public Proposal {
        
    public:
        MixtureAllocationProposal( StochasticNode<mixtureType> *n, size_t d=0 );                                                                    //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        MixtureAllocationProposal*          clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                              getProposalTuningParameter(void) const;
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                setProposalTuningParameter(double tp);
        void                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:

        // parameters
        StochasticNode<mixtureType>*        variable;                                                                           //!< The variable the Proposal is working on
        size_t                              delta;
        size_t                              newCategory;
        size_t                              oldCategory;
        
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
RevBayesCore::MixtureAllocationProposal<mixtureType>::MixtureAllocationProposal( StochasticNode<mixtureType> *n, size_t d ) : Proposal(),
    variable( n ),
    delta( d ),
    newCategory( 0 ),
    oldCategory( 0 )
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
void RevBayesCore::MixtureAllocationProposal<mixtureType>::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class mixtureType>
RevBayesCore::MixtureAllocationProposal<mixtureType>* RevBayesCore::MixtureAllocationProposal<mixtureType>::clone( void ) const
{
    
    return new MixtureAllocationProposal<mixtureType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class mixtureType>
const std::string& RevBayesCore::MixtureAllocationProposal<mixtureType>::getProposalName( void ) const
{
    static std::string name = "Mixture-Allocation";
    
    return name;
}


template <class mixtureType>
double RevBayesCore::MixtureAllocationProposal<mixtureType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * The reversible jump proposal switches the current "dimension".
 *
 * \return The hastings ratio.
 */
template <class mixtureType>
double RevBayesCore::MixtureAllocationProposal<mixtureType>::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    MixtureDistribution<mixtureType>& dist = static_cast<MixtureDistribution<mixtureType> &>( variable->getDistribution() );
    
    // get the number of categories
    size_t n = dist.getNumberOfMixtureElements();
    
    // get the current index
    oldCategory = dist.getCurrentIndex();
    
    // draw a new category uniformly at random
    if (delta == 0)
    {
        newCategory = (int)(rng->uniform01() * n);
    }
    // draw a new category from +/- delta (for ordered states)
    else
    {
     
        // new step size drawn (use hypergeometric in future?)
        int sign = ( rng->uniform01() < 0.5 ? -1 : 1 );
        int d = sign * ( (int)(rng->uniform01() * (delta-1)) + 1 );
        
        // tmp value used to support negative values
        int tmp = (int)oldCategory + d;
        
        /* reflect the new value */
        int min = 0;
        int max = (int)(n-1);
        do {
            if ( tmp < min )
            {
                tmp = 2 * min - tmp;
            }
            else if ( tmp > max )
            {
                tmp = 2 * max - tmp;
            }
        } while ( tmp < min || tmp > max );
        
        newCategory = tmp;
    }
    
    // set our new value
    dist.setCurrentIndex( newCategory );
    
    return 0.0;
    
}


/**
 *
 */
template <class mixtureType>
void RevBayesCore::MixtureAllocationProposal<mixtureType>::prepareProposal( void )
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
void RevBayesCore::MixtureAllocationProposal<mixtureType>::printParameterSummary(std::ostream &o, bool name_only) const
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
void RevBayesCore::MixtureAllocationProposal<mixtureType>::undoProposal( void )
{
    
    MixtureDistribution<mixtureType>& dist = static_cast<MixtureDistribution<mixtureType> &>( variable->getDistribution() );
    dist.setCurrentIndex( oldCategory );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class mixtureType>
void RevBayesCore::MixtureAllocationProposal<mixtureType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<mixtureType>* >(newN) ;
    
}


template <class mixtureType>
void RevBayesCore::MixtureAllocationProposal<mixtureType>::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
template <class mixtureType>
void RevBayesCore::MixtureAllocationProposal<mixtureType>::tune( double rate )
{
    // nothing to do here.
    
}



#endif


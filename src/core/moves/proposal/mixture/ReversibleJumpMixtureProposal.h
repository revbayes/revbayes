#ifndef ReversibleJumpMixtureProposal_H
#define ReversibleJumpMixtureProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The reversible jump proposal to switch between two elements of an RJ-Mixture.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    template <class mixtureType>
    class ReversibleJumpMixtureProposal : public Proposal {
        
    public:
        ReversibleJumpMixtureProposal( StochasticNode<mixtureType> *n);                                                                    //!<  constructor
        ReversibleJumpMixtureProposal( const ReversibleJumpMixtureProposal &p );
        virtual ~ReversibleJumpMixtureProposal();
        
        ReversibleJumpMixtureProposal&      operator=(const ReversibleJumpMixtureProposal& p);
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        ReversibleJumpMixtureProposal*      clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode<mixtureType>*        variable;                                                                           //!< The variable the Proposal is working on
        mixtureType*                        stored_value;                                                                        //!< The stored value of the Proposal used for rejections.
        size_t                              stored_index;

    };
    
}


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
RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::ReversibleJumpMixtureProposal( StochasticNode<mixtureType> *n ) : Proposal(),
    variable( n ),
    stored_value( NULL ),
    stored_index( 0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::ReversibleJumpMixtureProposal( const ReversibleJumpMixtureProposal &p ) : Proposal(),
    variable( p.variable ),
    stored_value( NULL ),
    stored_index( p.stored_index )
{
    // tell the base class to add the node
    addNode( variable );
    
    if ( p.stored_value != NULL )
    {
        Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( *p.stored_value );
    }
    
}


/**
 * Destructor
 *
 */
template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::~ReversibleJumpMixtureProposal( )
{
    delete stored_value;
    
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::cleanProposal( void )
{
    
    delete stored_value;
    stored_value = NULL;
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>& RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::operator=( const ReversibleJumpMixtureProposal &p )
{
    
    if ( this != &p )
    {
        Proposal::Cloneable::operator=( p );
        
        delete stored_value;
        
        stored_value    = NULL;
        variable        = p.variable;
        stored_index    = p.stored_index;

        if ( p.stored_value != NULL )
        {
            Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( *p.stored_value );
        }
        
    }
    
    return *this;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class mixtureType>
RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>* RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::clone( void ) const
{
    
    return new ReversibleJumpMixtureProposal<mixtureType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class mixtureType>
const std::string& RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::getProposalName( void ) const
{
    static std::string name = "RJ-Switch";
    
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
double RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::doProposal( void )
{
    
    
    mixtureType& v = variable->getValue();
    ReversibleJumpMixtureConstantDistribution<mixtureType> &d = static_cast< ReversibleJumpMixtureConstantDistribution<mixtureType>& >( variable->getDistribution() );
    
    // copy value
    stored_value = Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( v );
    stored_index = d.getCurrentIndex();
    
    double lnHastingsratio = 0.0;
    
    if ( stored_index == 0 )
    {
        // draw the new value
        d.redrawValueByIndex( 1 );

        // get the base distribution
        TypedDistribution<mixtureType> &baseDistribution = d.getBaseDistribution();
        
        // store the proposal ratio
        lnHastingsratio = - baseDistribution.computeLnProbability();
    }
    else
    {
        
        // get the base distribution
        TypedDistribution<mixtureType> &baseDistribution = d.getBaseDistribution();
    
        // store the proposal ratio
        lnHastingsratio = baseDistribution.computeLnProbability();
        
        // draw the new value
        d.redrawValueByIndex( 0 );
    }
    
    return lnHastingsratio;
}


/**
 *
 */
template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::prepareProposal( void )
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
void RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::printParameterSummary(std::ostream &o, bool name_only) const
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
void RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::undoProposal( void )
{
    
    // swap current value and stored value
    variable->setValue( Cloner<mixtureType, IsDerivedFrom<mixtureType, Cloneable>::Is >::createClone( *stored_value ) );
    
    // also reset the index
    ReversibleJumpMixtureConstantDistribution<mixtureType> &d = static_cast< ReversibleJumpMixtureConstantDistribution<mixtureType>& >( variable->getDistribution() );
    d.setCurrentIndex( stored_index );
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class mixtureType>
void RevBayesCore::ReversibleJumpMixtureProposal<mixtureType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<mixtureType>* >(newN) ;
    
}


#endif


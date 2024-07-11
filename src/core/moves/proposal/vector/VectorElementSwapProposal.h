#ifndef VectorElementSwapProposal_H
#define VectorElementSwapProposal_H

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
    template <class valueType>
    class VectorElementSwapProposal : public Proposal {
        
    public:
        VectorElementSwapProposal( const std::vector<StochasticNode<valueType> *>& n, bool neighors_only );                                                   //!< Constructor

        // Basic utility functions
        void                                        cleanProposal(void);                                                            //!< Clean up proposal
        VectorElementSwapProposal*                  clone(void) const;                                                              //!< Clone object
        double                                      doProposal(void);                                                               //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                    //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        prepareProposal(void);                                                          //!< Prepare the proposal
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                             //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                 //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:

        // parameters
        std::vector<StochasticNode<valueType> *>    variables;                                                                      //!< The vector to operate on

        bool                                        neighboring;                                                                    //!< Should we only swap neighbors?
        size_t                                      length;                                                                         //!< The index of the last modified element.
        size_t                                      stored_index_1;                                                                 //!< The stored value of the last modified element.
        size_t                                      stored_index_2;                                                                 //!< The stored value of the last modified element.
                
    };
    
}


#include "Cloner.h"
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
template <class valueType>
RevBayesCore::VectorElementSwapProposal<valueType>::VectorElementSwapProposal( const std::vector<StochasticNode<valueType> *>& n, bool neighors_only ) : Proposal(),
    variables( n ),
    neighboring( neighors_only ),
    length( variables.size() ),
    stored_index_1( 0 ),
    stored_index_2( 0 )
{
    
    // tell the base class to add the node
    typename std::vector< StochasticNode<valueType> *>::const_iterator it;
    for ( it = variables.begin(); it != variables.end(); it++)
    {
        addNode( *it );
    }
    
    // set the length variables for internal use
    length = variables.size();
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template <class valueType>
RevBayesCore::VectorElementSwapProposal<valueType>* RevBayesCore::VectorElementSwapProposal<valueType>::clone( void ) const
{
    
    return new VectorElementSwapProposal<valueType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template <class valueType>
const std::string& RevBayesCore::VectorElementSwapProposal<valueType>::getProposalName( void ) const
{
    static std::string name = "VectorElementSwap";
    
    return name;
}


template <class valueType>
double RevBayesCore::VectorElementSwapProposal<valueType>::getProposalTuningParameter( void ) const
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
template <class valueType>
double RevBayesCore::VectorElementSwapProposal<valueType>::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // either randomly draw to indices or pick to neighbors
    if ( neighboring == true )
    {
        // Randomly draw one index and set the second one to its neighbor
        stored_index_1 = size_t( floor(rng->uniform01() * double(length-1)) );
        stored_index_2 = stored_index_1 + 1;
    }
    else
    {
        // Randomly draw two indices (not necessarily unique)
        stored_index_1 = size_t( floor(rng->uniform01() * double(length)) );
        stored_index_2 = size_t( floor(rng->uniform01() * double(length)) );
    }
    
    // Swap the values located at the chosen indices
    valueType* tmp   = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( variables[stored_index_1]->getValue() );
    variables[stored_index_1]->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone(variables[stored_index_2]->getValue()) );
    variables[stored_index_2]->setValue( tmp );
    
    // This move is symmetric, so lnHastings = 0
    double ln_Hastings_ratio = 0.0;

    return ln_Hastings_ratio;
    
}


/**
 *
 */
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::prepareProposal( void )
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
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::printParameterSummary(std::ostream &o, bool name_only) const
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
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::undoProposal( void )
{
    
    valueType* tmp   = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( variables[stored_index_1]->getValue() );
    variables[stored_index_1]->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( variables[stored_index_2]->getValue() ) );
    variables[stored_index_2]->setValue( tmp );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    for (size_t i = 0; i < variables.size(); ++i)
    {
        if ( variables[i] == oldN )
        {
            variables[i] = static_cast<StochasticNode<valueType> *>(newN);
        }
    }
    
}


template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::setProposalTuningParameter(double tp)
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
template <class valueType>
void RevBayesCore::VectorElementSwapProposal<valueType>::tune( double rate )
{
    // nothing to do here.
    
}



#endif

#ifndef IndependentPriorProposal_H
#define IndependentPriorProposal_H

#include <string>

#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * An independent prior distribution proposal.
     *
     * Propose a new value from the prior distribution.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Pett)
     * @since 2012-07-12, version 1.0
     *
     */
    template <class valueType>
    class IndependentPriorProposal : public Proposal {
        
    public:
        IndependentPriorProposal( StochasticNode<valueType> *n, bool mh = true);                            //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        IndependentPriorProposal*               clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        StochasticNode<valueType>*              variable;
        valueType                               stored_value;

    };
    
}


template <class valueType>
RevBayesCore::IndependentPriorProposal<valueType>::IndependentPriorProposal( RevBayesCore::StochasticNode<valueType> *v, bool mh) : RevBayesCore::Proposal(),
    variable( v ),
    stored_value(v->getValue())
{

    addNode( variable );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class valueType>
RevBayesCore::IndependentPriorProposal<valueType>* RevBayesCore::IndependentPriorProposal<valueType>::clone( void ) const
{

    return new IndependentPriorProposal( *this );
}


/** Clean the proposal */
template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::cleanProposal( void )
{

    // nothing to do
}


template <class valueType>
const std::string& RevBayesCore::IndependentPriorProposal<valueType>::getProposalName( void ) const
{
    static std::string name = "IndependentPriorMove";

    return name;
}


template <class valueType>
double RevBayesCore::IndependentPriorProposal<valueType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/** Perform the move */
template <class valueType>
double RevBayesCore::IndependentPriorProposal<valueType>::doProposal( void )
{

    stored_value = variable->getValue();
    
    double ln_back_prob = variable->getLnProbability();

    variable->getDistribution().redrawValue();

    variable->touch( true );
    
    double ln_forw_prob = variable->getLnProbability();

    return ln_back_prob - ln_forw_prob;
}


/** Prepare the proposal */
template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::prepareProposal( void )
{

    // nothing to do
}



/** Print the parameter summary */
template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::printParameterSummary(std::ostream &o, bool name_only) const
{

    // nothing to do
}


template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::undoProposal( void )
{

    // reset to the old tree
    variable->setValue(new valueType(stored_value));

}

template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/** Tune the proposal */
template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::tune(double r)
{

    // nothing to do
}


template <class valueType>
void RevBayesCore::IndependentPriorProposal<valueType>::swapNodeInternal(RevBayesCore::DagNode *oldN, RevBayesCore::DagNode *newN)
{

    variable = static_cast<StochasticNode<valueType>* >( newN );

}


#endif

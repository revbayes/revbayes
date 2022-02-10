#ifndef SimpleProposal_H
#define SimpleProposal_H

#include "Proposal.h"
#include "StochasticNode.h"

#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    class Move;
    
    /**
     * Base class for all proposal distribution.
     *
     * The base class of all moves only provides the interface for the call to propose a new move.
     * Here the perform methods actually does the accept/reject step.
     * All specifics are implemented in the derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-27, version 1.0
     *
     */
    template <class valueType>
    class SimpleProposal : public Proposal {
        
    public:
        virtual                                                ~SimpleProposal(void);                                                                   //!< Destructor
        
        // public methods
        virtual double                                          doProposal(void);                                                                       //!< Actually do the proposal.
        
        // pure virtual public methods
        virtual void                                            cleanProposal(void) = 0;                                                                //!< Cleanup proposal
        virtual SimpleProposal<valueType>*                      clone(void) const = 0;                                                                  //!< Make a deep copy
        virtual const std::string&                              getProposalName(void) const = 0;                                                        //!< Get the name of this proposal used for printing out info.
        virtual double                                          getProposalTuningParameter(void) const = 0;
        virtual void                                            prepareProposal(void) = 0;                                                              //!< Propose a new state
        virtual void                                            printParameterSummary(std::ostream &o, bool name_only) const = 0;                       //!< Print the parameter summary
        virtual void                                            setProposalTuningParameter(double tp) = 0;
        virtual void                                            tune(double r) = 0;                                                                     //!< Tune the parameters of the proposal.
        virtual void                                            undoProposal(void) = 0;                                                                 //!< Undo the proposal.
        
        virtual double                                          propose(valueType &v) = 0;                                                              //!< Actually do the proposal.
        
        
    protected:
        SimpleProposal(StochasticNode<valueType> *v, double p=0.44);                                                                                                                           //!< Default constructor
        SimpleProposal(const SimpleProposal<valueType> &p);                                                                                             //!< Copy constructor
        SimpleProposal<valueType>&                              operator=(const SimpleProposal<valueType> &p);                                          //!< Assignment operator
        
        
        
        // pure virtual protected method that need to be implemented in derived classes
        virtual void                                            swapNodeInternal(DagNode *oldN, DagNode *newN) = 0;                                     //!< Swap the pointers to the variable on which the move works on.
        
        StochasticNode<valueType>*                              variable;                                                                               //!< The variable the Proposal is working on

    private:
        
        
    };
    
}


template <class valueType>
RevBayesCore::SimpleProposal<valueType>::SimpleProposal( StochasticNode<valueType> *v, double p ) : Proposal( p ),
    variable( v )
{
    
    // tell the base class to add the node
    addNode( variable );
}


template <class valueType>
RevBayesCore::SimpleProposal<valueType>::SimpleProposal(const SimpleProposal<valueType> &p) : Proposal( p ),
    variable( p.variable )
{
    
}


template <class valueType>
RevBayesCore::SimpleProposal<valueType>::~SimpleProposal( void )
{
    
}



template <class valueType>
RevBayesCore::SimpleProposal<valueType>& RevBayesCore::SimpleProposal<valueType>::operator=(const SimpleProposal<valueType> &p)
{
    
    // just delegate to the base class
    Proposal::operator=( p );
    
    if ( this != &p )
    {
        variable = p.variable;
    }
    
    return *this;
}


template <class valueType>
double RevBayesCore::SimpleProposal<valueType>::doProposal( void )
{
    
    valueType &val = variable->getValue();
    
    // copy value
    double ln_hastings_ratio = propose( val );

    return ln_hastings_ratio;
}




#endif

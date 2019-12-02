#ifndef TreeProposal_H
#define TreeProposal_H

#include "Proposal.h"

#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    class TopologyNode;
    
    /**
     * Base class for all tree proposals. 
     *
     * The base class of all tree proposal only provides and interface so that moves and other proposals
     * can be used joint with any tree proposal.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-04-08, version 1.0
     *
     */
    class TreeProposal : public Proposal {
        
    public:
        virtual                                                ~Proposal(void) {}                                                                       //!< Destructor
        
        // pure virtual public methods
        virtual Proposal*                                       clone(void) const = 0;                                                                  //!< Make a deep copy
        virtual double                                          doProposal(void) = 0;                                                                   //!< Actually do the proposal.
        virtual std::vector<TopologyNode*>                      getChosenTopologyNodes(void) const = 0;                                                 //!< Get the vector of chosen topology nodes that will be pruned.
        virtual const std::vector<DagNode*>&                    getNodes(void) const = 0;                                                               //!< Get the vector of nodes for which the proposal is drawing new values.
        virtual const std::string&                              getProposalName(void) const = 0;                                                        //!< Get the name of this proposal used for printing out info.
        virtual double                                          getProposalTuningParameter(void) const = 0;
        virtual void                                            prepareProposal(void) = 0;                                                              //!< Propose a new state
        virtual void                                            printParameterSummary(std::ostream &o, bool name_only) const = 0;                                       //!< Print the parameter summary
        virtual void                                            setProposalTuningParameter(double tp) = 0;
        virtual void                                            swapNode(DagNode *oldN, DagNode *newN) = 0;                                             //!< Swap the pointers to the variable on which the move works on.
        virtual void                                            tune(double r) = 0;                                                                     //!< Tune the parameters of the proposal.
        virtual void                                            undoProposal(void) = 0;                                                                 //!< Undo the proposal.
        
    protected:
        Proposal() {}                                                                                                                                   //!< Constructor
        
        
    };
    
}


#endif

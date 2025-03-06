#ifndef Proposal_H
#define Proposal_H

#include "Cloneable.h"

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
    class Proposal : public Cloneable {
        
    public:
        virtual                                                ~Proposal(void);                                                                         //!< Destructor
        
        // public methods
        virtual bool                                            allowClamped() const { return false; }                                                  //!< Allow proposal to operate on clamped nodes; see #600
        const std::vector<DagNode*>&                            getNodes(void) const;                                                                   //!< Get the vector of nodes for which the proposal is drawing new values.
        void                                                    swapNode(DagNode *oldN, DagNode *newN);                                                 //!< Swap the pointers to the variable on which the move works on.
        void                                                    setMove(Move *m);                                                                       //!< Set the pointer to move object holding this proposal
        const Move*                                             getMove(void) const;                                                                    //!< Get the pointer to move object holding this proposal

        // pure virtual public methods
        virtual void                                            cleanProposal(void) = 0;                                                                //!< Cleanup proposal
        virtual Proposal*                                       clone(void) const = 0;                                                                  //!< Make a deep copy
        virtual double                                          doProposal(void) = 0;                                                                   //!< Actually do the proposal.
        virtual const std::string&                              getProposalName(void) const = 0;                                                        //!< Get the name of this proposal used for printing out info.
        virtual std::string                                     getLongProposalName(void) const;                                                        //!< Get the name of this proposal used for printing out info.
        virtual double                                          getProposalTuningParameter(void) const = 0;
        virtual std::vector<DagNode*>                           identifyNodesToTouch(void);
        virtual void                                            prepareProposal(void) = 0;                                                              //!< Propose a new state
        virtual void                                            printParameterSummary(std::ostream &o, bool name_only) const = 0;                       //!< Print the parameter summary
        virtual void                                            setProposalTuningParameter(double tp) = 0;
        virtual void                                            tune(double r) = 0;                                                                     //!< Tune the parameters of the proposal.
        virtual void                                            undoProposal(void) = 0;                                                                 //!< Undo the proposal.
        
        
    protected:
        Proposal(double p=0.44);                                                                                                                                 //!< Default constructor
        Proposal(const Proposal &p);                                                                                                                    //!< Copy constructor
        Proposal&                                               operator=(const Proposal &p);                                                           //!< Assignment operator
        
        
        // helper methods for this and derived classes
        void                                                    addNode(DagNode* p);                                                                    //!< add a node to the proposal
        void                                                    removeNode(DagNode* p);                                                                 //!< remove a node from the proposal

        
        // pure virtual protected method that need to be implemented in derived classes
        virtual void                                            swapNodeInternal(DagNode *oldN, DagNode *newN) = 0;                                     //!< Swap the pointers to the variable on which the move works on.
        
        // members
        std::vector<DagNode*>                                   nodes;
        Move*                                                   move;
        double                                                  targetAcceptanceRate;

        
    };
    
}


#endif

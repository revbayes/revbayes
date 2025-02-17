#ifndef EventBirthDeathProposal_H
#define EventBirthDeathProposal_H

#include <set>
#include <string>

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "CharacterEvent.h"

namespace RevBayesCore {
    
    /**
     * The birth-death proposal for events along a tree.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since 2009-09-08, version 1.0
     *
     */
    class EventBirthDeathProposal : public Proposal {
        
    public:
        virtual EventBirthDeathProposal*        clone(void) const override = 0;                                                     //!< Clone object
        virtual const std::string&              getProposalName(void) const override = 0;                                           //!< Get the name of the proposal for summary printing
        
        
        // Basic utility functions
        bool                                    allowClamped() const override { return true; }                                      //!< Proposal doesn't change the tree, but changes parameters describing the process that generates the tree. See #600
        void                                    cleanProposal(void) override;                                                       //!< Clean up proposal
        double                                  doProposal(void) override;                                                          //!< Perform proposal
        virtual void                            initialize();                                                                       //!< Initialize the proposal
        double                                  getProposalTuningParameter(void) const override;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const override;              //!< Print the parameter summary
        void                                    prepareProposal(void) override;                                                     //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp) override;
        void                                    tune(double r) override;                                                            //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void) override;                                                        //!< Reject the proposal
        
    protected:
        EventBirthDeathProposal( StochasticNode<Tree> *n);                                                                          //!< Constructor

        // pure virtual methods
        virtual CharacterEvent*                 drawNewEvent(double event_time) = 0;
        virtual double                          computeEventProposalProbability( CharacterEvent* event ) = 0;

        double                                  doBirthProposal(void);
        double                                  doDeathProposal(void);
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN) override;                            //!< Swap the DAG nodes on which the Proposal is working on
        
        // parameters
        AbstractCharacterHistoryBirthDeathProcess* distribution;

    private:
        
        // parameters
        StochasticNode<Tree>*                   variable;                                                                           //!< The variable the Proposal is working on
        
        size_t accepted_birth;
        size_t trie_birth;
        size_t accepted_death;
        size_t trie_death;

        CharacterEvent*                         stored_value;                                                                       //!< The stored value of the Proposal used for rejections.
        size_t                                  stored_branch_index;
        bool                                    was_birth_proposal;                                                                 //!< The value we propose.
    };
    
}

#endif


#ifndef ContinuousEventBirthDeathProposal_H
#define ContinuousEventBirthDeathProposal_H

#include <set>
#include <string>

#include "EventBirthDeathProposal.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    class CharacterEvent;
    
    /**
     * The birth-death proposal for events along a tree.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since 2009-09-08, version 1.0
     *
     */
    class ContinuousEventBirthDeathProposal : public EventBirthDeathProposal {
        
    public:
        ContinuousEventBirthDeathProposal( StochasticNode<Tree> *n);                                                                //!<  constructor
        
        // Basic utility functions
        bool                                    allowClamped() const override { return true; }                                      //!< Proposal doesn't change the tree, but changes parameters describing the process that generates the tree. See #600
        ContinuousEventBirthDeathProposal*      clone(void) const;                                                                  //!< Clone object
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        
    protected:
        
        // pure virtual methods
        CharacterEvent*                         drawNewEvent(double event_time);
        double                                  computeEventProposalProbability( CharacterEvent* event );
        
    };
    
}

#endif



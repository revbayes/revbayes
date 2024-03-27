#ifndef VectorElementsSwapProposal_H
#define VectorElementsSwapProposal_H

#include <iosfwd>
#include <set>

#include "Proposal.h"
#include "Simplex.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The element swap simplex operator.
     *
     * An element swap simplex proposal simply selects two random elements of a simplex
     * and swaps them.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class VectorElementsSwapProposal : public Proposal {
        
    public:
        VectorElementsSwapProposal( std::vector<StochasticNode<double> *> n );                                                //!<  constructor
        
        // Basic utility functions
        void                                    addIndex(size_t index);
        
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        VectorElementsSwapProposal*             clone(void) const;                                          //!< Clone object
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
        
        // parameters
        std::vector<StochasticNode<double> *>   variables;
        std::set<size_t>                        indices;
        
        // stored objects to undo proposal
        size_t                                  stored_chosenidx1;
        size_t                                  stored_chosenidx2;
        double                                  stored_variable1;
        double                                  stored_variable2;
    };
    
}

#endif

#ifndef RandomCategoryWalkProposal_H
#define RandomCategoryWalkProposal_H

#include <ostream>
#include <vector>
#include <cstddef>
#include <set>

#include "Proposal.h"

#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Random category walk proposal.
     *
     *
     * This proposal is intended for a multinomial distribution. It takes one element of the vector,
     * decreases it by one (if possible) and either adds to the left or right neighboring category/element 1.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     *
     */
    class RandomCategoryWalkProposal : public Proposal {
        
    public:
        RandomCategoryWalkProposal( StochasticNode< RbVector<long> >* n);                                 //!< Constructor
        
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        RandomCategoryWalkProposal*                 clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        StochasticNode< RbVector<long> >*           variable;
        size_t                                      chosen_index;
        size_t                                      chosen_neighbour;
        bool                                        failed;
        
        
    };
    
}

#endif


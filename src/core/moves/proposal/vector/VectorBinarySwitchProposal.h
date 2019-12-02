#ifndef VectorBinarySwitchProposal_H
#define VectorBinarySwitchProposal_H

#include "RbVector.h"
#include "Proposal.h"
#include "StochasticNode.h"

#include <ostream>
#include <vector>
#include <string>

namespace RevBayesCore {
    
    /**
     * @brief Sliding proposal of a all elements of a vector.
     *
     *
     * This proposal randomly flips all binary elements of a vector.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     *
     */
    class VectorBinarySwitchProposal : public Proposal {
        
    public:
        VectorBinarySwitchProposal( std::vector< StochasticNode<long>* > n, const std::vector<long> &i, double sp);                                 //!< Constructor
        
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        VectorBinarySwitchProposal*                 clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        std::vector<DagNode*>                       identifyNodesToTouch(void);
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        std::vector<StochasticNode<long>* >         variables;
        std::vector<long>                           indices;
        size_t                                      length;
        double                                      switch_probability;
        std::set<size_t>                            update_set;
        
        
    };
    
}

#endif


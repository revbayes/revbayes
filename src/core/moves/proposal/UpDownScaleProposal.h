#ifndef UpDownScaleProposal_H
#define UpDownScaleProposal_H

#include "RbVector.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

#include <ostream>
#include <vector>
#include <string>

namespace RevBayesCore {
    
    /**
     * @brief Up-Down Scaling Proposal for several parameters jointly.
     *
     *
     * This proposal randomly scales all a set of parameter up while the other set of parameters
     * is scaled down by the same value. This should hopefully improve mixing in many cases.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-25, version 1.0
     *
     */
    class UpDownScaleProposal : public Proposal {
        
    public:
        UpDownScaleProposal(double l);                                 //!< Constructor
        
        void                                        addVariable(StochasticNode<double> *v, bool up);                                    //!< Add an up-scaling variable
        void                                        addVariable(StochasticNode<RbVector<double> > *v, bool up);                         //!< Add an up-scaling variable
        void                                        addVariable(StochasticNode<Tree> *v, bool up);                                      //!< Add an up-scaling variable
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        UpDownScaleProposal*                        clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        removeVariable(StochasticNode<double> *v, bool up);                                 //!< Add an up-scaling variable
        void                                        removeVariable(StochasticNode<RbVector<double> > *v, bool up);                      //!< Add an up-scaling variable
        void                                        removeVariable(StochasticNode<Tree> *v, bool up);                                   //!< Add an up-scaling variable
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        std::vector<StochasticNode<double> *>               up_scalar_variables;
        std::vector<StochasticNode<RbVector<double> > *>    up_vector_variables;
        std::vector<StochasticNode<Tree> *>                 up_tree_variables;
        std::vector<StochasticNode<double> *>               down_scalar_variables;
        std::vector<StochasticNode<RbVector<double> > *>    down_vector_variables;
        std::vector<StochasticNode<Tree> *>                 down_tree_variables;
        
        double                                              lambda;                                                                     //!< The scale parameter of the Proposal (larger lambda -> larger proposals).
        std::vector<double>                                 stored_up_scalar_values;
        std::vector<std::vector<double> >                   stored_up_vector_values;
        std::vector<std::vector<double> >                   stored_up_tree_values;
        std::vector<double>                                 stored_down_scalar_values;
        std::vector<std::vector<double> >                   stored_down_vector_values;
        std::vector<std::vector<double> >                   stored_down_tree_values;
        
    };
    
}

#endif


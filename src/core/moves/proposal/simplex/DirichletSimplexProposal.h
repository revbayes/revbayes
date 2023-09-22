#ifndef DirichletSimplexProposal_H
#define DirichletSimplexProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"
#include "Simplex.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The Dirichlet-simplex operator.
     *
     * A Dirichlet-simplex proposal randomly changes some values of a simplex, although the other values
     * change too because of the renormalization.
     * First, some random indices are drawn. Then, the proposal draws a new somplex 
     *   u ~ Dirichlet(val[index] * alpha)
     * where alpha is the tuning parameter.The new value is set to u.
     * The simplex is then renormalized.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class DirichletSimplexProposal : public Proposal {
        
    public:
        DirichletSimplexProposal( StochasticNode<Simplex> *n, double a, size_t nc, double o, double k=0.0, double p=0.234);                                                                    //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        DirichletSimplexProposal*               clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode<Simplex>*                variable;                                                                           //!< The variable the Proposal is working on
        Simplex                                 storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        double                                  alpha;                                                                             //!< The scaling parameter of the Proposal
        size_t                                  nCategories;
        double                                  offset;
        double									kappa;
//        double                                  proposedValue;                                                                      //!< The value we propose.
    };
    
}

#endif


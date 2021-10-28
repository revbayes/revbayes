#ifndef VectorSimplexSwapMove_H
#define VectorSimplexSwapMove_H

#include <stddef.h>
#include <ostream>

#include "Proposal.h"
#include "Simplex.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
//template <class variableType> class StochasticNode;
template <class variableType> class DeterministicNode;
    
    /**
     * @brief Swap 2 elements randomly from a vector of simplexes.
     *
     * This move randomly picks 2 indexes of a vector of simplexes, 
     * and swaps the the simplexes at these indices.
     */
    class VectorSimplexSwapProposal : public Proposal {
        
    public:
        //VectorSimplexSwapProposal( StochasticNode<RbVector<Simplex> >* n );                                 //!< Constructor
        VectorSimplexSwapProposal( DeterministicNode<RbVector<Simplex> >* n );                                 //!< Constructor
        
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        VectorSimplexSwapProposal*                  clone(void) const;                                                                  //!< Clone object
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
        
        //StochasticNode<RbVector<Simplex> >*         variable;
        DeterministicNode<RbVector<Simplex> >*         variable;
        
        size_t                                      length;                                                                              //!< The index of the last modified element.
        size_t                                      storedIndex1;                                                                        //!< The stored value of the last modified element.
        size_t                                      storedIndex2;                                                                        //!< The stored value of the last modified element.
        
    };
    
}

#endif


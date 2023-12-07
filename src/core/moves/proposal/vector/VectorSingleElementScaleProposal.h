#ifndef VectorSingleElementScaleMove_H
#define VectorSingleElementScaleMove_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Scaling proposal of a single element randomly picked from a Stochastic vector.
     *
     * This proposal randomly picks an element of a vector of positive real numbers,
     * proposes a scaling factor and then scales the value.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     * It generally makes more sense to apply the scaling proposal on a vector of positive
     * real numbers but technically it works on negative numbers too. However, 
     * the proposal will never change the sign of the value and thus is incomplete if applied
     * to a variable defined on the whole real line.
     *
     * @see SingleElementScaleProposal for the same applied to a vector of Stochastic nodes.
     * @see RevLanguage::Move_VectorSingleElementScale for the RL interface.
     *
     */
    class VectorSingleElementScaleProposal : public Proposal {
        
    public:
        VectorSingleElementScaleProposal(StochasticNode<RbVector<double> >* n, double l);
        
        void                                        cleanProposal(void);
        VectorSingleElementScaleProposal*           clone(void) const;
        double                                      doProposal(void);
        const std::string&                          getProposalName(void) const;
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;
        void                                        prepareProposal(void);
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);
        void                                        undoProposal(void);
        
    protected:
        
        void                                        swapNodeInternal(DagNode *, DagNode *newN);
        
        
    private:

        StochasticNode<RbVector<double> >*          variable; //!< The vector to operate on.
        
        double                                      lambda;  //!< The scale parameter of the proposal (larger lambda -> larger proposals).
        size_t                                      index;  //!< The index of the last modified element.
        double                                      storedValue;  //!< The stored value of the last modified element.

    };
    
}

#endif


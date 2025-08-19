#ifndef VectorScaleProposal_H
#define VectorScaleProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Scaling proposal of all elements of a Stochastic vector.
     *
     *
     * This proposal randomly scales all elements of a vector using the same scaling factor.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     * It generally makes more sense to apply the scaling proposal on a vector of positive
     * real numbers but technically it works on negative numbers too. However,
     * the proposal will never change the sign of the value and thus is incomplete if applied
     * to a variable defined on the whole real line.
     *
     * @see VectorSingleElementScaleProposal to scale a single element at a time.
     * @see RevLanguage::Move_VectorScale for the RL interface.
     */
    class VectorScaleProposal : public Proposal {
        
    public:
        VectorScaleProposal(StochasticNode<RbVector<double> >* n, double l);
        
        void                                        cleanProposal(void);
        VectorScaleProposal*                        clone(void) const;
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
        
        StochasticNode<RbVector<double> >*          variable;  //!< The vector to operate on.
        double                                      lambda;  //!< The scale parameter of the proposal (larger lambda -> larger proposals).
        size_t                                      length; //!< The size of the vector.
        double                                      storedScalingFactor; //!< The stored value of the last scaling factor.
                
    };
    
}

#endif


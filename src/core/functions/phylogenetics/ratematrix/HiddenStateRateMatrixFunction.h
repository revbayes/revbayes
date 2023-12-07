#ifndef HiddenStateRateMatrixFunction_H
#define HiddenStateRateMatrixFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Hidden state rate matrix function.
     *
     * The RevLanguage wrapper of the Free-K rate matrix connects
     * the variables/parameters of the function and creates the internal hiddenStateRateMatrixFunction object.
     * Please read the hiddenStateRateMatrixFunction.h for more info.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna & Will Freyman)
     *
     */
    class HiddenStateRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        HiddenStateRateMatrixFunction(bool r);
        virtual                                            ~HiddenStateRateMatrixFunction(void);                                      //!< Virtual destructor
        
        // public member functions
        HiddenStateRateMatrixFunction*                      clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        void                                                setObservedTransitionRates(const TypedDagNode< RbVector<double> > *tr);
        void                                                setObservedTransitionRates(const TypedDagNode< RbVector<RbVector<double> > > *tr);
        void                                                setObservedTransitionRates(const TypedDagNode< RevBayesCore::RateGenerator > *tr);
        void                                                setHiddenTransitionRates(const TypedDagNode< RbVector<double> > *tr);
        void                                                setHiddenTransitionRates(const TypedDagNode< RbVector<RbVector<double> > > *tr);
        void                                                setHiddenTransitionRates(const TypedDagNode< RevBayesCore::RateGenerator > *tr);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        bool                                                rescale;
        size_t                                              num_observed_states;
        size_t                                              num_hidden_states;
        const TypedDagNode<RbVector<RbVector<double> > >*   observed_transition_rates;
        const TypedDagNode<RbVector<double> >*              observed_transition_rates_flat;
        const TypedDagNode<RateGenerator>*                  observed_rate_generator;
        const TypedDagNode<RbVector<RbVector<double> > >*   hidden_transition_rates;
        const TypedDagNode<RbVector<double> >*              hidden_transition_rates_flat;
        const TypedDagNode<RateGenerator>*                  hidden_rate_generator;
        
    };
    
}

#endif

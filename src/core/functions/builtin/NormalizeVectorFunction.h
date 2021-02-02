#ifndef NormalizeVectorFunction_H
#define NormalizeVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /** @brief Function to normalize a vector of values so the sum of values equals a target
     *
     * @see RevLanguage::Func_fnNormalizedQuantile and RevLanguage::Func_normalize for the RL callers
     **/
    class NormalizeVectorFunction : public TypedFunction< RbVector<double> > {
        
    public:
        NormalizeVectorFunction(const TypedDagNode< RbVector<double> > * v, const TypedDagNode< double > * s);

        // public member functions
        NormalizeVectorFunction*                            clone(void) const;
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals; //!< vector to normalize
        const TypedDagNode< double >*                       sum; //!< target sum of values/size of the vector
        
    };
    
}


#endif

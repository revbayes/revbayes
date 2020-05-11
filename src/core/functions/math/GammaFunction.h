#ifndef GammaFunction_h
#define GammaFunction_h

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Calculates the gamma function for a real number.
     *
     * The gamma function is a generalization of the factorial function
     * to the complex numbers. 
     *
     * @note this should likely link to the gamma function to be clear
     * how this is calculated.
     *
     */
    class GammaFunction : public TypedFunction<double> {
        
    public:
        GammaFunction(const TypedDagNode<double> *a);
        
        GammaFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         a;
    };
}

#endif /* GammaFunction_h */

#ifndef LKJDistribution_H
#define LKJDistribution_H

#include <cstddef>

#include "TypedDistribution.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class LKJDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // LKJ distribution with parameter eta
        LKJDistribution(const TypedDagNode<double> *e, size_t d);
        
        virtual                                            ~LKJDistribution(void) {}
        
        // public member functions

        LKJDistribution*                                    clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // members
        
        const TypedDagNode<double>*                         eta;
        size_t                                              dim;
        
    };
    
}


#endif /* defined(LKJDistribution_H) */

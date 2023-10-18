#ifndef LKJPartialDistribution_H
#define LKJPartialDistribution_H

#include <cstddef>

#include "TypedDistribution.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class LKJPartialDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // LKJPartial distribution with parameter eta
        LKJPartialDistribution(const TypedDagNode<double> *e, size_t d);
        
        virtual                                            ~LKJPartialDistribution(void) {}
        
        // public member functions

        LKJPartialDistribution*                             clone(void) const;                                                          //!< Create an independent clone
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


#endif /* defined(LKJPartialDistribution_H) */

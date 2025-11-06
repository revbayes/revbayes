#ifndef __revbayes__GilbertGraphDistribution__
#define __revbayes__GilbertGraphDistribution__

#include <cstddef>

#include "TypedDistribution.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class GilbertGraphDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // LKJPartial distribution with parameter eta
        GilbertGraphDistribution(size_t n, const TypedDagNode<double> *p);
        
        virtual                                            ~GilbertGraphDistribution(void) {}
        
        // public member functions
        
        GilbertGraphDistribution*                           clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        
        const TypedDagNode<double>*                         edge_probability;
        size_t                                              dim;
        
    };
    
}


#endif /* defined(__revbayes__GilbertGraphDistribution__) */

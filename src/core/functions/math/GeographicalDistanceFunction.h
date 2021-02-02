#ifndef GeographicalDistanceFunction_hpp
#define GeographicalDistanceFunction_hpp

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Geographical distance be of a real number.
     *
     * Calculates the geographical distance between a pair of latitude
     * and longitudes.
     *
     * @note does this use the formula from pg. 792  of Landis et al.
     * 2013 paper?
     */
    class GeographicalDistanceFunction : public ContinuousFunction {
        
    public:
        GeographicalDistanceFunction(const TypedDagNode<double> *lat1, const TypedDagNode<double> *lon1, const TypedDagNode<double> *lat2, const TypedDagNode<double> *lon2);
        
        GeographicalDistanceFunction*       clone(void) const;                                                  //!< Create a clone.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         latitude1;
        const TypedDagNode<double>*         longitude1;
        const TypedDagNode<double>*         latitude2;
        const TypedDagNode<double>*         longitude2;
    };
}


#endif /* GeographicalDistanceFunction_hpp */

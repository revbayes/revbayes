#ifndef SortVectorFunction_H
#define SortVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Sort Vector Function
     *
     * This class implements a function that sorts a vector of double type values
     * This class has these member parameters:
     * @param vector The vector to be sorted
     * @param ascending a boolean denoting whether the vector should be sorted in ascending order
     *
     *
     */
    class SortVectorFunction : public TypedFunction< RbVector<double> > {
        
    public:
        SortVectorFunction(const TypedDagNode< RbVector<double> >* vec, bool ascending = true);                                              //!< Basic constructor
        
        virtual                                         ~SortVectorFunction(void);                                       //!< Virtual destructor for derived classes
        
        // public member functions
        SortVectorFunction*                             clone(void) const;                                                      //!< Create a clone
        void                                            update(void);                                                           //!< Update the value of the function
        
    protected:
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap parameters
        
    private:
        
        // members
        bool                                            ascending;
        const TypedDagNode< RbVector<double> >*         realVector;                                                          //!< Real positive vector parameter
    };
    
}


#endif

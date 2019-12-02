#ifndef SimplexFromVectorFunction_H
#define SimplexFromVectorFunction_H

#include "Simplex.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Simplex constructor from vector
     *
     * This class implements a function that constructs a simplex from a vector of
     * real positive values.
     */
    class SimplexFromVectorFunction : public TypedFunction< Simplex > {
        
    public:
        SimplexFromVectorFunction(const TypedDagNode< RbVector<double> >* vec);                                              //!< Basic constructor
        
        virtual                                         ~SimplexFromVectorFunction(void);                                       //!< Virtual destructor for derived classes
        
        // public member functions
        SimplexFromVectorFunction*                      clone(void) const;                                                      //!< Create a clone
        void                                            update(void);                                                           //!< Update the value of the function
        
    protected:
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*         realPosVector;                                                          //!< Real positive vector parameter
        
    };
    
}


#endif

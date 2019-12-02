#ifndef SortVectorFunction_H
#define SortVectorFunction_H

#include "RbVector.h"
#include "Simplex.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Simplex constructor from vector
     *
     * This class implements a function that constructs a simplex from a vector of
     * real positive values.
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

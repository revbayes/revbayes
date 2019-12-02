#ifndef ContinuousFunction_H
#define ContinuousFunction_H

#include "TypedFunction.h"

#include <iostream>

namespace RevBayesCore {
    
    /**
     * @brief Continuous function interface.
     *
     * A continuous function is defined to be continuously on the real line.
     * We use this interface if we want to allow only continuous functions.
     *
     */
    class ContinuousFunction : public TypedFunction<double> {
        
    public:
        // constructors and destructor
        virtual                            ~ContinuousFunction(void) {}
                
        // pure virtual public methors
        virtual ContinuousFunction*         clone(void) const = 0;                                  //!< Clone the function
        virtual void                        update(void) = 0;                                       //!< Clone the function
        
    protected:
        ContinuousFunction(double *v);
        
    };
        
}

#endif

#ifndef Func_absVectorInt_H
#define Func_absVectorInt_H

#include "ModelVector.h"
#include "Natural.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the absolute value function for vectors of integers.
     */
    class Func_absVectorInt : public TypedFunction< ModelVector<Natural> > {
        
    public:
        Func_absVectorInt( void );
        
        // Basic utility functions
        Func_absVectorInt*                                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<std::int64_t> >* createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
    
    
    protected:
        
    };
    
}

#endif

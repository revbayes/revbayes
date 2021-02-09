#ifndef Func__simplexIndexOperator_H
#define Func__simplexIndexOperator_H

#include "RlTypedFunction.h"
#include "Probability.h"

#include <string>

namespace RevLanguage {
    
    class Func__simplexIndexOperator : public TypedFunction<Probability> {
        
    public:
        Func__simplexIndexOperator( void );
        
        // Basic utility functions
        Func__simplexIndexOperator*                                     clone(void) const;                              //!< Clone the object
        static const std::string&                                       getClassType(void);                             //!< Get class name
        static const TypeSpec&                                          getClassTypeSpec(void);                         //!< Get class type spec
        std::string                                                     getFunctionName(void) const;                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                                 getTypeSpec(void) const;                        //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*                            createFunction(void) const ;                    //!< Create a new internal function object
        const ArgumentRules&                                            getArgumentRules(void) const;                   //!< Get argument rules
        
    private:
        
    };
    
}


#endif /* Func__simplexIndexOperator_H */

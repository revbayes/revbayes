#ifndef Func_meanPositive_H
#define Func_meanPositive_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the arithmetic meanPositive function.
     *
     * The RevLanguage wrapper of the meanPositive function connects
     * the variables/parameters of the function and creates the internal MeanFunction object.
     *
     * @copybrief RevBayesCore::MeanFunction
     * @see RevBayesCore::MeanFunction for the internal object
     */
    class Func_meanPositive : public TypedFunction<RealPos> {
        
    public:
        Func_meanPositive( void );
        
        // Basic utility functions
        Func_meanPositive*                              clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif

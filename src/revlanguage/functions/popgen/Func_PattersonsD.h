#ifndef Func_PattersonsD_H
#define Func_PattersonsD_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of Patterson's D function.
     *
     * The RevLanguage wrapper of Patterson's D function connects
     * creates the internal PattersonsDFunction object.
     * Please read the PattersonsDFunction.h for more info.
     *
     *
     * Will Freyman 8/25/15
     *
     */
    class Func_PattersonsD : public TypedFunction<Real> {
        
    public:
        Func_PattersonsD( void );
        
        // Basic utility functions
        Func_PattersonsD*                           clone(void) const;                                          //!< Clone the object
        static const std::string&                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< double >*      createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
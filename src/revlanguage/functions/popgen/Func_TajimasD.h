#ifndef Func_TajimasD_H
#define Func_TajimasD_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of Tajima's D function.
     *
     * The RevLanguage wrapper of Tajima's D function connects
     * the variables/parameters of the function and creates the internal TajimasDFunction object.
     * Please read the TajimasDFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-04-30, version 1.0
     *
     */
    class Func_TajimasD : public TypedFunction<Real> {
        
    public:
        Func_TajimasD( void );
        
        // Basic utility functions
        Func_TajimasD*                              clone(void) const;                                          //!< Clone the object
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

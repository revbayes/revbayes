#ifndef Func_TajimasPi_H
#define Func_TajimasPi_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of Tajima's pi function.
     *
     * The RevLanguage wrapper of Tajima's pi function connects
     * the variables/parameters of the function and creates the internal TajimasPiFunction object.
     * Please read the TajimasPiFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-04-30, version 1.0
     *
     */
    class Func_TajimasPi : public TypedFunction<RealPos> {
        
    public:
        Func_TajimasPi( void );
        
        // Basic utility functions
        Func_TajimasPi*                             clone(void) const;                                          //!< Clone the object
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

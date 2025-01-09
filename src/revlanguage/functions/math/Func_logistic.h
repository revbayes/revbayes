#ifndef Func_logistic_H
#define Func_logistic_H

#include "Probability.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the logistic function.
     *
     * The RevLanguage wrapper of the logistic function connects
     * the variables/parameters of the function and creates the internal LogisticFunction object.
     * Please read the LogisticFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2015-01-31, version 1.0
     *
     */
    class Func_logistic :  public TypedFunction<Probability> {
        
    public:
        Func_logistic( void );
        
        // Basic utility functions
        Func_logistic*                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
	std::vector<std::string>                        getFunctionNameAliases( void ) const;                       //!< Get aliases of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif

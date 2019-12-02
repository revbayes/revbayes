#ifndef Func_probability_H
#define Func_probability_H

#include "Probability.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage conversion function between a Real and a Probability.
     *
     * This Rev functions converts between Rev types from Real to Probability.
     * The Rev conversion functions only checks during instantiation if the conversion
     * is possible. We then simply return the same DAG node but in a different Rev wrapper.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-31, version 1.0
     *
     */
    class Func_probability :  public TypedFunction<Probability> {
        
    public:
        Func_probability( void );
        
        // Basic utility functions
        Func_probability*                               clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules                              
        std::vector<std::string>                        getHelpAuthor(void) const; 
        std::string                                     getHelpDescription(void) const; 
        std::string                                     getHelpExample(void) const;
        
    };
    
}

#endif

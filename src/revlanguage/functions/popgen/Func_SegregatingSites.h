#ifndef Func_SegregatingSites_H
#define Func_SegregatingSites_H

#include "Natural.h"
#include "RlTypedFunction.h"

#include <map>
#include <cstdint>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the Segregating-Sites function.
     *
     * The RevLanguage wrapper of the Segregating-Sites function connects
     * the variables/parameters of the function and creates the internal SegregatingSitesFunction object.
     * Please read the SegregatingSitesFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-04-30, version 1.0
     *
     */
    class Func_SegregatingSites : public TypedFunction<Natural> {
        
    public:
        Func_SegregatingSites( void );
        
        // Basic utility functions
        Func_SegregatingSites*                      clone(void) const;                                          //!< Clone the object
        static const std::string&                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< std::int64_t >*        createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif

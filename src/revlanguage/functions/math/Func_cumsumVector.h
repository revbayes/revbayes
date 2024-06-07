#ifndef Func_cumsumVector_H
#define Func_cumsumVector_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the cumulative sum function.
     *
     * The RevLanguage wrapper of the cumulative sum function connects
     * the variables/parameters of the function and creates the internal CumulativeSumVectorFunction object.
     * Please read the CumulativeSumVectorFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_cumsumVector : public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_cumsumVector( void );
        
        // Basic utility functions
        Func_cumsumVector*                                                  clone(void) const;                                          //!< Clone the object
        static const std::string&                                           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                     getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >*       createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                getArgumentRules(void) const;                               //!< Get argument rules
    
    
    protected:
        
    };
    
}

#endif

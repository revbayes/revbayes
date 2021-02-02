#ifndef Func_absVector_H
#define Func_absVector_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the absolute value function.
     *
     * The RevLanguage wrapper of the abs function connects
     * the variables/parameters of the function and creates the internal AbsoluteValueFunction object.
     * Please read the AbsoluteValueFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_absVector : public TypedFunction< ModelVector<RealPos> > {
        
    public:
        Func_absVector( void );
        
        // Basic utility functions
        Func_absVector*                                                     clone(void) const;                                          //!< Clone the object
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

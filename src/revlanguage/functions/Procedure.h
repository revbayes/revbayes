#ifndef Procedure_H
#define Procedure_H

#include "RlFunction.h"

namespace RevLanguage {
    
    class Procedure : public Function {
        
    public:
        virtual                                            ~Procedure(void);                                                                   //!< Destructor
        
        virtual RevPtr<RevVariable>                         execute(void) = 0;                                                                  //!< Create a random variable from this distribution
        virtual std::string                                 getFunctionName(void) const = 0;                                                    //!< Get the name of the function
        virtual const ArgumentRules&                        getArgumentRules(void) const = 0;                                                   //!< Get argument rules
        virtual const TypeSpec&                             getReturnType(void) const = 0;                                                      //!< Get type of return value

        // Basic utility functions you have to override
        virtual Procedure*                                  clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec
        
        
        
        
    protected:
        Procedure(void);                                                                                                                    //!< Basic constructor
        

    };
    
}


#endif



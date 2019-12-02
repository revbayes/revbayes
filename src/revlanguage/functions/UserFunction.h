#ifndef UserFunction_H
#define UserFunction_H

#include "ArgumentRules.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "UserFunctionDef.h"

#include <list>
#include <set>
#include <string>


namespace RevLanguage {

    /**
     * @brief Class holding user functions
     *
     * The UserFunction class holds information about a user function and is used both to execute
     * the function and make deterministic variables using the function.
     */
    class UserFunction :  public Function {
        
    public:
        UserFunction(UserFunctionDef* def);                                                                     //!< Standard concstructor
        
        virtual                                    ~UserFunction() {}                                           //!< Virtual destructor

        // Basic utility functions
        UserFunction*                               clone(void) const;                                          //!< Clone the object
        static const std::string&                   getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get language type of the object

        // Regular functions   
        virtual RevPtr<RevVariable>                 execute(void);                                              //!< Execute function (create det var)
        virtual RevPtr<RevVariable>                 executeCode(void);                                          //!< Execute function (execute code)
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        std::vector<const RevBayesCore::DagNode*>   getParameters(void) const;                                  //!< Get parameters
        const TypeSpec&                             getReturnType(void) const;                                  //!< Get type of return value
        
    protected:
        // Member variable
        RevPtr<UserFunctionDef>                     function_def;                                                //!< The function definition
    };
    
}

#endif


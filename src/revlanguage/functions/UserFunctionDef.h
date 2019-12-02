#ifndef UserFunctionDef_H
#define UserFunctionDef_H

#include "SyntaxElement.h"

#include <list>

namespace RevLanguage {

    /**
     * @brief Definitions of user functions or procedures
     *
     * The UserFunctionDef class holds the code, argument rules and return type of a user-defined Rev function 
     * or procedure. It supports reference counting so that UserFunction and UserProcedure can store smart
     * pointers to their definitions. The code in UserFunctionDef is used by UserFunction and UserProcedure
     * when executing the function, while the function execution frame is created and maintained by UserFunction
     * or UserProcedure.
     */
    class UserFunctionDef {
        
    public:
        UserFunctionDef(const ArgumentRules* argRules, const TypeSpec& retType, std::list<SyntaxElement*>* stmts, std::string name);  //!< Standard concstructor
        UserFunctionDef(const UserFunctionDef& x);                                                                  //!< Copy constructor

        virtual                                    ~UserFunctionDef();                                      //!< Virtual destructor
        
        // Assignment operator
        UserFunctionDef&                           operator=(const UserFunctionDef& f);                     //!< Assignment operator
        
        // Basic utility functions
        int                                         decrementReferenceCount(void);                          //!< Decrement reference count
        const ArgumentRules&                        getArgumentRules(void) const;                           //!< Get argument rules
        const std::list<SyntaxElement*>&            getCode(void) const;                                    //!< Get a pointer to the code
        const std::string&                          getName(void) const;
        const TypeSpec&                             getReturnType(void) const;                              //!< Get type of return value
        int                                         getReferenceCount(void) const;                          //!< Get reference count
        void                                        incrementReferenceCount(void);                          //!< Increment reference count
        
    protected:
        // Member variables
        const ArgumentRules*                        argumentRules;                                          //!< The argument rules
        std::list<SyntaxElement*>*                  code;                                                   //!< The code
        std::string                                 name;                                                   //!< The function name
        int                                         referenceCount;                                         //!< Reference count
        const TypeSpec                              return_type;                                             //!< The return type (complete specification)
    };
    
}

#endif


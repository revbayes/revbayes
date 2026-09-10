#ifndef SyntaxFormal_H
#define SyntaxFormal_H

#include "SyntaxElement.h"

#include <iostream>
#include <list>

namespace RevLanguage {

    class ArgumentRule;
    class ArgumentEnvironment;
    class TypeSpec;

    /**
     * @brief Formal arguments
     *
     * This syntax element class deals with formal arguments, that is, 
     * specifications of the type and name of arguments in function
     * definitions. Formals are also used to declare member variables in
     * class definitions.
     *
     * Note that formal argument specifications can only occur in the user
     * workspace. They cannot occur in local environments, such as inside
     * function definitions.
     */
    class SyntaxFormal : public SyntaxElement {

    public:
        SyntaxFormal(const std::string& lbl, SyntaxElement* defaultVal);                                //!< Implicit type
        SyntaxFormal(const std::string& type, const std::string& lbl, SyntaxElement* defaultVal);       //!< Explicit type or type modifier (const, mutable etc)
        SyntaxFormal(const SyntaxFormal& x);                                                            //!< Copy constructor
	    
        virtual                    ~SyntaxFormal();                                                     //!< Destructor

        // Assignment operator
        SyntaxFormal&               operator=(const SyntaxFormal& x);                                   //!< Assignment operator

        // Basic utility functions
        SyntaxFormal*               clone() const;                                                      //!< Clone object
        
        // Regular functions
        const ArgumentRule*         getArgumentRule(void) const;                                        //!< Get the argument rule
        ArgumentRule*               getArgumentRule(void);                                              //!< Get the argument rule (non-const to return non-const rule)
        const std::string&          getLabel(void) const;                                               //!< Get label
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value
        void                        setIsProtected( bool prot = true );                                 //!< Set protected flag of the formal argument specification
        
    protected:
        
        ArgumentRule*               argRule;                                                            //!< The argument rule version of the formal
    
    };
    
}

#endif

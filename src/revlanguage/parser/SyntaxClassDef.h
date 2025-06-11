#ifndef SyntaxClassDef_H
#define SyntaxClassDef_H

#include "SyntaxElement.h"
#include "SyntaxFormal.h"
#include "SyntaxVariable.h"

#include <ostream>
#include <list>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Class definition syntax element
     *
     * This syntax element is used to hold Rev class definitions.
     * These consist of a class name, a base class name, and a list
     * of definitions of member variables and member functions.
     *
     * Note that the Rev language specification only allows definition
     * of classes in the global namespace. That is, class definitions
     * within functions or within other class definitions are not
     * supported.
     */
    class SyntaxClassDef : public SyntaxElement {

    public:
        SyntaxClassDef(const std::string& name, const std::string& base, std::list<SyntaxElement*>* defs);          //!< Standard constructor
        SyntaxClassDef(const SyntaxClassDef& x);                                                                    //!< Copy constructor
    
	    virtual                    ~SyntaxClassDef();                                                               //!< Destructor

        // Assignment operator
        SyntaxClassDef&             operator=(const SyntaxClassDef& x);                                             //!< Assignment operator

        // Basic utility functions
        SyntaxElement*              clone() const;                                                                  //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);   //!< Get semantic value

    protected:
        std::string                 className;                                                                      //!< The name of the class
        std::string                 baseClass;                                                                      //!< The name of the base class
        std::list<SyntaxElement*>*  definitions;                                                                    //!< Defined member variables and functions
    
    };
    
}

#endif


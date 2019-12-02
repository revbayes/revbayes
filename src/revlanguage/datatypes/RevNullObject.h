#ifndef RevNullObject_H
#define RevNullObject_H

#include "RevVariable.h"
#include "WorkspaceObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {

class RevNullObject : public WorkspaceObject {
    
public:
    RevNullObject(void);    
    virtual                            ~RevNullObject(void) {}                                                                      //!< Virtual destructor
    
    bool                                operator==(void* ptr);                                                                      //!< Compare to NULL pointer
    bool                                operator==(const RevNullObject& obj);                                                       //!< Compare to Rev NULL object
    bool                                operator==(const RevObject& obj);                                                           //!< Compare to Rev object
    
    static RevNullObject&     getInstance(void)                                                                                     //!< Get the global instance
    {
        static RevNullObject theInstance = RevNullObject();
        return theInstance;
    }

    
    // Basic utility functions you have to override (also getClassTypeSpec()!)
    RevNullObject*                      clone(void) const;                                                                          //!< Clone object
    RevPtr<RevVariable>                 executeMethod(const std::string &name, const std::vector<Argument> &args, bool &found);     //!< Execute the member method
    static const std::string&           getClassType(void);                                                                         //!< Get Rev type
    static const TypeSpec&              getClassTypeSpec(void);                                                                     //!< Get class type spec
    const TypeSpec&                     getTypeSpec(void) const;                                                                    //!< Get the type spec for this instance
    void                                printValue(std::ostream& o, bool user) const;                                               //!< Print value for user
        
    
};

// Global functions using the class
std::ostream&                           operator<<(std::ostream& o, const RevNullObject& x);                    //!< Overloaded output operator

}

#endif

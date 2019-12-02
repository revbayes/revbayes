#ifndef RevMemberObject_H
#define RevMemberObject_H

#include <string>

#include "MethodTable.h"
//#include "RevVariable.h"

namespace RevLanguage {
    
    class Argument;
    class RevVariable;
    
    /**
     * @brief Rev interface for member object which provide methods.
     *
     * This interface specifies that classes implementing the interface will provide member methods.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-10-23
     *
     */
    class RevMemberObject {
        
    public:
        
        virtual RevPtr<RevVariable>             executeMethod(const std::string &name, const std::vector<Argument> &args, bool &found) = 0;     //!< Execute the member method
        virtual const MethodTable&              getMethods(void) const = 0;                                                                     //!< Get member methods
	virtual ~RevMemberObject() {}
    };
    
}

#endif


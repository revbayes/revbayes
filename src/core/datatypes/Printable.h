#ifndef Printable_H
#define Printable_H

#include <ostream>
#include <string>
#include "json.h"

namespace RevBayesCore {
    
    class Printable {
        
    public:
        virtual                         ~Printable(void) {}
        
	virtual json                    toJSON() const = 0;
        virtual void                    printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const = 0;                                         //!< print object for user (in user-formatted way)
        virtual void                    printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const = 0;           //!< print object with rounding
        virtual void                    printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const = 0;          //!< print object with maximum precision

    };
    
}

#endif


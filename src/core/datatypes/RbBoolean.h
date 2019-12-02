#ifndef RbBoolean_H
#define RbBoolean_H

#include <ostream>

namespace RevBayesCore {
    
    /**
     * RevBayes wrapper class for boolean.
     *
     * This class is essentially a wrapper of a bool.
     * We only need this as a workaround for std::vector<bool>.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-01, version 1.0
     */
    class Boolean {
        
    public:
                        Boolean() {}                                                                                      //!< Constructor requires character type
                        Boolean(bool b) : value( b ) {}                                                                   //!< Constructor requires character type
        virtual        ~Boolean(void) {}
        
                        operator bool() const { return value; }
        bool            operator==(const Boolean &rm) const { return value == rm.value; }
        bool            operator==(const bool &rm) const { return value == rm; }
        bool            operator!=(const Boolean &rm) const { return !operator==(rm); }
        bool            operator<(const Boolean &rm) const { return value < rm.value; }
        bool            operator<=(const Boolean &rm) const { return operator<(rm) || operator==(rm); }
        
        bool            value;
        
    protected:
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const Boolean& x);                    //!< Overloaded output operator
    
}

#endif

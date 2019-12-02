#ifndef RbHelpType_H
#define RbHelpType_H

#include "RbHelpEntry.h"
#include "RbHelpFunction.h"
#include "RbHelpReference.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for a type.
     *
     * A help entry for a type provides all the help information available on this type, e.g., 
     * about member methods.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     *
     */
    class RbHelpType : public RbHelpEntry {
        
    public:
        
        virtual RbHelpType*                         clone(void) const;                                                   //!< Clone this object
        
        // getter
        const std::vector<RbHelpFunction>&          getConstructors(void) const;
        const std::vector<RbHelpFunction>&          getMethods(void) const;
        
        // setters
        void                                        setConstructors(const std::vector<RbHelpFunction>& c);
        void                                        setMethods(const std::vector<RbHelpFunction>& c);
        
        
    private:
        
        std::vector<RbHelpFunction>                 constructors;
        std::vector<RbHelpFunction>                 methods;

        
    };
    
    
}

#endif


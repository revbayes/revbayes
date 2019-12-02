#ifndef RbHelpMonitor_H
#define RbHelpMonitor_H

#include "RbHelpType.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for a Monitor.
     *
     * A help entry for a Monitor provides all the help information available on this Monitor, e.g.,
     * about the constructor and member methods.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-10-10
     *
     */
    class RbHelpMonitor : public RbHelpType {
        
    public:
        
        RbHelpMonitor();
        
        virtual RbHelpMonitor*                      clone(void) const;                                                   //!< Clone this object

        
    private:
        
        
        
    };
    
    
}

#endif


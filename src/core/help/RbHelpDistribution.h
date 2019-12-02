#ifndef RbHelpDistribution_H
#define RbHelpDistribution_H

#include "RbHelpType.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for a distribution.
     *
     * A help entry for a distribution provides all the help information available on this distribution, e.g.,
     * about the constructor and member methods.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-10-10
     *
     */
    class RbHelpDistribution : public RbHelpType {
        
    public:
        
        RbHelpDistribution();
        
        virtual RbHelpDistribution*                     clone(void) const;                                                   //!< Clone this object

    private:
        
        
        
    };
    
    
}

#endif


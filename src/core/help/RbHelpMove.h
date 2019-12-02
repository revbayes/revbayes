#ifndef RbHelpMove_H
#define RbHelpMove_H

#include "RbHelpType.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for a Move.
     *
     * A help entry for a Move provides all the help information available on this Move, e.g.,
     * about the constructor and member methods.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-10-10
     *
     */
    class RbHelpMove : public RbHelpType {
        
    public:
        
        RbHelpMove();
        
        virtual RbHelpMove*                         clone(void) const;                                                   //!< Clone this object

        
    private:
        
        
        
    };
    
    
}

#endif


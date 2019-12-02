/**
 * @file
 * This file contains the declaration of Argument, which is
 * used to hold a potentially labeled argument passed to a
 * function. The Argument class just holds a DAG node pointer,
 * it does not manage the DAG node itself.
 *
 * @brief Declaration of Argument
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#ifndef Argument_H
#define Argument_H

#include "RevPtr.h"
#include "RevVariable.h"

#include <ostream>
#include <string>

namespace RevLanguage {

class Argument {

    public:

    Argument(const RevPtr<RevVariable> &arg, const std::string& argLabel = "", bool isConst = true);                   //!< Constructor
    
    // Basic utility functions
    Argument*                           clone(void) const;                                                          //!< Clone object
    void                                printValue(std::ostream& o) const;                                          //!< Complete info about object
    
    // Regular functions
    const std::string&                  getLabel(void) const;
    RevPtr<RevVariable>&                getVariable(void);                                                          //!< Get the variable contained in this argument
    RevPtr<const RevVariable>           getVariable(void) const;                                                    //!< Get the variable contained in this argument
    const RevPtr<RevVariable>&          getReferenceVariable(void) const;
    bool                                isConstant(void) const;
    void                                setLabel(const std::string&);

    private:
    
    RevPtr<RevVariable>                 var;                                                                        //!< Pointer to the variable slot containing the variable (and value)
    bool                                isConst;                                                                    //!< Is this a constant variable?
    std::string                         label;                                                                      //!< Label of argument
    
};
    
}

#endif


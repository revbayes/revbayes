/**
 * @file
 * This file contains the implementation of Argument, which is
 * used to hold a potentially labeled argument passed to a
 * function.
 *
 * @brief Implementation of Argument
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

#include <sstream>

#include "Argument.h"
#include "RbException.h"
#include "RevPtr.h"
#include "RevVariable.h"

using namespace RevLanguage;


/** Construct from argument label and DAG node */
Argument::Argument(const RevPtr<RevVariable> &v, const std::string& argLabel, bool c) :
    var( v ),
    isConst( c ),
    label( argLabel )
{
    
}


const std::string& Argument::getLabel( void ) const
{
    return label;
}


Argument* Argument::clone( void ) const
{
    return new Argument( *this );
}


const RevPtr<RevVariable>& Argument::getReferenceVariable(void) const
{
    if ( isConst )
    {
        throw RbException("Illegal attempt to get a non-const reference to a constant variable.");
    }
    return var;
}


RevPtr<RevVariable>& Argument::getVariable(void)
{
    return var;
}


RevPtr<const RevVariable> Argument::getVariable(void) const
{
    return RevPtr<const RevVariable>( var );
}


bool Argument::isConstant( void ) const
{
    return isConst;
}



/** Complete info about object */
void Argument::printValue(std::ostream &o) const
{
    
    o << label << " = ";
    var->printValue(o, true);
    
}


void Argument::setLabel( const std::string& l)
{
    label = l;
}

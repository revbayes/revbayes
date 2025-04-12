/**
 * @file
 * This file contains the implementatin of OptionRule, which is
 * used to describe argument rules corresponding to the
 * selection of one of several RlString options.
 *
 * @brief Implementation of OptionRule
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#include "OptionRule.h"

#include <cstddef>
#include <sstream>
#include <string>

#include "RbException.h"
#include "RlString.h"
#include "Argument.h"
#include "DagNode.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

using namespace RevLanguage;

/** Construct rule without default value; use "" for no label. */
OptionRule::OptionRule( const std::string& argName, const std::vector<std::string>& optVals, const std::string& argDesc ) : ArgumentRule( argName, RlString::getClassTypeSpec(), argDesc, BY_VALUE, ANY ),
    options( optVals )
{

    if ( areOptionsUnique( optVals ) == false )
    {
        throw RbException( "Options are not unique" );
    }
    
}


/** Construct rule with default value; use "" for no label. */
OptionRule::OptionRule(const std::string& argName, RlString* defVal, const std::vector<std::string>& optVals, const std::string& argDesc  ) : ArgumentRule( argName, RlString::getClassTypeSpec(), argDesc, BY_VALUE, ANY, defVal ),
    options( optVals )
{

    if ( areOptionsUnique( optVals ) == false )
    {
        throw RbException( "Options are not unique" );
    }
    
}


/** Help function to test whether a std::vector of std::string contains unique string values */
bool OptionRule::areOptionsUnique( const std::vector<std::string>& optVals ) const
{

    for ( size_t i = 0; i < optVals.size(); i++ )
    {
        for ( size_t j = i + 1; j < optVals.size(); j++ )
        {
            if ( optVals[i] == optVals[j] )
            {
                return false;
            }
        }
    }

    return true;
}



OptionRule* OptionRule::clone( void ) const
{
    
    return new OptionRule( *this );
}



const std::vector<std::string>& OptionRule::getOptions( void ) const
{
    // return a const reference to the internal value
    return options;
}


double OptionRule::isArgumentValid( Argument &arg ) const
{
    RevPtr<RevVariable> the_var = arg.getVariable();
    if ( the_var == NULL )
    {
        return -1;
    }
    
    bool once = false;
    if (not the_var->getRevObject().hasDagNode() or the_var->getRevObject().getDagNode()->isConstant())
    {
        once = true;
    }
    else if ( evalType == BY_VALUE || the_var->isWorkspaceVariable() || the_var->getRevObject().isConstant())
    {
        once = true;
    }
    
    RlString *revObj = dynamic_cast<RlString *>( &the_var->getRevObject() );
    if ( revObj != NULL )
    {
        
        const std::string &argValue = revObj->getValue();
        for ( std::vector<std::string>::const_iterator it = options.begin(); it != options.end(); ++it)
        {
            
            if ( argValue == *it )
            {
                return 0.0;
            }
        }
        return -1;
    }
    else
    {
        return -1;
    }

    
}


/** Print value for user */
void OptionRule::printValue(std::ostream& o) const
{

    ArgumentRule::printValue(o);

    o << " {valid options: ";
    for (std::vector<std::string>::const_iterator it = options.begin(); it != options.end(); ++it)
    {
        
        if ( it != options.begin() )
        {
            o << "|";
        }
        
        o << "\"" << *it << "\"";

    }
    o << "}"; // << std::endl;
}


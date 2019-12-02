#include "RbHelpArgument.h"

using namespace RevBayesCore;

const std::string& RbHelpArgument::getArgumentDagNodeType( void ) const
{
    // return a const reference to the internal value
    return argumentDagNodeType;
}


const std::string& RbHelpArgument::getArgumentPassingMethod( void ) const
{
    // return a const reference to the internal value
    return argumentPassingMethod;
}


const std::string& RbHelpArgument::getDefaultValue( void ) const
{
    // return a const reference to the internal value
    return defaultValue;
}


const std::string& RbHelpArgument::getDescription( void ) const
{
    // return a const reference to the internal value
    return description;
}


const std::string& RbHelpArgument::getLabel( void ) const
{
    // return a const reference to the internal value
    return label;
}


const std::vector<std::string>& RbHelpArgument::getOptions( void ) const
{
    // return a const reference to the internal value
    return options;
}


const std::string& RbHelpArgument::getValueType( void ) const
{
    // return a const reference to the internal value
    return valueType;
}


void RbHelpArgument::setArgumentDagNodeType(const std::string &dt)
{
    // overwrite the internal value with the new value
    argumentDagNodeType = dt;
}


void RbHelpArgument::setArgumentPassingMethod(const std::string &pm)
{
    // overwrite the internal value with the new value
    argumentPassingMethod = pm;
}



void RbHelpArgument::setDefaultValue(const std::string &df)
{
    // overwrite the internal value with the new value
    defaultValue = df;
}


void RbHelpArgument::setDescription(const std::string &d)
{
    // overwrite the internal value with the new value
    description = d;
}


void RbHelpArgument::setLabel(const std::string &l)
{
    // overwrite the internal value with the new value
    label = l;
}


void RbHelpArgument::setOptions(const std::vector<std::string> &o)
{
    // overwrite the internal value with the new value
    options = o;
}


void RbHelpArgument::setValueType(const std::string &vt)
{
    // overwrite the internal value with the new value
    valueType = vt;
}

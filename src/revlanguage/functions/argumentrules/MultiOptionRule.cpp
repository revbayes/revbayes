#include "MultiOptionRule.h"

#include <algorithm>

#include "Argument.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlString.h"

using namespace RevLanguage;


MultiOptionRule::MultiOptionRule(const std::string& argName, RlString* defVal, const std::vector<std::string>& optVals, const std::string& sep, const std::string& argDesc ) : OptionRule( argName, defVal, optVals, argDesc ),
    separator( sep )
{
    argTypeSpecs.push_back( ModelVector<RlString>::getClassTypeSpec() );

    // defaults are not fitted, so store the default as the String[] that fitArgument would make
    std::vector<std::string> tokens;
    selected( *defVal, tokens );
    RevBayesCore::RbVector<std::string> v( tokens );
    defaultVar = new RevVariable( new ModelVector<RlString>( v ) );
}


MultiOptionRule* MultiOptionRule::clone( void ) const
{
    return new MultiOptionRule( *this );
}


Argument MultiOptionRule::fitArgument( Argument& arg ) const
{
    std::vector<std::string> tokens;
    selected( arg.getVariable()->getRevObject(), tokens );
    RevBayesCore::RbVector<std::string> v( tokens );
    RevPtr<RevVariable> var = new RevVariable( new ModelVector<RlString>( v ), arg.getVariable()->getName() );

    return Argument( var, arg.getLabel(), true );
}


double MultiOptionRule::isArgumentValid( Argument &arg ) const
{
    RevPtr<RevVariable> the_var = arg.getVariable();
    std::vector<std::string> tokens;
    if ( the_var == NULL || selected( the_var->getRevObject(), tokens ) == false || tokens.empty() )
    {
        return -1;
    }

    for ( size_t i = 0; i < tokens.size(); ++i )
    {
        if ( std::find( options.begin(), options.end(), tokens[i] ) == options.end() )
        {
            return -1;
        }
        if ( std::find( tokens.begin(), tokens.begin() + i, tokens[i] ) != tokens.begin() + i )
        {
            return -1;
        }
    }

    return 0.0;
}


bool MultiOptionRule::selected(const RevObject& value, std::vector<std::string>& tokens) const
{
    tokens.clear();
    std::vector<std::string> strings;
    if ( value.isType( RlString::getClassTypeSpec() ) )
    {
        strings.push_back( static_cast<const RlString&>( value ).getValue() );
    }
    else if ( value.isType( ModelVector<RlString>::getClassTypeSpec() ) )
    {
        const RevBayesCore::RbVector<std::string>& v = static_cast<const ModelVector<RlString>&>( value ).getValue();
        for ( size_t i = 0; i < v.size(); ++i )
        {
            strings.push_back( v[i] );
        }
    }
    else
    {
        return false;
    }

    for ( const std::string& s : strings )
    {
        size_t start = 0;
        size_t pos;
        while ( (pos = s.find(separator, start)) != std::string::npos )
        {
            tokens.push_back( s.substr(start, pos - start) );
            start = pos + separator.size();
        }
        tokens.push_back( s.substr(start) );
    }

    return true;
}


#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Move.h"
#include "Natural.h"
#include "Probability.h"
#include "RevObject.h"
#include "RealPos.h"
#include "RlMove.h"
#include "TypeSpec.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "StringUtilities.h"
#include "WorkspaceObject.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

Move::Move(void) : WorkspaceToCoreWrapperObject<RevBayesCore::Move>()
{
    
}


/** Get Rev type of object */
const std::string& Move::getClassType(void)
{
    
    static std::string rev_type = "Move";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::Move>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the aliases of the Rev name for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Move::getConstructorFunctionAliases( void ) const
{
    // create a constructor function name alias variable that is the same for all instance of this class
    std::vector<std::string> aliases;
    
    std::vector<std::string> move_aliases = getMoveAliases();
    for (size_t i=0; i < move_aliases.size(); ++i)
    {
        std::string tmp = move_aliases[i];
        std::string c_name = "mv" + StringUtilities::firstCharToUpper( tmp );
        
        aliases.push_back( c_name );
    }
    
    return aliases;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string tmp = getMoveName();
    std::string c_name = "mv" + StringUtilities::firstCharToUpper( tmp );
    
    return c_name;
}


/** Return member rules */
const MemberRules& Move::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "weight", RealPos::getClassTypeSpec(), "The weight how often on average this move will be used per iteration.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos( 1.0 ) ) );
        move_member_rules.push_back( new ArgumentRule( "delay", Natural::getClassTypeSpec(), "The of iterations before using.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( 0L ) ) );
        move_member_rules.push_back( new ArgumentRule( "weight", RealPos::getClassTypeSpec(), "The weight determines the relative frequency with which this move will be attempted. For details, see the description of the 'moveschedule' parameter on the documentation page for 'mcmc()'.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos( 1.0 ) ) );
        move_member_rules.push_back( new ArgumentRule( "tuneTarget", Probability::getClassTypeSpec(), "The acceptance probability targeted by auto-tuning.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability( 0.44 ) ) );

        
        /* Inherit variables from WorkspaceObject, put them last */
        const MemberRules& inheritedRules = WorkspaceObject::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() ); 
        
        rules_set = true;
    }
    
    return move_member_rules;
}


/**
 * Print the value for the user.
 */
void Move::printValue(std::ostream &o, bool user) const
{
    printValue(o);
}

/**
 * Print the value for the user.
 */
void Move::printValue(std::ostream &o) const
{
    
    o << getMoveName() << "(";
    o << "?";
    o << ")";
    
}

/** Set a member variable */
void Move::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "weight" ) 
    {
        weight = var;
    }
    else if ( name == "tuneTarget" )
    {
        tune_target = var;
    }
    else if ( name == "delay" )
    {
        delay = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}


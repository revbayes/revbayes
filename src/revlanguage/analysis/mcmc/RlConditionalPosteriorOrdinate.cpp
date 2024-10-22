#include "RlConditionalPosteriorOrdinate.h"

#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Delimiter.h"
#include "RevObject.h"
#include "Real.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "ConditionalPosteriorOrdinate.h"
#include "WorkspaceToCoreWrapperObject.h"
#include "WorkspaceVector.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

ConditionalPosteriorOrdinate::ConditionalPosteriorOrdinate() : WorkspaceToCoreWrapperObject<RevBayesCore::ConditionalPosteriorOrdinate>()
{

    ArgumentRules* cpo_arg_rules = new ArgumentRules();
    cpo_arg_rules->push_back( new ArgumentRule("counts", ModelVector<RealPos>::getClassTypeSpec(), "The number of observation for each site (if a site pattern was observed multiple times and we compressed the data).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<RealPos>(  ) ) );
    cpo_arg_rules->push_back( new ArgumentRule("log", RlBoolean::getClassTypeSpec(), "Whether the probabilities were stored as log values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
    cpo_arg_rules->push_back( new ArgumentRule("folded", RlBoolean::getClassTypeSpec(), "Whether the actual observed data are folded.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
    methods.addFunction(new MemberProcedure( "predictiveProbability", Real::getClassTypeSpec(), cpo_arg_rules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
ConditionalPosteriorOrdinate* ConditionalPosteriorOrdinate::clone(void) const
{
    
    return new ConditionalPosteriorOrdinate(*this);
}


void ConditionalPosteriorOrdinate::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the parameter values
    const std::string&              fn              = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&              del             = static_cast<const RlString &>( delimiter->getRevObject() ).getValue();
    const std::vector<std::string>& cols_to_skip    = static_cast<const ModelVector<RlString> &>( col_names_to_skip->getRevObject() ).getValue();

    value = new RevBayesCore::ConditionalPosteriorOrdinate(fn, del, cols_to_skip);
    
}


/* Map calls to member methods */
RevPtr<RevVariable> ConditionalPosteriorOrdinate::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "predictiveProbability")
    {
        found = true;
        
        const std::vector<double>& counts   = static_cast<const ModelVector<RealPos> &>( args[0].getVariable()->getRevObject() ).getValue();
        bool as_log                         = static_cast<const RlBoolean            &>( args[1].getVariable()->getRevObject() ).getValue();
        bool use_folded                     = static_cast<const RlBoolean            &>( args[2].getVariable()->getRevObject() ).getValue();

        double p = value->predictiveProbability( counts, as_log, use_folded );
        
        return new RevVariable( new Real( p ) );
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& ConditionalPosteriorOrdinate::getClassType(void)
{
    
    static std::string rev_type = "ConditionalPosteriorOrdinate";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& ConditionalPosteriorOrdinate::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::ConditionalPosteriorOrdinate>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string ConditionalPosteriorOrdinate::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ConditionalPosteriorOrdinate";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& ConditionalPosteriorOrdinate::getParameterRules(void) const
{
    
    static MemberRules cpo_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        cpo_member_rules.push_back( new ArgumentRule("filename"            , RlString::getClassTypeSpec(), "The name of the file where the likelhood samples are stored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        RevBayesCore::RbVector<std::string> default_columns_names_to_skip;
        default_columns_names_to_skip.push_back( "Iteration" );
        default_columns_names_to_skip.push_back( "Posterior" );
        default_columns_names_to_skip.push_back( "Likelihood" );
        default_columns_names_to_skip.push_back( "Prior" );
        default_columns_names_to_skip.push_back( "Replicate_ID" );
        cpo_member_rules.push_back( new ArgumentRule("columnNamesToSkip"   , ModelVector<RlString>::getClassTypeSpec(), "The names of the columns that we are going to skip.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<RlString>( default_columns_names_to_skip ) ) );
        cpo_member_rules.push_back( new Delimiter() );

        rules_set = true;
    }
    
    return cpo_member_rules;
}

/** Get type spec */
const TypeSpec& ConditionalPosteriorOrdinate::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void ConditionalPosteriorOrdinate::printValue(std::ostream &o) const
{
    
    o << "ConditionalPosteriorOrdinate";
}


/** Set a member variable */
void ConditionalPosteriorOrdinate::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "filename")
    {
        filename = var;
    }
    else if ( name == "columnNamesToSkip" )
    {
        col_names_to_skip = var;
    }
    else if ( name == "delimiter" || name == "separator" )
    {
        delimiter = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}

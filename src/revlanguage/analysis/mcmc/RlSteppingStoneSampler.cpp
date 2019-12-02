#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "RevObject.h"
#include "RbException.h"
#include "Real.h"
#include "RevNullObject.h"
#include "RlSteppingStoneSampler.h"
#include "RlString.h"
#include "TypeSpec.h"


using namespace RevLanguage;

SteppingStoneSampler::SteppingStoneSampler() : WorkspaceToCoreWrapperObject<RevBayesCore::SteppingStoneSampler>()
{

    ArgumentRules* marginalArgRules = new ArgumentRules();
    methods.addFunction(new MemberProcedure( "marginal", Real::getClassTypeSpec(), marginalArgRules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
SteppingStoneSampler* SteppingStoneSampler::clone(void) const
{
    
	return new SteppingStoneSampler(*this);
}


void SteppingStoneSampler::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the parameter values
    const std::string&   fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&   pn      = static_cast<const RlString &>( powerColumnName->getRevObject() ).getValue();
    const std::string&   ln      = static_cast<const RlString &>( likelihoodColumnName->getRevObject() ).getValue();
    const std::string&   del     = static_cast<const RlString &>( delimmiter->getRevObject() ).getValue();
    
    value = new RevBayesCore::SteppingStoneSampler(fn, pn, ln, del);
    
}


/* Map calls to member methods */
RevPtr<RevVariable> SteppingStoneSampler::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "marginal")
    {
        found = true;
        
        double ml = value->marginalLikelihood();
        
        return new RevVariable( new Real( ml ) );
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& SteppingStoneSampler::getClassType(void)
{
    
    static std::string rev_type = "SteppingStoneSampler";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& SteppingStoneSampler::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::SteppingStoneSampler>::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string SteppingStoneSampler::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "steppingStoneSampler";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& SteppingStoneSampler::getParameterRules(void) const
{
    
    static MemberRules samplerMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        samplerMemberRules.push_back( new ArgumentRule("filename"            , RlString::getClassTypeSpec(), "The name of the file where the likelhood samples are stored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        samplerMemberRules.push_back( new ArgumentRule("powerColumnName"     , RlString::getClassTypeSpec(), "The name of the column of the powers.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        samplerMemberRules.push_back( new ArgumentRule("likelihoodColumnName", RlString::getClassTypeSpec(), "The name of the column of the likelihood samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        samplerMemberRules.push_back( new ArgumentRule("delimiter"           , RlString::getClassTypeSpec(), "The column delimiter.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "\t" ) ) );
        
        rules_set = true;
    }
    
    return samplerMemberRules;
}

/** Get type spec */
const TypeSpec& SteppingStoneSampler::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void SteppingStoneSampler::printValue(std::ostream &o) const
{
    
    o << "SteppingStoneSampler";
}


/** Set a member variable */
void SteppingStoneSampler::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "likelihoodColumnName")
    {
        likelihoodColumnName = var;
    }
    else if ( name == "powerColumnName")
    {
        powerColumnName = var;
    }
    else if ( name == "filename")
    {
        filename = var;
    }
    else if ( name == "delimiter")
    {
        delimmiter = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}

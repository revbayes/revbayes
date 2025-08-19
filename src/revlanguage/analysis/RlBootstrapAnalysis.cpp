#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "BootstrapAnalysis.h"
#include "RealPos.h"
#include "RlMaximumLikelihoodAnalysis.h"
#include "RlBootstrapAnalysis.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevBayesCore { class MaximumLikelihoodAnalysis; }


using namespace RevLanguage;

BootstrapAnalysis::BootstrapAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::BootstrapAnalysis>()
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule( "epsilon", RealPos::getClassTypeSpec(), "The minimum improvement in the last interval.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.001) ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
        
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
BootstrapAnalysis* BootstrapAnalysis::clone(void) const
{
    
    return new BootstrapAnalysis(*this);
}


void BootstrapAnalysis::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const RevBayesCore::MaximumLikelihoodAnalysis&  s   = static_cast<const MaximumLikelihoodAnalysis &>( sampler->getRevObject() ).getValue();
    int                                             n   = static_cast<const Natural &>( simulations->getRevObject() ).getValue();
    
    value = new RevBayesCore::BootstrapAnalysis( s, size_t(n) );
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string BootstrapAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "BootstrapAnalysis";
    
    return c_name;
}


/* Map calls to member methods */
RevPtr<RevVariable> BootstrapAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "run")
    {
        
        found = true;
        
        // get the member with give index
        double e = static_cast<const RealPos &>( args[0].getVariable()->getRevObject() ).getValue();
        value->runAll( e );
        
        return NULL;
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& BootstrapAnalysis::getClassType(void)
{
    
    static std::string rev_type = "BootstrapAnalysis";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& BootstrapAnalysis::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::BootstrapAnalysis>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules (no members) */
const MemberRules& BootstrapAnalysis::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("estimator"     , MaximumLikelihoodAnalysis::getClassTypeSpec(), "The maximum likelihood estimation object.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("simulations"   , Natural::getClassTypeSpec()                  , "How many replicate simulations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& BootstrapAnalysis::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void BootstrapAnalysis::printValue(std::ostream &o) const
{
    
    o << "BootstrapAnalysis";
}


/** Set a member variable */
void BootstrapAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "estimator")
    {
        sampler = var;
    }
    else if ( name == "simulations")
    {
        simulations = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}

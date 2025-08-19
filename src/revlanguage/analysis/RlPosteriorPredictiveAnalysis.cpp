#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "PosteriorPredictiveAnalysis.h"
#include "RlMonteCarloAnalysis.h"
#include "RlPosteriorPredictiveAnalysis.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevBayesCore { class MonteCarloAnalysis; }


using namespace RevLanguage;

PosteriorPredictiveAnalysis::PosteriorPredictiveAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::PosteriorPredictiveAnalysis>()
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule("generations", Natural::getClassTypeSpec(), "The number of generation to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
    
//    ArgumentRules* burninArgRules = new ArgumentRules();
//    burninArgRules->push_back( new ArgumentRule("generations"   , Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//    burninArgRules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The number of iterations after which we tune the parameters of the moves.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//    methods.addFunction( new MemberProcedure( "burnin", RlUtils::Void, burninArgRules) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
PosteriorPredictiveAnalysis* PosteriorPredictiveAnalysis::clone(void) const
{
    
    return new PosteriorPredictiveAnalysis(*this);
}


void PosteriorPredictiveAnalysis::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const RevBayesCore::MonteCarloAnalysis&         s     = static_cast<const MonteCarloAnalysis &>( sampler->getRevObject() ).getValue();
    const std::string&                              dir   = static_cast<const RlString &>( directory->getRevObject() ).getValue();
    
    value = new RevBayesCore::PosteriorPredictiveAnalysis( s, dir );
    
}


/* Map calls to member methods */
RevPtr<RevVariable> PosteriorPredictiveAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "run")
    {
        
        found = true;
        
        // get the member with give index
        int gen = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        value->runAll( size_t(gen) );
        
        return NULL;
    }
//    else if (name == "burnin")
//    {
//        found = true;
//        
//        // get the member with give index
//        int gen = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
//        int tuningInterval = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
//        value->burnin( size_t(gen), size_t(tuningInterval) );
//        
//        return NULL;
//    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& PosteriorPredictiveAnalysis::getClassType(void)
{
    
    static std::string rev_type = "PosteriorPredictiveAnalysis";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& PosteriorPredictiveAnalysis::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::PosteriorPredictiveAnalysis>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string PosteriorPredictiveAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "posteriorPredictiveAnalysis";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& PosteriorPredictiveAnalysis::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("sampler"    , MonteCarloAnalysis::getClassTypeSpec(), "The template Monte Carlo sampler instance.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("directory"  , RlString::getClassTypeSpec()          , "The name of the directory where the simulated data are.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& PosteriorPredictiveAnalysis::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void PosteriorPredictiveAnalysis::printValue(std::ostream &o) const
{
    
    o << "PowerPosterior";
}


/** Set a member variable */
void PosteriorPredictiveAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "sampler")
    {
        sampler = var;
    }
    else if ( name == "directory")
    {
        directory = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}

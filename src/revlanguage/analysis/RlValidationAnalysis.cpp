#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "Probability.h"
#include "ValidationAnalysis.h"
#include "RlMonteCarloAnalysis.h"
#include "RlValidationAnalysis.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlString.h"
#include "RlUtils.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevBayesCore { class MonteCarloAnalysis; }


using namespace RevLanguage;

ValidationAnalysis::ValidationAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::ValidationAnalysis>()
{
    initializeMethods();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the analysis.
 */
ValidationAnalysis* ValidationAnalysis::clone(void) const
{
    
    return new ValidationAnalysis(*this);
}

/** Construct a new internal object and sets value to it **/
void ValidationAnalysis::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new validation analysis
    const RevBayesCore::MonteCarloAnalysis&         s   = static_cast<const MonteCarloAnalysis &>( sampler->getRevObject() ).getValue();
    int                                             n   = (int)static_cast<const Natural &>( simulations->getRevObject() ).getValue();
    const std::string&                              d   = static_cast<const RlString &>( output_directory->getRevObject() ).getValue();

    value = new RevBayesCore::ValidationAnalysis( s, size_t(n), d );
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string ValidationAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "validationAnalysis";
    
    return c_name;
}


/** Map calls to member methods
* @return result of the call
**/
RevPtr<RevVariable> ValidationAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "run")
    {
        
        found = true;
        
        // get the member with given index
        
        int gen = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        value->runAll( size_t(gen) );
        
        return NULL;
    }
    else if (name == "burnin")
    {
        found = true;
        
        // get the member with given index
        int gen = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        int tuningInterval = (int)static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
        value->burnin( size_t(gen), size_t(tuningInterval) );
        
        return NULL;
    }
    else if (name == "summarize")
    {
        found = true;
        
        double coverage = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        
        value->summarizeAll(coverage);
        
        return NULL;
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& ValidationAnalysis::getClassType(void)
{
    
    static std::string rev_type = "ValidationAnalysis";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& ValidationAnalysis::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::ValidationAnalysis>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules */
const MemberRules& ValidationAnalysis::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("sampler"       , MonteCarloAnalysis::getClassTypeSpec(), "The template Monte Carlo sampler instance.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("simulations"   , Natural::getClassTypeSpec()           , "How many replicate simulations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("directory"     , RlString::getClassTypeSpec()          , "The directory for the output of the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "output" ) ) );

        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& ValidationAnalysis::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Initialize the member methods */
void ValidationAnalysis::initializeMethods()
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule( "generations", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
    
    ArgumentRules* burninArgRules = new ArgumentRules();
    burninArgRules->push_back( new ArgumentRule( "generations"   , Natural::getClassTypeSpec(), "The number of generation to run this burnin simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    burninArgRules->push_back( new ArgumentRule( "tuningInterval", Natural::getClassTypeSpec(), "The interval when to update the tuning parameters of the moves.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    methods.addFunction( new MemberProcedure( "burnin", RlUtils::Void, burninArgRules) );

    ArgumentRules* summarizeArgRules = new ArgumentRules();
    summarizeArgRules->push_back( new ArgumentRule( "coverageProbability", Probability::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.9) ) );
    methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarizeArgRules) );
    
}


void ValidationAnalysis::printValue(std::ostream &o) const
{
    
    o << "ValidationAnalysis";
}


/** Set member variable */
void ValidationAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "sampler")
    {
        sampler = var;
    }
    else if ( name == "simulations")
    {
        simulations = var;
    }
    else if ( name == "directory")
    {
        output_directory = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}

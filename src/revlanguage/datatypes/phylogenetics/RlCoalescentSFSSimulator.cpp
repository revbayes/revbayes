#include "CoalescentSFSSimulator.h"

#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Delimiter.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "Real.h"
#include "RlCoalescentSFSSimulator.h"
#include "RlDemographicFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlMatrixRealPos.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

CoalescentSFSSimulator::CoalescentSFSSimulator() : WorkspaceToCoreWrapperObject<RevBayesCore::CoalescentSFSSimulator>()
{
    
    // simulating an SFS
    
    ArgumentRules* sfs_arg_rules = new ArgumentRules();
    
    sfs_arg_rules->push_back( new ArgumentRule( "mutationRate"  , RealPos::getClassTypeSpec(), "The mutation rate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    sfs_arg_rules->push_back( new ArgumentRule( "sampleSize"    , Natural::getClassTypeSpec(), "The sample size, i.e., number of individuals at present.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    sfs_arg_rules->push_back( new ArgumentRule( "reps"          , Natural::getClassTypeSpec(), "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction(new MemberProcedure( "simulateSFS", ModelVector< Natural >::getClassTypeSpec(), sfs_arg_rules) );
    
    // simulating coalescent trees and storing statistics in a file
    
    ArgumentRules* coal_arg_rules = new ArgumentRules();
    
    coal_arg_rules->push_back( new ArgumentRule( "sampleSize"    , Natural::getClassTypeSpec(),  "The sample size, i.e., number of individuals at present.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    coal_arg_rules->push_back( new ArgumentRule( "reps"          , Natural::getClassTypeSpec(),  "The number of replicate to simulate the frequencies.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    coal_arg_rules->push_back( new ArgumentRule( "statsFilename" , RlString::getClassTypeSpec(), "The filename in which to store the coalescent summary statistics (empty if you don't want it).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
    coal_arg_rules->push_back( new ArgumentRule( "treesFilename" , RlString::getClassTypeSpec(), "The filename in which to store the coalescent trees (empty if you don't want it).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );

    methods.addFunction(new MemberProcedure( "simulateCoalescent", ModelVector< Natural >::getClassTypeSpec(), coal_arg_rules) );

    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
CoalescentSFSSimulator* CoalescentSFSSimulator::clone(void) const
{
    
    return new CoalescentSFSSimulator(*this);
}


void CoalescentSFSSimulator::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the parameter values
    std::vector<double> cp  = static_cast<const ModelVector<RealPos> &>( change_points->getRevObject() ).getValue();
    double gt               = static_cast<const RealPos &>( generation_time->getRevObject() ).getValue();
    const std::string& p    = static_cast<const RlString &>( ploidy->getRevObject() ).getValue();
    
    // demographic functions
    const WorkspaceVector<DemographicFunction> &ws_vec_df   = static_cast<const WorkspaceVector<DemographicFunction> &>( demographies->getRevObject() );
    RevBayesCore::RbVector<RevBayesCore::DemographicFunction> df;
    for ( size_t i = 0; i < ws_vec_df.size(); ++i )
    {
        df.push_back( ws_vec_df[i].getValue() );
    }
    
    value = new RevBayesCore::CoalescentSFSSimulator(df, cp, gt, p);

    
}


/* Map calls to member methods */
RevPtr<RevVariable> CoalescentSFSSimulator::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "simulateSFS")
    {
        found = true;
        
        size_t args_index   = 0;
        
        double mu           = static_cast<const RealPos &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        long sample_size    = static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        long reps           = static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();

        RevBayesCore::RbVector<long>* m = value->simulateSFS(mu, sample_size, reps);
        
        return new RevVariable( new ModelVector<Natural>( *m ) );
    }
    else if (name == "simulateCoalescent")
    {
        found = true;
        
        size_t args_index   = 0;
            
        long                sample_size     = static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        long                reps            = static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        const std::string&  f_name_stats    = static_cast<const RlString &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        const std::string&  f_name_trees    = static_cast<const RlString &>( args[args_index++].getVariable()->getRevObject() ).getValue();

        value->simulateCoalescent(sample_size, reps, f_name_stats, f_name_trees);
            
        return NULL;
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::CoalescentSFSSimulator>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& CoalescentSFSSimulator::getClassType(void)
{
    
    static std::string rev_type = "CoalescentSFSSimulator";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& CoalescentSFSSimulator::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::CoalescentSFSSimulator>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string CoalescentSFSSimulator::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CoalescentSFSSimulator";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& CoalescentSFSSimulator::getParameterRules(void) const
{
    
    static MemberRules argument_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "demographies"      , WorkspaceVector<DemographicFunction>::getClassTypeSpec(), "The vector of demographic functions.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "changePoints"      , ModelVector<RealPos>::getClassTypeSpec(), "The start times of the intervals (the first interval is implicit and starts at 0).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "generationTime"    , RealPos::getClassTypeSpec(), "The generation time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> ploidy_options;
        ploidy_options.push_back( "haploid" );
        ploidy_options.push_back( "diploid" );
        argument_rules.push_back( new OptionRule( "ploidy", new RlString("diploid"), ploidy_options, "The ploidy type of the individuals in the population, to scale N if necessary." ) );
        
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get type spec */
const TypeSpec& CoalescentSFSSimulator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void CoalescentSFSSimulator::printValue(std::ostream &o) const
{
    
    o << "CoalescentSFSSimulator";
}


/** Set a member variable */
void CoalescentSFSSimulator::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "generationTime")
    {
        generation_time = var;
    }
    else if ( name == "demographies")
    {
        demographies = var;
    }
    else if ( name == "changePoints")
    {
        change_points = var;
    }
    else if ( name == "ploidy")
    {
        ploidy = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}

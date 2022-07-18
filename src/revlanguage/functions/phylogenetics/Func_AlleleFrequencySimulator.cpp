#include "Func_AlleleFrequencySimulator.h"

#include <map>
#include <string>
#include <vector>

#include "AlleleFrequencySimulator.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTaxon.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "HomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/** default constructor */
Func_AlleleFrequencySimulator::Func_AlleleFrequencySimulator( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_AlleleFrequencySimulator* Func_AlleleFrequencySimulator::clone( void ) const
{
    
    return new Func_AlleleFrequencySimulator( *this );
}


RevPtr<RevVariable> Func_AlleleFrequencySimulator::execute()
{
    
    RevBayesCore::Tree* tree = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getValue().clone();
    
    std::vector<long> population_sizes      = static_cast<const ModelVector<Natural> &>( this->args[1].getVariable()->getRevObject() ).getValue();
    double generation_time                  = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getValue();
    bool moran_generations                  = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getValue();
    long num_sites                          = static_cast<const Natural &>( this->args[4].getVariable()->getRevObject() ).getValue();
    std::vector<double> mutation_rates      = static_cast<const ModelVector<RealPos> &>( this->args[5].getVariable()->getRevObject() ).getValue();
    std::vector<long> samples_per_species   = static_cast<const ModelVector<Natural> &>( this->args[6].getVariable()->getRevObject() ).getValue();
    double root_branch                      = static_cast<const RealPos &>( this->args[7].getVariable()->getRevObject() ).getValue();
    bool variable                           = static_cast<const RlBoolean &>( this->args[8].getVariable()->getRevObject() ).getValue();
    const std::string& fn                   = static_cast<const RlString &>( this->args[9].getVariable()->getRevObject() ).getValue();


    RevBayesCore::AlleleFrequencySimulator sim = RevBayesCore::AlleleFrequencySimulator( generation_time, mutation_rates, moran_generations );
    sim.simulateAlleleFrequencies(tree, population_sizes, num_sites, samples_per_species, root_branch, fn, variable);
    
    return NULL;
}


/* Get argument rules */
const ArgumentRules& Func_AlleleFrequencySimulator::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "tree"              , TimeTree::getClassTypeSpec(), "The tree along which we want to simulate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "populationSizes"   , ModelVector<Natural>::getClassTypeSpec(), "The population sizes for all branches, including the root branch.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "generationTime"    , RealPos::getClassTypeSpec(), "The generation time for the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "moranGenerations"  , RlBoolean::getClassTypeSpec(), "Is the generation time in Moran generation, i.e., scaled by N.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "numSites"          , Natural::getClassTypeSpec(), "The number of sites to simulate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mutationRates"     , ModelVector<RealPos>::getClassTypeSpec(), "The mutation rates from 0 to 1 and 1 to 0.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "samplesPerSpecies" , ModelVector<Natural>::getClassTypeSpec(), "The observed number of samples per species.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "rootBranch"        , RealPos::getClassTypeSpec(), "The length of the root branch for simulating polymorphisms at the root.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "variable"          , RlBoolean::getClassTypeSpec(), "Do we condition on only observing variable sites?", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "filename"          , RlString::getClassTypeSpec(), "The filename for the counts file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


const std::string& Func_AlleleFrequencySimulator::getClassType(void)
{
    
    static std::string rev_type = "Func_AlleleFrequencySimulator";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_AlleleFrequencySimulator::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_AlleleFrequencySimulator::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "AlleleFrequencySimulator";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_AlleleFrequencySimulator::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RevNullObject::getClassTypeSpec();
    return return_typeSpec;
}


const TypeSpec& Func_AlleleFrequencySimulator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

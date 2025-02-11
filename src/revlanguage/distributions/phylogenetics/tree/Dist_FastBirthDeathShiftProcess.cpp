#include <math.h>
#include <cstddef>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_FastBirthDeathShiftProcess.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDistributionMemberFunction.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "FastBirthDeathShiftProcess.h"
#include "StochasticNode.h"
#include "Integer.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "RlTypedDistribution.h"
#include "RlUtils.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;


Dist_FastBirthDeathShiftProcess::Dist_FastBirthDeathShiftProcess() : TypedDistribution<TimeTree>()
{
    
}


Dist_FastBirthDeathShiftProcess::~Dist_FastBirthDeathShiftProcess()
{
    
}



Dist_FastBirthDeathShiftProcess* Dist_FastBirthDeathShiftProcess::clone( void ) const
{
    
    return new Dist_FastBirthDeathShiftProcess( *this );
}


RevBayesCore::TypedDistribution<RevBayesCore::Tree>* Dist_FastBirthDeathShiftProcess::createDistribution( void ) const
{
    
    // get process start age
    RevBayesCore::TypedDagNode<double>* ra   = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();
    
    // condition off origin age or root age?
    bool uo = ( start_condition == "originAge" ? true : false );
   
    // get speciation scale 
    RevBayesCore::TypedDagNode<double>* sp_scale = nullptr;
    sp_scale = static_cast<const RealPos &>( speciation_scale->getRevObject() ).getDagNode();

    // extinction scale
    RevBayesCore::TypedDagNode<double>* ex_scale = nullptr;
    ex_scale = static_cast<const RealPos &>( extinction_scale->getRevObject() ).getDagNode();

    // get spread parameters
    RevBayesCore::TypedDagNode<double>* sp_sd = nullptr;
    sp_sd = static_cast<const RealPos &>( speciation_sd->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode<double>* ex_sd = nullptr;
    ex_sd = static_cast<const RealPos &>( extinction_sd->getRevObject() ).getDagNode();

    // get alpha and beta
    RevBayesCore::TypedDagNode<double>* r_sp = static_cast<const RealPos &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* r_ext = static_cast<const RealPos &>( beta->getRevObject() ).getDagNode();

    // get number of rate classes
    //size_t num_classes = static_cast<const Integer &>( num_rate_classes->getRevObject() ).getValue();
    size_t num_speciation_classes = static_cast<const Integer &>( num_sp_classes->getRevObject() ).getValue();
    size_t num_extinction_classes = static_cast<const Integer &>( num_ex_classes->getRevObject() ).getValue();
        
    // condition
    const std::string& cond                  = static_cast<const RlString &>( condition->getRevObject() ).getValue();
   
    // condition for simulating under
    const std::string& simulate_cond         = static_cast<const RlString &>( simulation_condition->getRevObject() ).getValue();
    
    bool cond_tip_states = false;
    bool cond_num_tips = false;
    bool cond_tree = false;
    if (simulate_cond == "tipStates")
    {
        cond_tip_states = true;
    }
    else if (simulate_cond == "numTips")
    {
        cond_num_tips = true;
    }
    else if (simulate_cond == "tree")
    {
        cond_tree = true;
    }
    
    size_t max_l = static_cast<const Integer &>( max_lineages->getRevObject() ).getValue();
    size_t min_l = static_cast<const Integer &>( min_lineages->getRevObject() ).getValue();
    size_t exact_l = static_cast<const Integer &>( exact_lineages->getRevObject() ).getValue();
    double max_t = static_cast<const RealPos &>( max_time->getRevObject() ).getValue();
    
    size_t prune = static_cast<const RlBoolean &>( prune_extinct_lineages->getRevObject() ).getValue();
    
    // finally make the distribution 
    RevBayesCore::FastBirthDeathShiftProcess*   d = new RevBayesCore::FastBirthDeathShiftProcess( 
            ra, 
            sp_scale, 
            ex_scale, 
            sp_sd, 
            ex_sd, 
            r_sp, 
            r_ext, 
            num_speciation_classes, 
            num_extinction_classes, 
            cond, 
            uo, 
            min_l,
            max_l,
            exact_l,
            max_t,
            prune, 
            cond_tip_states,
            cond_num_tips,
            cond_tree);
   
    

    // set sampling probabilities/fractions
    RevBayesCore::TypedDagNode<double>* rh = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    d->setSamplingFraction( rh );
    
   
    return d;
}



/* Get Rev type of object */
const std::string& Dist_FastBirthDeathShiftProcess::getClassType(void)
{
    
    static std::string rev_type = "Dist_FastBirthDeathShiftProcess";
    
    return rev_type;
}


/* Get class type spec describing type of object. TODO: Check if the correct parent is TypedDistribution or Distribution */
const TypeSpec& Dist_FastBirthDeathShiftProcess::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_FastBirthDeathShiftProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "FastBirthDeathShift";
    
    return d_name;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_FastBirthDeathShiftProcess::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "dist_alias_here" );

    return a_names;
}


MethodTable Dist_FastBirthDeathShiftProcess::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();
    
    ArgumentRules* clampCharDataArgRules = new ArgumentRules();
    clampCharDataArgRules->push_back( new ArgumentRule( "value", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "clampCharData", RlUtils::Void, clampCharDataArgRules ) );
   
    ArgumentRules* getCharDataArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getCharData", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), getCharDataArgRules ) ); 
    
    ArgumentRules* avgSpeciationArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FastBirthDeathShiftProcess, ModelVector<RealPos> >( "averageSpeciationRate", variable, avgSpeciationArgRules   ) );

    ArgumentRules* avgExtinctionArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FastBirthDeathShiftProcess, ModelVector<RealPos> >( "averageExtinctionRate", variable, avgExtinctionArgRules   ) );
    
    ArgumentRules* num_events_arg_rules = new ArgumentRules();
    //    parentArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_FastBirthDeathShiftProcess, ModelVector<Natural> >( "numberEvents", variable, num_events_arg_rules   ) );
    
    ArgumentRules* timeInStateArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FastBirthDeathShiftProcess, ModelVector<RealPos> >( "getTimeInState", variable, timeInStateArgRules   ) );
    
    return methods;
}


/* Return member rules */
const MemberRules& Dist_FastBirthDeathShiftProcess::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        memberRules.push_back( new ArgumentRule( aliases, RealPos::getClassTypeSpec()    , "The start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );


        std::vector<std::string> slabels;
        slabels.push_back("speciationScale");
        std::vector<TypeSpec> speciationTypes;
        speciationTypes.push_back( RealPos::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule( slabels     , speciationTypes ,                          "The scale (in the log-normal distribution) of the speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );

        std::vector<std::string> elabels;
        elabels.push_back("extinctionScale");
        memberRules.push_back( new ArgumentRule( elabels     , RealPos::getClassTypeSpec() , "The scale (in the log-normal distribution) of the speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );

        memberRules.push_back( new ArgumentRule( "speciationSD", RealPos::getClassTypeSpec(), "The spread (sigma-parameter in the log-normal distribution) of the speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.587)));
        memberRules.push_back( new ArgumentRule( "extinctionSD", RealPos::getClassTypeSpec(), "The spread (sigma-parameter in the log-normal distribution) of the extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.587)));

        memberRules.push_back( new ArgumentRule( "alpha", RealPos::getClassTypeSpec()        , "the rate of shifts in speciation rate", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "beta", RealPos::getClassTypeSpec()        , "the rate of shifts in extinction rate", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "numSpeciationClasses", Natural::getClassTypeSpec(), "The number of speciation rate classes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(6) ) );
        memberRules.push_back( new ArgumentRule( "numExtinctionClasses", Natural::getClassTypeSpec(), "The number of extinction rate classes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(6) ) );
        
        std::vector<TypeSpec> sampling_fraction_types;
        sampling_fraction_types.push_back( Probability::getClassTypeSpec() );
        sampling_fraction_types.push_back( ModelVector<Probability>::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule( "rho"       , sampling_fraction_types, "The taxon sampling probability."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        memberRules.push_back( new OptionRule( "condition"    , new RlString("time"), optionsCondition, "The condition of the birth-death process." ) );
        std::vector<std::string> optionsSimulateCondition;
        optionsSimulateCondition.push_back("startTime");
        optionsSimulateCondition.push_back("numTips");
        optionsSimulateCondition.push_back("tipStates");
        optionsSimulateCondition.push_back("tree");
        memberRules.push_back( new OptionRule("simulateCondition", new RlString("startTime"), optionsSimulateCondition, "The conditions under which to simulate." ) );
        memberRules.push_back( new ArgumentRule("minNumLineages", Natural::getClassTypeSpec(), "The minimum number of lineages to simulate; applied under the startTime condition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural() ) );
        memberRules.push_back( new ArgumentRule("maxNumLineages", Natural::getClassTypeSpec(), "The maximum number of lineages to simulate; applied under the startTime condition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(500) ) );
        memberRules.push_back( new ArgumentRule("exactNumLineages", Natural::getClassTypeSpec(), "The exact number of lineages to simulate; applied under the numTips condition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
        memberRules.push_back( new ArgumentRule("maxTime", RealPos::getClassTypeSpec(), "Maximum time for lineages to coalesce when simulating; applied under the numTips and tipStates condition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1000.0) ) );
        memberRules.push_back( new ArgumentRule("pruneExtinctLineages", RlBoolean::getClassTypeSpec(), "When simulating should extinct lineages be pruned off?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_FastBirthDeathShiftProcess::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_FastBirthDeathShiftProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "rootAge" || name == "originAge" )
    {
        start_age = var;
        start_condition = name;
    }
    else if ( name == "speciationScale" ) 
    {
        speciation_scale = var;
    }
    else if ( name == "extinctionScale" )
    {
        extinction_scale = var;
    }
    else if ( name == "speciationSD" )
    {
        speciation_sd = var;
    }
    else if ( name == "extinctionSD" )
    {
        extinction_sd = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "beta" )
    {
        beta = var;
    }
    else if ( name == "numSpeciationClasses" )
    {
        num_sp_classes = var;
    }
    else if ( name == "numExtinctionClasses" )
    {
        num_ex_classes = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "minNumLineages" )
    {
        min_lineages = var;
    }
    else if ( name == "maxNumLineages" )
    {
        max_lineages = var;
    }
    else if ( name == "exactNumLineages" )
    {
        exact_lineages = var;
    }
    else if ( name == "maxTime" )
    {
        max_time = var;
    }
    else if ( name == "pruneExtinctLineages" )
    {
        prune_extinct_lineages = var;
    }
    else if ( name == "simulateCondition" )
    {
        simulation_condition = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}


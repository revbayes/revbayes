#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_FastBirthDeathShiftProcess.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RateGenerator.h"
#include "RealPos.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDistributionMemberFunction.h"
#include "RlRateGenerator.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "FastBirthDeathShiftProcess.h"
#include "StochasticNode.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "Integer.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "RlUtils.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

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
   
    // extinction rates
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ex  = static_cast<const ModelVector<RealPos> &>( extinction_rates->getRevObject() ).getDagNode();
    
    // get speciation rates or cladogenetic speciation rate event map
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sp = nullptr;
    sp  = static_cast<const ModelVector<RealPos> &>( speciation_rates->getRevObject() ).getDagNode();
        
    // rate matrix
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* q      = NULL;
    if ( event_rate_matrix->getRevObject() != RevNullObject::getInstance() )
    {
        q = static_cast<const RateGenerator &>( event_rate_matrix->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<double>* r = NULL;
    if ( event_rate->getRevObject() != RevNullObject::getInstance() )
    {
        r = static_cast<const RealPos &>( event_rate->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<RevBayesCore::Simplex >* bf = NULL;
    if ( root_frequencies->getRevObject() != RevNullObject::getInstance() )
    {
        bf = static_cast<const Simplex &>( root_frequencies->getRevObject() ).getDagNode();
    }
    
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
    
    bool allow_shifts_extinct = static_cast<const RlBoolean &>( allow->getRevObject() ).getValue();
    
    // finally make the distribution 
    RevBayesCore::FastBirthDeathShiftProcess*   d = new RevBayesCore::FastBirthDeathShiftProcess( ra, ex, q, r, bf, cond, uo, min_l, max_l, exact_l, max_t, prune, cond_tip_states, cond_num_tips, cond_tree, allow_shifts_extinct );
   
    
    size_t ex_size = ex->getValue().size();
    size_t q_size = q->getValue().size();
    size_t sp_size = 0;
    
    // set speciation/cladogenetic event rates
    d->setSpeciationRates( sp );
    sp_size = sp->getValue().size();

    std::stringstream ss_err;
    if (ex_size != q_size)
    {
        ss_err << "State count mismatch between extinction rates (" << ex_size << ") and Q (" << q_size << ")";
        throw RbException(ss_err.str());
    }
    if (ex_size != sp_size)
    {
        ss_err << "State count mismatch between extinction rates (" << ex_size << ") and speciation rates (" << sp_size << ")";
        throw RbException(ss_err.str());
    }
    if (q_size != sp_size)
    {
        ss_err << "State count mismatch between speciation rates (" << sp_size << ") and Q (" << q_size << ")";
        throw RbException(ss_err.str());
    }
    

    // set sampling probabilities/fractions
    if (rho->getRevObject().isType( Probability::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode<double>* rh   = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
        d->setSamplingFraction( rh );
    }
    else if (rho->getRevObject().isType( ModelVector<Probability>::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* rh   = static_cast<const ModelVector<Probability> &>( rho->getRevObject() ).getDagNode();
        
        if (rh->getValue().size() != ex_size)
        {
            ss_err << "State count mismatch between extinction rates (" << ex_size << ") and sampling probabilities (" << rh->getValue().size() << ")";
            throw RbException(ss_err.str());
        }
        
        d->setSamplingFraction( rh );
    }
    
    // set the number of time slices for the numeric ODE
    double n = static_cast<const RealPos &>( num_time_slices->getRevObject() ).getValue();
    d->setNumberOfTimeSlices( n );
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ps = NULL;
    if ( psi->getRevObject() != RevNullObject::getInstance() )
    {
        ps  = static_cast<const ModelVector<RealPos> &>( psi->getRevObject() ).getDagNode();

        d->setSerialSamplingRates( ps );
    }

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
    
    ArgumentRules* getCharHistoryArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getCharHistory", RlString::getClassTypeSpec(), getCharHistoryArgRules   ) );
    
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
        slabels.push_back("speciationRates");
        slabels.push_back("lambda");
        std::vector<TypeSpec> speciationTypes;
        speciationTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule( slabels     , speciationTypes ,                          "The vector of speciation rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );

        std::vector<std::string> elabels;
        elabels.push_back("extinctionRates");
        elabels.push_back("mu");
        memberRules.push_back( new ArgumentRule( elabels     , ModelVector<RealPos>::getClassTypeSpec() , "The vector of extinction rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        std::vector<std::string> flabels;
        flabels.push_back("psi");
        flabels.push_back("phi");
        memberRules.push_back( new ArgumentRule( flabels     , ModelVector<RealPos>::getClassTypeSpec() , "The vector of serial sampling rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "Q"         , RateGenerator::getClassTypeSpec()        , "The rate matrix of jumping between rate categories.", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "delta"     , RealPos::getClassTypeSpec()              , "The rate-factor of jumping between rate categories.", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, new RealPos(1.0) ) );
        memberRules.push_back( new ArgumentRule( "pi"        , Simplex::getClassTypeSpec()              , "State frequencies at the root."              , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        std::vector<TypeSpec> sampling_fraction_types;
        sampling_fraction_types.push_back( Probability::getClassTypeSpec() );
        sampling_fraction_types.push_back( ModelVector<Probability>::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule( "rho"       , sampling_fraction_types, "The taxon sampling probability."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        memberRules.push_back( new OptionRule( "condition"    , new RlString("time"), optionsCondition, "The condition of the birth-death process." ) );
        memberRules.push_back( new ArgumentRule("nTimeSlices", RealPos::getClassTypeSpec(), "The number of time slices for the numeric ODE.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(500.0) ) );
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
        memberRules.push_back( new ArgumentRule("allowRateShiftsAtExtinctLineages", RlBoolean::getClassTypeSpec(), "Should we allow rate shifts to occur on extinct lineages?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

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
    else if ( name == "pi" )
    {
        root_frequencies = var;
    }
    else if ( name == "speciationRates" || name == "lambda" || name == "cladoEventMap" )
    {
        speciation_rates = var;
    }
    else if ( name == "extinctionRates" || name == "mu" )
    {
        extinction_rates = var;
    }
    else if ( name == "psi" || name == "phi" )
    {
        psi = var;
    }
    else if ( name == "Q" )
    {
        event_rate_matrix = var;
    }
    else if ( name == "delta" )
    {
        event_rate = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "nTimeSlices" )
    {
        num_time_slices = var;
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
    else if ( name == "allowRateShiftsAtExtinctLineages" )
    {
        allow = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}


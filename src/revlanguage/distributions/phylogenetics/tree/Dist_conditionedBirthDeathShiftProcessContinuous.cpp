#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_conditionedBirthDeathShiftProcessContinuous.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlDistributionMemberFunction.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "ConditionedBirthDeathShiftProcessContinuous.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
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
#include "Taxon.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;


Dist_conditionedBirthDeathShiftProcessContinuous::Dist_conditionedBirthDeathShiftProcessContinuous() : TypedDistribution<TimeTree>()
{
    
}


Dist_conditionedBirthDeathShiftProcessContinuous::~Dist_conditionedBirthDeathShiftProcessContinuous()
{
    
}



Dist_conditionedBirthDeathShiftProcessContinuous* Dist_conditionedBirthDeathShiftProcessContinuous::clone( void ) const
{
    
    return new Dist_conditionedBirthDeathShiftProcessContinuous( *this );
}


RevBayesCore::ConditionedBirthDeathShiftProcessContinuous* Dist_conditionedBirthDeathShiftProcessContinuous::createDistribution( void ) const
{
    
    // Get the parameters
    RevBayesCore::TypedDagNode<double>* ra   = static_cast<const RealPos &>( root_age->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* r_sp = static_cast<const RealPos &>( root_speciation->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* r_ex = static_cast<const RealPos &>( root_extinction->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDistribution<double>* distribution_speciation = NULL;
    if ( speciation->getRevObject() != RevNullObject::getInstance() )
    {
        const Distribution& rl_distribution_speciation  = static_cast<const Distribution &>( speciation->getRevObject() );
        distribution_speciation = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_distribution_speciation.createDistribution() );
    }
    
    RevBayesCore::TypedDistribution<double>* distribution_extinction = NULL;
    if ( extinction->getRevObject() != RevNullObject::getInstance() )
    {
        const Distribution& rl_distribution_extinction  = static_cast<const Distribution &>( extinction->getRevObject() );
        distribution_extinction = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_distribution_extinction.createDistribution() );
    }

    RevBayesCore::TypedDagNode<double>* er   = static_cast<const RealPos &>( event_rate->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rh   = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    // condition
    const std::string& cond                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    RevBayesCore::ConditionedBirthDeathShiftProcessContinuous*   d = new RevBayesCore::ConditionedBirthDeathShiftProcessContinuous( ra, r_sp, r_ex, distribution_speciation, distribution_extinction, er, rh, cond, t );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_conditionedBirthDeathShiftProcessContinuous::getClassType(void)
{
    
    static std::string rev_type = "Dist_condBDS";
    
    return rev_type;
}


/* Get class type spec describing type of object. TODO: Check if the correct parent is TypedDistribution or Distribution */
const TypeSpec& Dist_conditionedBirthDeathShiftProcessContinuous::getClassTypeSpec(void)
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
std::string Dist_conditionedBirthDeathShiftProcessContinuous::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "CBDSP";
    
    return d_name;
}


/**
 * Get the method table for this distribution.
 * We need to implement this function when a random variable drawn from this distribution
 * has specific member method. Here, these are:
 * - x.numberEvents()
 * - x.averageSpeciationRate()
 * - x.averageExtinctionRate()
 * - x.probabilitySpeciationRateIncrease()
 * - x.probabilitySpeciationRateDecrease()
 * - x.probabilityExtinctionRateIncrease()
 * - x.probabilityExtinctionRateDecrease()
 *
 * \return Rev name of constructor function.
 */
MethodTable Dist_conditionedBirthDeathShiftProcessContinuous::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();
    
    // member functions
    ArgumentRules* numEventsArgRules = new ArgumentRules();
    //    parentArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<Natural> >( "numberEvents", variable, numEventsArgRules   ) );
    
    ArgumentRules* avgSpeciationArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "averageSpeciationRate", variable, avgSpeciationArgRules   ) );
    
    ArgumentRules* avgExtinctionArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "averageExtinctionRate", variable, avgExtinctionArgRules   ) );

    ArgumentRules* prob_speciation_inc_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "probabilitySpeciationRateIncrease", variable, prob_speciation_inc_arg_rules   ) );

    ArgumentRules* prob_speciation_dec_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "probabilitySpeciationRateDecrease", variable, prob_speciation_dec_arg_rules   ) );

    ArgumentRules* prob_extinction_inc_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "probabilityExtinctionRateIncrease", variable, prob_extinction_inc_arg_rules   ) );

    ArgumentRules* prob_extinction_dec_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_conditionedBirthDeathShiftProcessContinuous, ModelVector<RealPos> >( "probabilityExtinctionRateDecrease", variable, prob_extinction_dec_arg_rules   ) );

    return methods;
}



/* Return member rules */
const MemberRules& Dist_conditionedBirthDeathShiftProcessContinuous::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "rootAge"      , RealPos::getClassTypeSpec()                   , "The age of the root."                                , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rootLambda"   , RealPos::getClassTypeSpec()                   , "The speciation rate at the root."                    , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rootMu"       , RealPos::getClassTypeSpec()                   , "The extinction rate at the root."                    , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "lambda"       , TypedDistribution<RealPos>::getClassTypeSpec(), "The prior distribution for the speciation rates."    , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "mu"           , TypedDistribution<RealPos>::getClassTypeSpec(), "The prior distribution for the extinction rates."    , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule( "delta"        , RealPos::getClassTypeSpec()                   , "The rate of jumping between rate categories."        , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rho"          , Probability::getClassTypeSpec()               , "The taxon sampling probability."                     , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        memberRules.push_back( new OptionRule( "condition"      , new RlString("survival"), optionsCondition    , "The condition of the birth-death process." ) );
        memberRules.push_back( new ArgumentRule( "taxa"         , ModelVector<Taxon>::getClassTypeSpec()        , "The taxon names used for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_conditionedBirthDeathShiftProcessContinuous::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_conditionedBirthDeathShiftProcessContinuous::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "rootAge" )
    {
        root_age = var;
    }
    else if ( name == "rootLambda" )
    {
        root_speciation = var;
    }
    else if ( name == "rootMu" )
    {
        root_extinction = var;
    }
    else if ( name == "lambda" )
    {
        speciation = var;
    }
    else if ( name == "mu" )
    {
        extinction = var;
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
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}



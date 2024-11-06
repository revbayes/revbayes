#include <math.h>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_heterogeneousRateBirthDeath.h"
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
#include "HeterogeneousRateBirthDeath.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
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


Dist_heterogeneousRateBirthDeath::Dist_heterogeneousRateBirthDeath() : TypedDistribution<TimeTree>()
{
    
}


Dist_heterogeneousRateBirthDeath::~Dist_heterogeneousRateBirthDeath()
{
    
}



Dist_heterogeneousRateBirthDeath* Dist_heterogeneousRateBirthDeath::clone( void ) const
{
    
    return new Dist_heterogeneousRateBirthDeath( *this );
}


RevBayesCore::HeterogeneousRateBirthDeath* Dist_heterogeneousRateBirthDeath::createDistribution( void ) const
{
    
    // Get the parameters
    RevBayesCore::TypedDagNode<double>* ra   = static_cast<const RealPos &>( root_age->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>*    rs   = static_cast<const Natural &>( root_state->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sp   = static_cast<const ModelVector<RealPos> &>( speciation->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ex   = static_cast<const ModelVector<RealPos> &>( extinction->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* er   = static_cast<const RealPos &>( event_rate->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rh   = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    bool allow_same   = static_cast<const RlBoolean &>( allow_same_category->getRevObject() ).getValue();
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    // condition
    const std::string& cond                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    RevBayesCore::HeterogeneousRateBirthDeath*   d = new RevBayesCore::HeterogeneousRateBirthDeath( ra, rs, sp, ex, er, rh, cond, allow_same, t );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_heterogeneousRateBirthDeath::getClassType(void)
{
    
    static std::string rev_type = "Dist_heterogeneousRateBirthDeath";
    
    return rev_type;
}


/* Get class type spec describing type of object. TODO: Check if the correct parent is TypedDistribution or Distribution */
const TypeSpec& Dist_heterogeneousRateBirthDeath::getClassTypeSpec(void)
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
std::string Dist_heterogeneousRateBirthDeath::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "HBDP";
    
    return d_name;
}


MethodTable Dist_heterogeneousRateBirthDeath::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();
    
    // member functions
    ArgumentRules* numEventsArgRules = new ArgumentRules();
    //    parentArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_heterogeneousRateBirthDeath, ModelVector<Natural> >( "numberEvents", variable, numEventsArgRules   ) );

    ArgumentRules* avgSpeciationArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_heterogeneousRateBirthDeath, ModelVector<RealPos> >( "averageSpeciationRate", variable, avgSpeciationArgRules   ) );

    ArgumentRules* avgExtinctionArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_heterogeneousRateBirthDeath, ModelVector<RealPos> >( "averageExtinctionRate", variable, avgExtinctionArgRules   ) );

    return methods;
}



/* Return member rules */
const MemberRules& Dist_heterogeneousRateBirthDeath::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "rootAge"   , RealPos::getClassTypeSpec()              , "The age of the root."                        , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rootState" , Natural::getClassTypeSpec()              , "The index of the rate category at the root." , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "lambda"    , ModelVector<RealPos>::getClassTypeSpec() , "The vector of speciation rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "mu"        , ModelVector<RealPos>::getClassTypeSpec() , "The vector of extinction rates."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "delta"     , RealPos::getClassTypeSpec()              , "The rate of jumping between rate categories.", ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rho"       , Probability::getClassTypeSpec()          , "The taxon sampling probability."             , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        memberRules.push_back( new OptionRule( "condition"    , new RlString("survival"), optionsCondition, "The condition of the birth-death process." ) );
        memberRules.push_back( new ArgumentRule( "allowSameCategory", RlBoolean::getClassTypeSpec()          , "Do we allow shifts into the same category", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule( "taxa"      , ModelVector<Taxon>::getClassTypeSpec(), "The taxon names used for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_heterogeneousRateBirthDeath::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_heterogeneousRateBirthDeath::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "rootAge" )
    {
        root_age = var;
    }
    else if ( name == "rootState" )
    {
        root_state = var;
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
    else if ( name == "allowSameCategory" )
    {
        allow_same_category = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}


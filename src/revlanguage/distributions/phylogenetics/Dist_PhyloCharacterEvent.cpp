#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_PhyloCharacterEvent.h"
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
#include "PhyloCharacterEventDistribution.h"
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


Dist_PhyloCharacterEvent::Dist_PhyloCharacterEvent() : TypedDistribution<CharacterHistory>()
{
    
}


Dist_PhyloCharacterEvent::~Dist_PhyloCharacterEvent()
{
    
}



Dist_PhyloCharacterEvent* Dist_PhyloCharacterEvent::clone( void ) const
{
    
    return new Dist_PhyloCharacterEvent( *this );
}


RevBayesCore::PhyloCharacterEventDistribution* Dist_PhyloCharacterEvent::createDistribution( void ) const
{
    
    // Get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>*             t   = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    const std::vector<std::string>&                             n   = static_cast<const ModelVector<RlString> &>( names->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double>>* r   = static_cast<const ModelVector<Real> &>( root_values->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>*                         s   = static_cast<const RealPos &>( event_rate->getRevObject() ).getDagNode();

    const WorkspaceVector<Distribution>& rl_bd = static_cast<const WorkspaceVector<Distribution > &>( base_distribution->getRevObject() );
    std::vector<RevBayesCore::TypedDistribution<double>* > bd;
    for (size_t i=0; i<rl_bd.size();++i)
    {
        RevBayesCore::TypedDistribution<double>* tmp = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_bd[i].createDistribution() );
        bd.push_back( tmp );
    }

    RevBayesCore::PhyloCharacterEventDistribution*   d = new RevBayesCore::PhyloCharacterEventDistribution( r, bd, t, s, n );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_PhyloCharacterEvent::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloCharacterEvent";
    
    return rev_type;
}


/* Get class type spec describing type of object. */
const TypeSpec& Dist_PhyloCharacterEvent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<CharacterHistory>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloCharacterEvent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloCharacterEvent";
    
    return d_name;
}


/**
 * Get the method table for this distribution.
 * We need to implement this function when a random variable drawn from this distribution
 * has specific member method. Here, these are:
 * - x.numberEvents()
 *
 * \return Rev name of constructor function.
 */
MethodTable Dist_PhyloCharacterEvent::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<CharacterHistory>::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_num_events_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_PhyloCharacterEvent, ModelVector<Natural> >( "numberEvents", variable, get_num_events_arg_rules   ) );
    
    // member functions
    ArgumentRules* get_real_values_arg_rules = new ArgumentRules();
    get_real_values_arg_rules->push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_PhyloCharacterEvent, ModelVector<Real> >( "getRealValues", this->variable, get_real_values_arg_rules, true ) );

    ArgumentRules* get_real_pos_values_arg_rules = new ArgumentRules();
    get_real_pos_values_arg_rules->push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_PhyloCharacterEvent, ModelVector<RealPos> >( "getRealPosValues", this->variable, get_real_pos_values_arg_rules, true ) );

    return methods;
}



/* Return member rules */
const MemberRules& Dist_PhyloCharacterEvent::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        member_rules.push_back( new ArgumentRule( "tree"                , Tree::getClassTypeSpec()                                          , "The tree on which we let the character events occur."    , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "eventRate"           , RealPos::getClassTypeSpec()                                       , "The rate of character events."                           , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "rootValues"          , ModelVector<Real>::getClassTypeSpec()                             , "The values at the root (also for character 0)."          , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "baseDistributions"   , WorkspaceVector<TypedDistribution<Real> >::getClassTypeSpec()     , "The prior distribution for the speciation rates."        , ArgumentRule::BY_CONSTANT_REFERENCE   , ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "names"               , ModelVector<RlString>::getClassTypeSpec()                         , "The names of the variables used for retrieval."          , ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return member_rules;
}


const TypeSpec& Dist_PhyloCharacterEvent::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_PhyloCharacterEvent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "eventRate" )
    {
        event_rate = var;
    }
    else if ( name == "rootValues" )
    {
        root_values = var;
    }
    else if ( name == "baseDistributions" )
    {
        base_distribution = var;
    }
    else if ( name == "names" )
    {
        names = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}



#include <math.h>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_AutocorrelatedEvent.h"
#include "ModelVector.h"
#include "AutocorrelatedEventDistribution.h"
#include "Natural.h"
#include "RlDistributionMemberFunction.h"
#include "RlString.h"
#include "StochasticNode.h"
#include "WorkspaceVector.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "MultiValueEvent.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlMultiValueEvent.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "WorkspaceToCoreWrapperObject.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_AutocorrelatedEvent::Dist_AutocorrelatedEvent() : TypedDistribution<MultiValueEvent>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_AutocorrelatedEvent* Dist_AutocorrelatedEvent::clone( void ) const
{
    
    return new Dist_AutocorrelatedEvent(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::AutocorrelatedEventDistribution* Dist_AutocorrelatedEvent::createDistribution( void ) const
{
    // get the parameters
//    const WorkspaceVector<TypedDistribution<Real> >& rl_vp                   = static_cast<const WorkspaceVector<TypedDistribution<Real> > &>( value_priors->getRevObject() );
    const WorkspaceVector<Distribution>& rl_vp                   = static_cast<const WorkspaceVector<Distribution > &>( value_priors->getRevObject() );
    std::vector<RevBayesCore::TypedDistribution<double>* > vp;
    for (size_t i=0; i<rl_vp.size();++i)
    {
        RevBayesCore::TypedDistribution<double>* tmp = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp[i].createDistribution() );
        vp.push_back( tmp );
    }
    
    const Distribution& rl_ep                                           = static_cast<const Distribution &>( event_number_prior->getRevObject() );
    RevBayesCore::TypedDistribution<long>* ep                           = static_cast<RevBayesCore::TypedDistribution<long>* >( rl_ep.createDistribution() );

    const std::vector<std::string>& n = static_cast<const ModelVector<RlString> &>( names->getRevObject() ).getValue();
    // offsets
    const std::vector<long> & m = static_cast<const ModelVector<Natural> &>( min_elements->getRevObject() ).getValue();

    
    const std::vector<std::string>& ac_type_names = static_cast<const ModelVector<RlString> &>( autocorrelation_types->getRevObject() ).getValue();
    std::vector< RevBayesCore::AutocorrelatedEventDistribution::Autocorrelation > ac_types = std::vector< RevBayesCore::AutocorrelatedEventDistribution::Autocorrelation >( ac_type_names.size(), RevBayesCore::AutocorrelatedEventDistribution::NONE );
    for ( size_t i=0; i<ac_type_names.size(); ++i )
    {
        if ( ac_type_names[i] == "NONE" )
        {
            ac_types[i] = RevBayesCore::AutocorrelatedEventDistribution::NONE;
        }
        else if ( ac_type_names[i] == "ACN" )
        {
            ac_types[i] = RevBayesCore::AutocorrelatedEventDistribution::ACN;
        }
        else if ( ac_type_names[i] == "ACLN" )
        {
            ac_types[i] = RevBayesCore::AutocorrelatedEventDistribution::ACLN;
        }
    }
    
    const std::vector<std::string>& ac_dep_var = static_cast<const ModelVector<RlString> &>( autocorrelation_dep_var->getRevObject() ).getValue();
    
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ac_sd = static_cast<const ModelVector<RealPos> &>( autocorrelation_sigmas->getRevObject() ).getDagNode();
    
    const std::string& sort_by = static_cast<const RlString &>( name_of_var_to_sort_by->getRevObject() ).getValue();

    RevBayesCore::AutocorrelatedEventDistribution* d    = new RevBayesCore::AutocorrelatedEventDistribution(ep, vp, ac_types, ac_dep_var, ac_sd, n, m, sort_by);
    
    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_AutocorrelatedEvent::getClassType(void)
{
    
    static std::string rev_type = "Dist_AutocorrelatedEvent";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_AutocorrelatedEvent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<MultiValueEvent>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_AutocorrelatedEvent::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_AutocorrelatedEvent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "AutocorrelatedEvent";
    
    return d_name;
}



MethodTable Dist_AutocorrelatedEvent::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<MultiValueEvent>::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_real_values_arg_rules = new ArgumentRules();
    get_real_values_arg_rules->push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_AutocorrelatedEvent, ModelVector<Real> >( "getRealValues", this->variable, get_real_values_arg_rules, true ) );

    ArgumentRules* get_real_pos_values_arg_rules = new ArgumentRules();
    get_real_pos_values_arg_rules->push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_AutocorrelatedEvent, ModelVector<RealPos> >( "getRealPosValues", this->variable, get_real_pos_values_arg_rules, true ) );

    ArgumentRules* get_num_events_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_AutocorrelatedEvent, Natural >( "getNumberOfEvents", this->variable, get_num_events_arg_rules, true ) );
    
    
    return methods;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the AutocorrelatedEvent distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_AutocorrelatedEvent::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "eventDistribution", TypedDistribution<Natural>::getClassTypeSpec(), "The prior on the number of events.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "valueDistribution", WorkspaceVector< Distribution >::getClassTypeSpec(), "The base distribution for the per category values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "names", ModelVector<RlString>::getClassTypeSpec(), "The names of the values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "minNumberEvents", ModelVector<Natural>::getClassTypeSpec(), "The minum number of values per value category.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "autocorrelationTypes", ModelVector<RlString>::getClassTypeSpec(), "The autocorrelation types per variables. Allowed values are NONE|ACN|ACLN.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "autocorrelationDependencies", ModelVector<RlString>::getClassTypeSpec(), "The autocorrelation dependencies per variables. This is important if you want autocorrelation to be time-dependent. Allowed values are NONE or the name of another variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "autocorrelationSigmas", ModelVector<RealPos>::getClassTypeSpec(), "The standard deviation parameters of the autocorrelated variables.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "sort", RlString::getClassTypeSpec(), "The name of the variable to sort by. For example, if you want to enforce sorting by time the you need to specify the name of your time variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlString("none") ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_AutocorrelatedEvent::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_AutocorrelatedEvent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
     
    if ( name == "eventDistribution" )
    {
        event_number_prior = var;
    }
    else if ( name == "valueDistribution" )
    {
        value_priors = var;
    }
    else if ( name == "names" )
    {
        names = var;
    }
    else if ( name == "minNumberEvents" )
    {
        min_elements = var;
    }
    else if ( name == "autocorrelationTypes" )
    {
        autocorrelation_types = var;
    }
    else if ( name == "autocorrelationDependencies" )
    {
        autocorrelation_dep_var = var;
    }
    else if ( name == "autocorrelationSigmas" )
    {
        autocorrelation_sigmas = var;
    }
    else if ( name == "sort" )
    {
        name_of_var_to_sort_by = var;
    }
    else
    {
        TypedDistribution<MultiValueEvent>::setConstParameter(name, var);
    }
    
}

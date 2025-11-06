#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BirthDeathBurstProcess.h"
#include "Dist_BirthDeathBurstProcess.h"
#include "RlDistributionMemberFunction.h"
#include "DeterministicNode.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "DagMemberFunction.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "MethodTable.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBirthDeathProcess.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "RlTaxon.h" // IWYU pragma: keep

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_BirthDeathBurstProcess::Dist_BirthDeathBurstProcess() : BirthDeathProcess()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_BirthDeathBurstProcess* Dist_BirthDeathBurstProcess::clone( void ) const
{
    return new Dist_BirthDeathBurstProcess(*this);
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
RevBayesCore::BirthDeathBurstProcess* Dist_BirthDeathBurstProcess::createDistribution( void ) const
{
    
    // get the parameters
    
    // the root age
    RevBayesCore::TypedDagNode<double>* ra          = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();
    
    // speciation rate
    RevBayesCore::TypedDagNode<double>* birth       = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();
    
    // extinction rate
    RevBayesCore::TypedDagNode<double>* death       = static_cast<const RealPos &>( mu->getRevObject() ).getDagNode();
    
    // burst probability
    RevBayesCore::TypedDagNode<double>* burst_prob  = static_cast<const Probability &>( beta->getRevObject() ).getDagNode();
    
    // burst time
    RevBayesCore::TypedDagNode<double>* bt          = static_cast<const RealPos &>( burst_time->getRevObject() ).getDagNode();
    
    
    // sampling probability
    RevBayesCore::TypedDagNode<double>* r           = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();
    
    // sampling strategy
//    const std::string &strategy                     = static_cast<const RlString &>( sampling_strategy->getRevObject() ).getValue();

    // condition
    const std::string& cond                         = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    
    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t              = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    
    // create the internal distribution object
    RevBayesCore::BirthDeathBurstProcess* d         = new RevBayesCore::BirthDeathBurstProcess(ra, birth, death, burst_prob, bt, r, cond, t);
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_BirthDeathBurstProcess::getClassType( void )
{
    
    static std::string rev_type = "Dist_BirthDeathBurstProcess";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_BirthDeathBurstProcess::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( BirthDeathProcess::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_BirthDeathBurstProcess::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "BDBP" );
    
    return a_names;
}


MethodTable Dist_BirthDeathBurstProcess::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();
    
    ArgumentRules* num_spec_events_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_BirthDeathBurstProcess, Natural >( "numBurstSpeciations", variable, num_spec_events_arg_rules   ) );

    return methods;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_BirthDeathBurstProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "BirthDeathBurstProcess";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the constant-rate birth-death process are:
 * (1) the speciation rate lambda which must be a positive real.
 * (2) the extinction rate mu that must be a positive real.
 * (3) all member rules specified by BirthDeathProcess.
 *
 * \return The member rules.
 */
const MemberRules& Dist_BirthDeathBurstProcess::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "lambda"     , RealPos::getClassTypeSpec(), "The speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu"         , RealPos::getClassTypeSpec(), "The piecewise-constant extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "beta"       , Probability::getClassTypeSpec(), "The constant extinction rate change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "burstTime"  , RealPos::getClassTypeSpec(), "The constant extinction rate change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "rho"        , Probability::getClassTypeSpec(), "The taxon sampling fraction(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = BirthDeathProcess::getParameterRules();
        dist_member_rules.insert(dist_member_rules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_BirthDeathBurstProcess::getTypeSpec( void ) const
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
void Dist_BirthDeathBurstProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "beta" )
    {
        beta = var;
    }
    else if ( name == "burstTime" )
    {
        burst_time = var;
    }
    else
    {
        BirthDeathProcess::setConstParameter(name, var);
    }
    
}

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_ConstrainedNodeAge.h"
#include "ModelVector.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "Distribution.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "NodeAgeConstrainedTreeDistribution.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Clade; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_ConstrainedNodeAge::Dist_ConstrainedNodeAge() : TypedDistribution<TimeTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_ConstrainedNodeAge* Dist_ConstrainedNodeAge::clone( void ) const
{
    return new Dist_ConstrainedNodeAge(*this);
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
RevBayesCore::NodeAgeConstrainedTreeDistribution* Dist_ConstrainedNodeAge::createDistribution( void ) const
{
    
    // get the parameters
    const RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::Clade> >& c               = static_cast<const ModelVector< ModelVector< Clade > > &>( constraints->getRevObject() ).getValue();
    const Distribution&                                                         rlDistribution  = static_cast<const Distribution &>( baseDistribution->getRevObject() );
    RevBayesCore::TypedDistribution<RevBayesCore::Tree>*                        base            = dynamic_cast<RevBayesCore::TypedDistribution<RevBayesCore::Tree>* >( rlDistribution.createDistribution() );
    
    // create the internal distribution object
    RevBayesCore::NodeAgeConstrainedTreeDistribution* d = new RevBayesCore::NodeAgeConstrainedTreeDistribution(base, c);
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_ConstrainedNodeAge::getClassType( void )
{
    
    static std::string rev_type = "Dist_ConstrainedNodeAge";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_ConstrainedNodeAge::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_ConstrainedNodeAge::getDistributionFunctionAliases( void ) const
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
std::string Dist_ConstrainedNodeAge::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "ConstrainedNodeAge";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 *
 * \return The member rules.
 */
const MemberRules& Dist_ConstrainedNodeAge::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "treeDistribution", TypedDistribution<TimeTree>::getClassTypeSpec(), "The base distribution for the tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "constraints", ModelVector< ModelVector<Clade> >::getClassTypeSpec(),"The node Age constraints.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_ConstrainedNodeAge::getTypeSpec( void ) const
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
void Dist_ConstrainedNodeAge::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "constraints" )
    {
        constraints = var;
    }
    else if ( name == "treeDistribution" )
    {
        baseDistribution = var;
    }
    else
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }
    
}

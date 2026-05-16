#include <math.h>
#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_empiricalTopology.h"
#include "EmpiricalTopologyDistribution.h"
#include "RlTraceTree.h"
#include "RlBranchLengthTree.h"
#include "TypeSpec.h"
#include "TypedDistribution.h"



using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_empiricalTopology::Dist_empiricalTopology() : TypedDistribution<BranchLengthTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_empiricalTopology* Dist_empiricalTopology::clone( void ) const
{

    return new Dist_empiricalTopology(*this);
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
RevBayesCore::EmpiricalTopologyDistribution* Dist_empiricalTopology::createDistribution( void ) const
{
    
    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    const RevBayesCore::TraceTree& tt  = static_cast<const TraceTree &>( tree_trace->getRevObject() ).getValue();

        
    RevBayesCore::EmpiricalTopologyDistribution* d = new RevBayesCore::EmpiricalTopologyDistribution( tt );
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_empiricalTopology::getClassType(void)
{
    
    static std::string rev_type = "Dist_empiricalTopology";
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_empiricalTopology::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<BranchLengthTree>::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_empiricalTopology::getDistributionFunctionName( void ) const
{

    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "EmpiricalTopology";
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the empirical topology distribution are:
 * (1) the number of taxa.
 * (2) the names of the taxa.
 *
 * \return The member rules.
 */
const MemberRules& Dist_empiricalTopology::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "treeTrace", TraceTree::getClassTypeSpec(), "The tree samples for computing the topology probabilities.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_empiricalTopology::getTypeSpec( void ) const
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
void Dist_empiricalTopology::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "treeTrace" )
    {
        tree_trace = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

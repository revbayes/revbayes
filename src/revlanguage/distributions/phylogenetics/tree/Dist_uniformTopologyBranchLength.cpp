#include <math.h>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_uniformTopologyBranchLength.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "StochasticNode.h"
#include "UniformTopologyBranchLengthDistribution.h"
#include "Clade.h"
#include "RlClade.h"
#include "RlTaxon.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBranchLengthTree.h"
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



using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_uniformTopologyBranchLength::Dist_uniformTopologyBranchLength() : TypedDistribution<BranchLengthTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_uniformTopologyBranchLength* Dist_uniformTopologyBranchLength::clone( void ) const
{
    
    return new Dist_uniformTopologyBranchLength(*this);
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
RevBayesCore::UniformTopologyBranchLengthDistribution* Dist_uniformTopologyBranchLength::createDistribution( void ) const
{
    
    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();
    
    RevBayesCore::Clade og;
    if ( outgroup != NULL && outgroup->getRevObject() != RevNullObject::getInstance())
    {
        og = static_cast<const Clade &>( outgroup->getRevObject() ).getValue();
    }
    
    const Distribution& rlDistribution              = static_cast<const Distribution &>( branch_length_prior->getRevObject() );
    RevBayesCore::TypedDistribution<double>* blp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rlDistribution.createDistribution() );
    
    
    bool r = false;
//    r = static_cast<const RlBoolean &>( rooted->getRevObject() ).getValue();
    
    RevBayesCore::UniformTopologyBranchLengthDistribution* d = new RevBayesCore::UniformTopologyBranchLengthDistribution( t, og, blp, r);
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_uniformTopologyBranchLength::getClassType(void)
{
    
    static std::string rev_type = "Dist_uniformTopologyBranchLength";
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_uniformTopologyBranchLength::getClassTypeSpec(void)
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
std::string Dist_uniformTopologyBranchLength::getDistributionFunctionName( void ) const
{
    
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "UniformTopologyBranchLength";
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the uniform TopologyBranchLength distribution are:
 * (1) the number of taxa.
 * (2) the names of the taxa.
 *
 * \return The member rules.
 */
const MemberRules& Dist_uniformTopologyBranchLength::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        member_rules.push_back( new ArgumentRule( "taxa"       , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa that will be used for the tips.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "outgroup"   , Clade::getClassTypeSpec(), "The clade (consisting of one or more taxa) used as an outgroup.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        member_rules.push_back( new ArgumentRule( "branchLengthDistribution", TypedDistribution<RealPos>::getClassTypeSpec(), "The base distribution for the branch lengths.",   ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//        memberRules.push_back( new ArgumentRule( "rooted",         RlBoolean::getClassTypeSpec(), "Is the distribution over rooted topologies?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        rules_set = true;
    }
    
    return member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_uniformTopologyBranchLength::getTypeSpec( void ) const
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
void Dist_uniformTopologyBranchLength::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "outgroup" )
    {
        outgroup = var;
    }
    else if ( name == "branchLengthDistribution" )
    {
        branch_length_prior = var;
    }
    else if ( name == "rooted" )
    {
        rooted = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}


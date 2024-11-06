#include <math.h>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Dist_ConstrainedUnrootedTopology.h"
#include "ModelVector.h"
#include "RlBoolean.h"
#include "RlClade.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
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
#include "StochasticNode.h"
#include "TopologyConstrainedTreeDistribution.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_ConstrainedUnrootedTopology::Dist_ConstrainedUnrootedTopology() : TypedDistribution<BranchLengthTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_ConstrainedUnrootedTopology* Dist_ConstrainedUnrootedTopology::clone( void ) const
{

    return new Dist_ConstrainedUnrootedTopology(*this);
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
RevBayesCore::TopologyConstrainedTreeDistribution* Dist_ConstrainedUnrootedTopology::createDistribution( void ) const
{

    // get the parameters

    const RevBayesCore::RbVector<RevBayesCore::Clade>& c      = static_cast<const ModelVector<Clade> &>( constraints->getRevObject() ).getValue();
    const Distribution& rlDistribution                        = static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<RevBayesCore::Tree>* base = static_cast<RevBayesCore::TypedDistribution<RevBayesCore::Tree>* >( rlDistribution.createDistribution() );
    
    // tree for initialization
    RevBayesCore::Tree* init = nullptr;
    if ( initial_tree->getRevObject() != RevNullObject::getInstance() )
    {
        init = static_cast<const BranchLengthTree &>( initial_tree->getRevObject() ).getDagNode()->getValue().clone();
    }
    
    /* Since an unrooted topology is not a time tree, its tips are not associated with times that would require checking.
     * Therefore, this distribution has no user-facing (= Rev-level) "age check precision" argument, but because
     * TopologyConstrainedTreeDistribution expects one, here we set it to 0 and pass it to the constructor.
     */
    std::int64_t pr = 0;
    
    // create the internal distribution object
    RevBayesCore::TopologyConstrainedTreeDistribution* dist = new RevBayesCore::TopologyConstrainedTreeDistribution(base, c, init, pr); // , bb);
    
    if (backbone == NULL && backbone->getRevObject() != RevNullObject::getInstance())
    {
        ; // do nothing
    }
    else if ( backbone->getRevObject().isType( ModelVector<BranchLengthTree>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::Tree> >* bb = NULL;
        bb = static_cast<const ModelVector<BranchLengthTree> &>( backbone->getRevObject() ).getDagNode();
        
        dist->setBackbone( NULL, bb );
    }
    else if ( backbone->getRevObject().isType( BranchLengthTree::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* bb = NULL;
        bb = static_cast<const BranchLengthTree&>( backbone->getRevObject() ).getDagNode();

        dist->setBackbone( bb, NULL );
    }
    
    return dist;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_ConstrainedUnrootedTopology::getClassType( void )
{
    
    static std::string rev_type = "Dist_ConstrainedUnrootedTopology";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_ConstrainedUnrootedTopology::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<BranchLengthTree>::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_ConstrainedUnrootedTopology::getDistributionFunctionAliases( void ) const
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
std::string Dist_ConstrainedUnrootedTopology::getDistributionFunctionName( void ) const
{

    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "ConstrainedTopology";
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
const MemberRules& Dist_ConstrainedUnrootedTopology::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        member_rules.push_back( new ArgumentRule( "treeDistribution", TypedDistribution<BranchLengthTree>::getClassTypeSpec(), "The base distribution for the tree.",   ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "constraints",      ModelVector<Clade>::getClassTypeSpec(),          "The topological constraints.",          ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<Clade>() ) );
        member_rules.push_back( new ArgumentRule( "initialTree" ,     BranchLengthTree::getClassTypeSpec() , "Instead of drawing a tree from the distribution, initialize distribution with this tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        std::vector<TypeSpec> backboneTypes;
        backboneTypes.push_back( BranchLengthTree::getClassTypeSpec() );
        backboneTypes.push_back( ModelVector<BranchLengthTree>::getClassTypeSpec() );
        member_rules.push_back( new ArgumentRule( "backbone", backboneTypes, "The backbone topological constraints.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        member_rules.push_back( new ArgumentRule( "inverse",          RlBoolean::getClassTypeSpec(),                   "Should the constraint be inverted?",    ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        rules_set = true;
    }
    return member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_ConstrainedUnrootedTopology::getTypeSpec( void ) const
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
void Dist_ConstrainedUnrootedTopology::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var)
{
    
    if ( name == "constraints" )
    {
        constraints = var;
    }
    else if ( name == "treeDistribution" )
    {
        base_distribution = var;
    }
    else if ( name == "backbone" )
    {
        backbone = var;
    }
    else if ( name == "inverse" )
    {
        invert_constraint = var;
    }
    else if ( name == "initialTree" )
    {
        initial_tree = var;
    }
    else
    {
        TypedDistribution<BranchLengthTree>::setConstParameter(name, var);
    }
}

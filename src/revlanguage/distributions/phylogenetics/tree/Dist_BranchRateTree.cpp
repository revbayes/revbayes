#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_BranchRateTree.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDistributionMemberFunction.h"
#include "ModelVector.h"
#include "RlTimeTree.h"
#include "RlTrace.h"
#include "RlTraceTree.h"
#include "BranchRateTreeDistribution.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "Distribution.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbVector.h"
#include "Real.h"
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
#include "TraceNumeric.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class TraceTree; }
namespace RevBayesCore { template <class valueType> class Trace; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_BranchRateTree::Dist_BranchRateTree() : TypedDistribution<BranchLengthTree>()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_BranchRateTree* Dist_BranchRateTree::clone( void ) const
{

    return new Dist_BranchRateTree(*this);
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
RevBayesCore::BranchRateTreeDistribution* Dist_BranchRateTree::createDistribution( void ) const
{

    // get the parameters
    const Distribution& rl_rate_prior                           = static_cast<const Distribution &>( branch_rate_prior->getRevObject() );
    RevBayesCore::TypedDistribution<double>* brp                = dynamic_cast<RevBayesCore::TypedDistribution<double>* >( rl_rate_prior.createDistribution() );

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tt          = static_cast<const TimeTree &>( time_tree->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode<double>* rbf = NULL;
    if( root_branch_fraction->getRevObject() != RevNullObject::getInstance() )
    {
        rbf = static_cast<const Probability &>( root_branch_fraction->getRevObject() ).getDagNode();
    }

    // create the internal distribution object
    RevBayesCore::BranchRateTreeDistribution* dist              = new RevBayesCore::BranchRateTreeDistribution(tt, brp, rbf);


    return dist;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_BranchRateTree::getClassType( void )
{

    static std::string rev_type = "Dist_BranchRateTree";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_BranchRateTree::getClassTypeSpec( void )
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<BranchLengthTree>::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_BranchRateTree::getDistributionFunctionAliases( void ) const
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
std::string Dist_BranchRateTree::getDistributionFunctionName( void ) const
{

    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "BranchRateTree";
    return d_name;
}


MethodTable Dist_BranchRateTree::getDistributionMethods( void ) const
{

//    const Distribution& rlDistribution = static_cast<const Distribution &>( baseDistribution->getRevObject() );

    MethodTable methods = TypedDistribution<BranchLengthTree>::getDistributionMethods();

    // member functions
    // ArgumentRules* sample_prob_arg_rules = new ArgumentRules();
    // sample_prob_arg_rules->push_back( new ArgumentRule( "log", RlBoolean::getClassTypeSpec(), "If we should return the log-transformed probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );
    // methods.addFunction( new DistributionMemberFunction<Dist_BranchRateTree, ModelVector<Real> >( "getSampleProbabilities", this->variable, sample_prob_arg_rules   ) );
    //
    // // member functions
    // ArgumentRules* branch_rates_arg_rules = new ArgumentRules();
    // branch_rates_arg_rules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the tree in the trace for which we want to get the branch rates.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    // methods.addFunction( new DistributionMemberFunction<Dist_BranchRateTree, ModelVector<RealPos> >( "getBranchRates", this->variable, branch_rates_arg_rules   ) );

    return methods;
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
const MemberRules& Dist_BranchRateTree::getParameterRules(void) const
{

    static MemberRules member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        member_rules.push_back( new ArgumentRule( "branchRatePrior", TypedDistribution<RealPos>::getClassTypeSpec(), "The prior distribution for the branch rates.",   ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "timeTree", TimeTree::getClassTypeSpec(), "The time tree", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "rootBranchFraction", Probability::getClassTypeSpec(), "The fraction of how much of the root branch is assigned to the left subtree.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );

        rules_set = true;
    }
    return member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_BranchRateTree::getTypeSpec( void ) const
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
void Dist_BranchRateTree::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var)
{

    if ( name == "branchRatePrior" )
    {
        branch_rate_prior = var;
    }
    else if ( name == "timeTree" )
    {
        time_tree = var;
    }
    else if ( name == "rootBranchFraction" )
    {
        root_branch_fraction = var;
    }
    else
    {
        TypedDistribution<BranchLengthTree>::setConstParameter(name, var);
    }
}

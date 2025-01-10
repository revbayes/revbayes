#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_empiricalTree.h"
#include "Natural.h"
#include "RlDistributionMemberFunction.h"
#include "RlTrace.h"
#include "RlTraceTree.h"
#include "StochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "EmpiricalDistribution.h"
#include "MethodTable.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTree.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TraceNumeric.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class Trace; }

using namespace RevLanguage;


Dist_empiricalTree::Dist_empiricalTree() : TypedDistribution<Tree>()
{
    
}


Dist_empiricalTree::~Dist_empiricalTree()
{
    
}


Dist_empiricalTree* Dist_empiricalTree::clone( void ) const
{
    
    return new Dist_empiricalTree( *this );
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
RevBayesCore::EmpiricalDistribution<RevBayesCore::Tree>* Dist_empiricalTree::createDistribution( void ) const
{
    
    // get the parameters
    
    // tree trace
    RevBayesCore::TraceTree* tt = &static_cast<const TraceTree &>( trace->getRevObject() ).getValue();

    RevBayesCore::Trace<double>* nt = NULL;
    if( density->getRevObject() != RevNullObject::getInstance() )
    {
        nt = &static_cast<const Trace &>( density->getRevObject() ).getValue();
    }

    // create the internal distribution object
    RevBayesCore::EmpiricalDistribution<RevBayesCore::Tree>* d = new RevBayesCore::EmpiricalDistribution<RevBayesCore::Tree>( tt, nt );
    
    return d;
}


/* Get Rev type of object */
const std::string& Dist_empiricalTree::getClassType(void)
{
    
    static std::string rev_type = "Dist_empiricalTree";
    
    return rev_type;
}


/* Get class type spec describing type of object.  */
const TypeSpec& Dist_empiricalTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Tree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_empiricalTree::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "EmpiricalTree";
    
    return d_name;
}


MethodTable Dist_empiricalTree::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<Tree>::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_tree_index_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_empiricalTree, Natural >( "getTreeIndex", variable, get_tree_index_arg_rules   ) );
    
    
    return methods;
}


/* Return member rules */
const MemberRules& Dist_empiricalTree::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        member_rules.push_back( new ArgumentRule( "trace", TraceTree::getClassTypeSpec(), "The trace of tree samples.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "density", Trace::getClassTypeSpec(), "Optional trace of probability density values for each tree.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );

        rules_set = true;
    }
    
    return member_rules;
}


const TypeSpec& Dist_empiricalTree::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_empiricalTree::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "trace" )
    {
        trace = var;
    }
    else if ( name == "density" )
    {
        density = var;
    }
    else {
        Distribution::setConstParameter(name, var);
    }
}

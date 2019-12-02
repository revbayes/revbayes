#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_empiricalTree.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDistributionMemberFunction.h"
#include "RlString.h"
#include "RlTrace.h"
#include "RlTraceTree.h"
#include "StochasticNode.h"

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_empiricalTree::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Will Freyman" );
    authors.push_back( "Sebastian Hoehna" );
    authors.push_back( "Will Pett" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_empiricalTree::getHelpDescription(void) const
{
    std::string description = "";
    description += "Creates a distribution of trees from a trace of trees.";
    
    return description;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_empiricalTree::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Read in tree trace\n";
    example += "tree_trace = readTreeTrace(\"output/my.trees\", burnin=0.25)\n";
    example += "\n";
    example += "# Create a distribution of trees\n";
    example += "tree ~ dnEmpiricalTree(tree_trace)\n";
    example += "\n";
    example += "# Add an MCMC move\n";
    example += "moves[1] = mvEmpiricalTree(tree)\n";
    
    return example;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_empiricalTree::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mvEmpiricalTree" );
    see_also.push_back( "treeTrace" );
    see_also.push_back( "readTreeTrace" );
    
    
    return see_also;
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

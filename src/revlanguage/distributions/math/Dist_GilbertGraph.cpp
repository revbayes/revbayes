#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_GilbertGraph.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "GilbertGraphDistribution.h"
#include "ArgumentRules.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;

Dist_GilbertGraph::Dist_GilbertGraph() : TypedDistribution<MatrixRealSymmetric>()
{
    
}


Dist_GilbertGraph::~Dist_GilbertGraph()
{
    
}



Dist_GilbertGraph* Dist_GilbertGraph::clone( void ) const
{
    return new Dist_GilbertGraph(*this);
}


RevBayesCore::GilbertGraphDistribution* Dist_GilbertGraph::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* p = static_cast<const RealPos&>( edge_probability->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( dim->getRevObject() ).getValue() );
    
    RevBayesCore::GilbertGraphDistribution* f = new RevBayesCore::GilbertGraphDistribution( n, p );
    
    return f;
    
}



/* Get Rev type of object */
const std::string& Dist_GilbertGraph::getClassType(void)
{
    
    static std::string rev_type = "Dist_GilbertGraph";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_GilbertGraph::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_GilbertGraph::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "GilbertGraph";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_GilbertGraph::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "p", RealPos::getClassTypeSpec(), "The probability that two vertices share an edge.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "n", Natural::getClassTypeSpec(), "The dimensions of the correlation matrix.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_GilbertGraph::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_GilbertGraph::printValue(std::ostream& o) const {
    
    o << "GilbertGraph(";
    o << ")";
}


/** Set a member variable */
void Dist_GilbertGraph::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "p" )
    {
        edge_probability = var;
    }
    else if ( name == "n" )
    {
        dim = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

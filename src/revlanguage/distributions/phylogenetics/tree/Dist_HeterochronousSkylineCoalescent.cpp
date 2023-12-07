#include <math.h>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Dist_HeterochronousSkylineCoalescent.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "PiecewiseConstantHeterochronousCoalescent.h"
#include "RealPos.h"
#include "RlClade.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
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
Dist_HeterochronousSkylineCoalescent::Dist_HeterochronousSkylineCoalescent() : TypedDistribution<TimeTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_HeterochronousSkylineCoalescent* Dist_HeterochronousSkylineCoalescent::clone( void ) const
{
    return new Dist_HeterochronousSkylineCoalescent(*this);
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
RevBayesCore::PiecewiseConstantHeterochronousCoalescent* Dist_HeterochronousSkylineCoalescent::createDistribution( void ) const
{
    
    // get the parameters
    
    // theta
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* th       = static_cast<const ModelVector<RealPos> &>( theta->getRevObject() ).getDagNode();
    // theta
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* ti       = static_cast<const ModelVector<RealPos> &>( times->getRevObject() ).getDagNode();
    // method
//    const std::string &m                        = static_cast<const RlString &>( method->getRevObject() ).getValue();
    // taxon names
    const std::vector<RevBayesCore::Taxon> tn   = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getDagNode()->getValue();
    // clade constraints
    const std::vector<RevBayesCore::Clade> &c   = static_cast<const ModelVector<Clade> &>( constraints->getRevObject() ).getValue();
    // create the internal distribution object
    RevBayesCore::PiecewiseConstantHeterochronousCoalescent*   d = new RevBayesCore::PiecewiseConstantHeterochronousCoalescent(th, ti, tn, c);
    
    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_HeterochronousSkylineCoalescent::getClassType( void )
{
    
    static std::string rev_type = "Dist_HeterochronousSkylineCoalescent";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_HeterochronousSkylineCoalescent::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_HeterochronousSkylineCoalescent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "HeterochronousCoalescentSkyline";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the diversity-dependent pure-birth process are:
 * (1) the initial speciation rate lambda which must be a positive real.
 * (2) the carrying capacity that must be a natural number.
 * (3) all member rules specified by BirthDeathProcess.
 *
 * \return The member rules.
 */
const MemberRules& Dist_HeterochronousSkylineCoalescent::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "theta"      , ModelVector<RealPos>::getClassTypeSpec(), "A vector of per interval population sizes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "times"      , ModelVector<RealPos>::getClassTypeSpec(), "A vector of times for the intervals, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "events" );
        optionsCondition.push_back( "uniform" );
        optionsCondition.push_back( "specified" );
        memberRules.push_back( new OptionRule( "method", new RlString("events"), optionsCondition, "The method how intervals are defined." ) );
        memberRules.push_back( new ArgumentRule( "taxa"       , ModelVector<Taxon>::getClassTypeSpec(), "The taxon names used for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "constraints", ModelVector<Clade>::getClassTypeSpec(), "The strictly enforced topology constraints.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<Clade>() ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = TypedDistribution<TimeTree>::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_HeterochronousSkylineCoalescent::getTypeSpec( void ) const
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
void Dist_HeterochronousSkylineCoalescent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "theta" )
    {
        theta = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "constraints" )
    {
        constraints = var;
    }
    else if ( name == "times" )
    {
        times = var;
    }
    else if ( name == "method" )
    {
        method = var;
    }
    else
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }
    
}

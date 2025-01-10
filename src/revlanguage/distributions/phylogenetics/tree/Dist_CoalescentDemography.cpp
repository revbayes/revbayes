#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Dist_CoalescentDemography.h"
#include "DemographyCoalescent.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlClade.h"
#include "RlDemographicFunction.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DemographicFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbException.h"
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
#include "WorkspaceToCoreWrapperObject.h"
#include "WorkspaceVector.h"

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_CoalescentDemography::Dist_CoalescentDemography() : TypedDistribution<TimeTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_CoalescentDemography* Dist_CoalescentDemography::clone( void ) const
{
    return new Dist_CoalescentDemography(*this);
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
RevBayesCore::DemographyCoalescent* Dist_CoalescentDemography::createDistribution( void ) const
{
    
    // get the parameters
    
    // intervals
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* iv      = static_cast<const ModelVector<RealPos> &>( change_points->getRevObject() ).getDagNode();
    // demographic functions
    const WorkspaceVector<DemographicFunction> &ws_vec_df   = static_cast<const WorkspaceVector<DemographicFunction> &>( demographies->getRevObject() );
    RevBayesCore::RbVector<RevBayesCore::DemographicFunction> df;
    for ( size_t i = 0; i < ws_vec_df.size(); ++i )
    {
        df.push_back( ws_vec_df[i].getValue() );
    }
    
    if ( (iv->getValue().size()+1) != df.size() )
    {
        throw RbException("You need to provide 1 more demographic function than change points.");
    }
    
    // taxon names
    const std::vector<RevBayesCore::Taxon> tn               = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getDagNode()->getValue();
    // clade constraints
    const std::vector<RevBayesCore::Clade> &c               = static_cast<const ModelVector<Clade> &>( constraints->getRevObject() ).getValue();
    
    
    // create the internal distribution object
    RevBayesCore::DemographyCoalescent*   d             = new RevBayesCore::DemographyCoalescent(iv, df, tn, c);
    
    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_CoalescentDemography::getClassType( void )
{
    
    static std::string rev_type = "Dist_CoalescentDemography";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_CoalescentDemography::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_CoalescentDemography::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "HeterochronousCoalescentDemography" );
//    a_names.push_back( "CoalescentDemography" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_CoalescentDemography::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "CoalescentDemography";
    
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
const MemberRules& Dist_CoalescentDemography::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "df", WorkspaceVector<DemographicFunction>::getClassTypeSpec(), "The vector of demographic functions.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "changePoints", ModelVector<RealPos>::getClassTypeSpec(), "The times when the demographic function changes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"       , ModelVector<Taxon>::getClassTypeSpec(), "The taxon names used for initialization.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "constraints", ModelVector<Clade>::getClassTypeSpec()   , "The topological constraints strictly enforced.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<Clade>() ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = TypedDistribution<TimeTree>::getParameterRules();
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
const TypeSpec& Dist_CoalescentDemography::getTypeSpec( void ) const
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
void Dist_CoalescentDemography::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "changePoints" )
    {
        change_points = var;
    }
    else if ( name == "df" )
    {
        demographies = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "constraints" )
    {
        constraints = var;
    }
    else
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }
    
}

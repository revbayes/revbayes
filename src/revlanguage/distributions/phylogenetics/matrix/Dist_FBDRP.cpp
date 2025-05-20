#include "Dist_FBDRP.h"

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "FossilizedBirthDeathRangeProcess.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_FBDRP::Dist_FBDRP() : FossilizedBirthDeathRangeProcess<MatrixReal>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_FBDRP* Dist_FBDRP::clone( void ) const
{
    return new Dist_FBDRP(*this);
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
RevBayesCore::FossilizedBirthDeathRangeProcess* Dist_FBDRP::createDistribution( void ) const
{
    // get the parameters
    
    // sampling condition
    const std::string& cond  = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    
    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // speciation rate
    RevBayesCore::DagNode* l = lambda->getRevObject().getDagNode();
    // extinction rate
    RevBayesCore::DagNode* m = mu->getRevObject().getDagNode();
    // fossilization rate
    RevBayesCore::DagNode* p = psi->getRevObject().getDagNode();

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // rate change times
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    bool c  = static_cast<const RlBoolean &>( complete->getRevObject() ).getValue();
    bool use_bds = static_cast<const RlBoolean &>( bds->getRevObject() ).getValue();
    bool re = static_cast<const RlBoolean &>( resample->getRevObject() ).getValue();

    RevBayesCore::FossilizedBirthDeathRangeProcess* d = new RevBayesCore::FossilizedBirthDeathRangeProcess(l, m, p, r, rt, cond, t, c, re, use_bds);
    
    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_FBDRP::getClassType( void )
{
    
    static std::string rev_type = "Dist_FBDRP";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_FBDRP::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<ModelVector<ModelVector<RealPos> > >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_FBDRP::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "FBDRMatrix" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_FBDRP::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "FossilizedBirthDeathRange";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the fossilized birth-death process are:
 * (1) the speciation rate lambda which must be a positive real.
 * (2) the extinction rate mu that must be a positive real.
 * (3) the fossil sampling rate psi that must be a positive real.
 * (4) the extant sampling rate rho that must be a positive real.
 *
 * \return The member rules.
 */
const MemberRules& Dist_FBDRP::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "BDS", RlBoolean::getClassTypeSpec(), "Assume complete lineage sampling? (BDS model of Silvestro et al. 2019)", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        // add the rules from the base class
        const MemberRules &parentRules = FossilizedBirthDeathRangeProcess<MatrixReal>::getParameterRules();
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
const TypeSpec& Dist_FBDRP::getTypeSpec( void ) const
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
void Dist_FBDRP::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "BDS" )
    {
        bds = var;
    }
    else
    {
       FossilizedBirthDeathRangeProcess<MatrixReal>::setConstParameter(name,var);
    }

}

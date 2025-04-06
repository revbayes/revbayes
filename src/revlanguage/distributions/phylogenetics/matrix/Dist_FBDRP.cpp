#include "Dist_FBDRP.h"

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_FBDRP.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "FossilizedBirthDeathRangeProcess.h"
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
Dist_FBDRP::Dist_FBDRP() : TypedDistribution<MatrixReal >()
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
    const std::string& cond                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    
    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // speciation rate
    RevBayesCore::DagNode* l = lambda->getRevObject().getDagNode();
    // extinction rate
    RevBayesCore::DagNode* m = mu->getRevObject().getDagNode();
    // fossilization rate
    RevBayesCore::DagNode* p = psi->getRevObject().getDagNode();

    // fossil counts
    RevBayesCore::DagNode* c = NULL;
    if ( fossil_counts->getRevObject() != RevNullObject::getInstance() )
    {
        c = fossil_counts->getRevObject().getDagNode();
    }

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r       = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // rate change times
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    bool pa = static_cast<const RlBoolean &>( presence_absence->getRevObject() ).getValue();

    RevBayesCore::FossilizedBirthDeathRangeProcess* d = new RevBayesCore::FossilizedBirthDeathRangeProcess(l, m, p, c, r, rt, cond, t, pa);

    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_FBDRP::getClassType( void )
{
    
    static std::string rev_type = "Dist_FBDRMatrix";
    
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
    std::string d_name = "FossilizedBirthDeathRangeMatrix";
    
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

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( RealPos::getClassTypeSpec() );
        paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "lambda",  paramTypes, "The speciation rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu",      paramTypes, "The extinction rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "psi",     paramTypes, "The fossil sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "rho",     Probability::getClassTypeSpec(), "The extant sampling fraction.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

        dist_member_rules.push_back( new ArgumentRule( "timeline",   ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process (from oldest to youngest).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<TypeSpec> intTypes;
        intTypes.push_back( Natural::getClassTypeSpec() );
        intTypes.push_back( ModelVector<Natural>::getClassTypeSpec() );
        intTypes.push_back( ModelVector<ModelVector<Natural> >::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "k",   intTypes, "The fossil observation counts (total or [interval] or [interval][species]).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        dist_member_rules.push_back( new OptionRule( "condition", new RlString("time"), optionsCondition, "The condition of the process." ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa with stratigraphic ranges used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        dist_member_rules.push_back( new ArgumentRule( "binary" , RlBoolean::getClassTypeSpec() , "Treat fossil counts as binary presence/absence data?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

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
    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "psi" )
    {
        psi = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else if ( name == "k" )
    {
        fossil_counts = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "binary" )
    {
        presence_absence = var;
    }
    else
    {
        TypedDistribution<MatrixReal >::setConstParameter(name, var);
    }
    
}

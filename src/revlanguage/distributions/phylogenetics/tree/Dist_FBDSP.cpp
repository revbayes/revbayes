#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MethodTable.h"
#include "StochasticNode.h"
#include "RlStochasticNode.h"
#include "DistributionMemberFunction.h"
#include "RlDistributionMemberFunction.h"
#include "Dist_FBDSP.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "AbstractBirthDeathProcess.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBirthDeathProcess.h"
#include "RlBoolean.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_FBDSP::Dist_FBDSP() : FossilizedBirthDeathRangeProcess<TimeTree>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_FBDSP* Dist_FBDSP::clone( void ) const
{
    return new Dist_FBDSP(*this);
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
RevBayesCore::FossilizedBirthDeathSpeciationProcess* Dist_FBDSP::createDistribution( void ) const
{
    // get the parameters

    // the start age
    RevBayesCore::TypedDagNode<double>* sa = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

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
    // anagnetic speciation rate
    RevBayesCore::DagNode* la = lambda_a->getRevObject().getDagNode();
    // symmetric speciation probability
    RevBayesCore::DagNode* b = beta->getRevObject().getDagNode();

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // rate change times
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    // the range process carries no reporting term; a dnFossilRecord
    // node supplies the reporting model
    bool ext = static_cast<const RlBoolean &>( extended->getRevObject() ).getValue();

    RevBayesCore::FossilizedBirthDeathSpeciationProcess* d = new RevBayesCore::FossilizedBirthDeathSpeciationProcess(sa, l, m, p, r, la, b, rt, cond, t, true, ext);

    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_FBDSP::getClassType( void )
{
    
    static std::string rev_type = "Dist_FBDSP";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_FBDSP::getClassTypeSpec( void )
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( FossilizedBirthDeathRangeProcess<TimeTree>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_FBDSP::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "FBDSP" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_FBDSP::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "FossilizedBirthDeathSpeciation";
    
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
const MemberRules& Dist_FBDSP::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "origin", RealPos::getClassTypeSpec(), "The origin time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // add the bare range-process rules from the base class (no reporting args; the reporting
        // model is supplied by a dnFossilRecord node conditioned on this range process)
        const MemberRules &parentRules = FossilizedBirthDeathRangeProcess<TimeTree>::getCoreParameterRules();
        dist_member_rules.insert(dist_member_rules.end(), parentRules.begin(), parentRules.end());

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( RealPos::getClassTypeSpec() );
        paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "lambda_a",  paramTypes, "The anagenetic speciation rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

        std::vector<TypeSpec> betaParamTypes;
        betaParamTypes.push_back( Probability::getClassTypeSpec() );
        betaParamTypes.push_back( ModelVector<Probability>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "beta",  betaParamTypes, "The probability of symmetric speciation.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "extended",  RlBoolean::getClassTypeSpec(), "Are the tips extinction times? If false the extinction times are marginalized out and each range ends at min(tau_K, youngest child birth).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_FBDSP::getTypeSpec( void ) const
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
void Dist_FBDSP::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "lambda_a" )
    {
        lambda_a = var;
    }
    else if ( name == "beta" )
    {
        beta = var;
    }
    else if ( name == "extended" )
    {
        extended = var;
    }
    else if ( name == "origin" )
    {
        start_age = var;
    }
    else
    {
        FossilizedBirthDeathRangeProcess<TimeTree>::setConstParameter(name, var);
    }
    
}


/**
 * The first (tau_1) and last (tau_K) appearances, and the times the range originates and ends. They
 * are internal to the distribution, so a deterministic node is the only way a monitor reaches them.
 */
RevLanguage::MethodTable Dist_FBDSP::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<TimeTree>::getDistributionMethods();

    ArgumentRules* first_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDSP, ModelVector<RealPos> >( "getFirstAppearances", variable, first_app_arg_rules, true ) );

    ArgumentRules* last_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDSP, ModelVector<RealPos> >( "getLastAppearances", variable, last_app_arg_rules, true ) );

    ArgumentRules* origin_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDSP, RealPos >( "getOrigin", variable, origin_arg_rules, true ) );

    ArgumentRules* origination_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDSP, ModelVector<RealPos> >( "getOriginationTimes", variable, origination_arg_rules, true ) );

    ArgumentRules* extinction_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDSP, ModelVector<RealPos> >( "getExtinctionTimes", variable, extinction_arg_rules, true ) );

    return methods;
}

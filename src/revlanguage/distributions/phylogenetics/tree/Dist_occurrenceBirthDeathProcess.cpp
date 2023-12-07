#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_occurrenceBirthDeathProcess.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "OccurrenceBirthDeathProcess.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "AbstractBirthDeathProcess.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBirthDeathProcess.h"
#include "RlConstantNode.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_occurrenceBirthDeathProcess::Dist_occurrenceBirthDeathProcess() : BirthDeathProcess()
{}





/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_occurrenceBirthDeathProcess* Dist_occurrenceBirthDeathProcess::clone( void ) const
{
    return new Dist_occurrenceBirthDeathProcess(*this);
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
RevBayesCore::AbstractBirthDeathProcess* Dist_occurrenceBirthDeathProcess::createDistribution( void ) const
{
    // get the parameters

    // the start age
    RevBayesCore::TypedDagNode<double>* sa          = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // the start condition
    bool                                uo          = ( start_condition == "originAge" ? true : false );

    // sampling condition
    const std::string& cond                         = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon>    tree        = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // tree for initialization
    RevBayesCore::Tree* init = NULL;
    if ( initial_tree->getRevObject() != RevNullObject::getInstance() )
    {
                                        init        = static_cast<const TimeTree &>( initial_tree->getRevObject() ).getDagNode()->getValue().clone();
    }
    // maximum number of hidden lineages
    RevBayesCore::TypedDagNode< long >* n           = static_cast<const Natural &>( maxHiddenLin->getRevObject() ).getDagNode();
    // birth rate
    RevBayesCore::DagNode*              l           = lambda->getRevObject().getDagNode();
    // death rate
    RevBayesCore::DagNode*              m           = mu->getRevObject().getDagNode();
    // serial sampling rate
    RevBayesCore::DagNode*              p           = psi->getRevObject().getDagNode();
    // treatment probability
    RevBayesCore::DagNode*              rm          = r->getRevObject().getDagNode();
    // occurrence rate
    RevBayesCore::DagNode*              o           = omega->getRevObject().getDagNode();
    // event sampling
    RevBayesCore::DagNode*              rh          = NULL;

    if ( rho->getRevObject() != RevNullObject::getInstance() )
    {
                                        rh          = rho->getRevObject().getDagNode();
    }
    // global rate change timeline
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* gt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
                                        gt          = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    // parameter-specific timelines
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* lt = NULL;
    if ( lambda_timeline->getRevObject() != RevNullObject::getInstance() )
    {
        lt = static_cast<const ModelVector<RealPos> &>( lambda_timeline->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* mt = NULL;
    if ( mu_timeline->getRevObject() != RevNullObject::getInstance() )
    {
        mt = static_cast<const ModelVector<RealPos> &>( mu_timeline->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* pt = NULL;
    if ( psi_timeline->getRevObject() != RevNullObject::getInstance() )
    {
        pt = static_cast<const ModelVector<RealPos> &>( psi_timeline->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( r_timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( r_timeline->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ot = NULL;
    if ( omega_timeline->getRevObject() != RevNullObject::getInstance() )
    {
        ot = static_cast<const ModelVector<RealPos> &>( omega_timeline->getRevObject() ).getDagNode();
    }


    // occurrence ages
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*  occAges = NULL;
    if ( occurrence_ages->getRevObject() != RevNullObject::getInstance() )
    {
                                        occAges     = static_cast<const ModelVector<Real> &>( occurrence_ages->getRevObject() ).getDagNode();
    }
    // boolean : use Mt, otherwise use Lt
    bool                                Mt          = static_cast<const RlBoolean &>( useMt->getRevObject() ).getValue();
    // boolean : verbose
    bool                                vb          = static_cast<const RlBoolean &>( verbose->getRevObject() ).getValue();

    RevBayesCore::AbstractBirthDeathProcess* d;

    d = new RevBayesCore::OccurrenceBirthDeathProcess(sa, l, m, p, rm, o, rh, gt, lt, mt, pt, rt, ot, cond, tree, uo, init, n, occAges, Mt, vb);

    return d;
}





/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_occurrenceBirthDeathProcess::getClassType( void )
{
    static std::string rev_type = "Dist_occurrenceBirthDeathProcess";

    return rev_type;
}





/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_occurrenceBirthDeathProcess::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( BirthDeathProcess::getClassTypeSpec() ) );

    return rev_type_spec;
}





/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_occurrenceBirthDeathProcess::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "OBDP" );

    return a_names;
}





/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_occurrenceBirthDeathProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PiecewiseOccurrenceBirthDeath";

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
const MemberRules& Dist_occurrenceBirthDeathProcess::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        dist_member_rules.push_back( new ArgumentRule( aliases, RealPos::getClassTypeSpec()    , "Start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( RealPos::getClassTypeSpec() );
        paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "lambda",  paramTypes, "Speciation/birth rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu",      paramTypes, "Extinction/death rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "psi",     paramTypes, "Serial sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "r",       paramTypes, "Probabilit(y|ies) of death upon sampling (treatment).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "omega",   paramTypes, "Occurrence sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

        std::vector<TypeSpec> other_event_paramTypes;
        other_event_paramTypes.push_back( ModelVector<Probability>::getClassTypeSpec() );

        std::vector<TypeSpec> event_sampling_paramTypes;
        event_sampling_paramTypes.push_back( Probability::getClassTypeSpec() );
        event_sampling_paramTypes.push_back( ModelVector<Probability>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "rho",               event_sampling_paramTypes, "Sampling fraction at present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "timeline",          ModelVector<RealPos>::getClassTypeSpec(), "Rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        dist_member_rules.push_back( new ArgumentRule( "lambdaTimeline",    ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the speciation/birth rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "muTimeline",        ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the extinction/death rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "psiTimeline",       ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the serial sampling rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "omegaTimeline",     ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the occurrence sampling rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "rTimeline",         ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the (serial) treatment probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "survival" );
        optionsCondition.push_back( "survival2" );
        dist_member_rules.push_back( new OptionRule( "condition",           new RlString("none"), optionsCondition, "Condition of the process on the survival of either 1 (survival) or 2 lineages (survival2) to the present" ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  ,            ModelVector<Taxon>::getClassTypeSpec(), "Taxa used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "initialTree" ,      TimeTree::getClassTypeSpec() , "Instead of drawing a tree from the distribution, initialize distribution with this tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "maxHiddenLin",      Natural::getClassTypeSpec(), "Maximum number of hidden lineages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(30) ) );
        dist_member_rules.push_back( new ArgumentRule( "occurrence_ages" ,  ModelVector<Real>::getClassTypeSpec() , "Fixed occurrence ages", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "useMt",             RlBoolean::getClassTypeSpec(), "If true computes densities with the Mt forward traversal algorithm otherwise uses Lt backward one.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        dist_member_rules.push_back( new ArgumentRule( "verbose",           RlBoolean::getClassTypeSpec(), "If true displays warnings and information messages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        rules_set = true;
    }

    return dist_member_rules;
}





/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_occurrenceBirthDeathProcess::getTypeSpec( void ) const
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
void Dist_occurrenceBirthDeathProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "r" )
    {
        r = var;
    }
    else if ( name == "omega" )
    {
        omega = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "rootAge" || name == "originAge" )
    {
        start_age = var;
        start_condition = name;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else if ( name == "lambdaTimeline" )
    {
        lambda_timeline = var;
    }
    else if ( name == "muTimeline" )
    {
        mu_timeline = var;
    }
    else if ( name == "psiTimeline" )
    {
        psi_timeline = var;
    }
    else if ( name == "rTimeline" )
    {
        r_timeline = var;
    }
    else if ( name == "omegaTimeline" )
    {
        omega_timeline = var;
    }
    else if ( name == "initialTree" )
    {
        initial_tree = var;
    }
    else if ( name == "maxHiddenLin" )
    {
        maxHiddenLin = var;
    }
    else if ( name == "occurrence_ages" )
    {
        occurrence_ages = var;
    }
    else if ( name == "useMt" )
    {
        useMt = var;
    }
    else if ( name == "verbose" )
    {
        verbose = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else
    {
        BirthDeathProcess::setConstParameter(name, var);
    }
}

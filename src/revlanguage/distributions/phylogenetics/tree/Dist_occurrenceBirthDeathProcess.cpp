#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "OccurrenceBirthDeathProcess.h"
#include "Dist_occurrenceBirthDeathProcess.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Probability.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "Real.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "AbstractBirthDeathProcess.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBirthDeathProcess.h"
//#include "Rltree.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbVector; }

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_occurrenceBirthDeathProcess::Dist_occurrenceBirthDeathProcess() : BirthDeathProcess()
{
}


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
    //
    // bool piecewise = false;
    //
    // if ( lambda->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    // {
    //     piecewise = true;
    // }
    //
    // if ( mu->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    // {
    //     piecewise = true;
    // }
    //
    // if ( psi->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    // {
    //     piecewise = true;
    // }
    //
    // if ( rho->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    // {
    //     piecewise = true;
    // }
    //
    // if ( piecewise )
    // {
    //     //throw(RbException("Piecewise Constant FBDP currently unsupported. Please use constant-rate process by specifying only constant-rate parameters."));
    //     // speciation rate
    //     RevBayesCore::DagNode* l = lambda->getRevObject().getDagNode();
    //     // extinction rate
    //     RevBayesCore::DagNode* m = mu->getRevObject().getDagNode();
    //     // serial sampling rate
    //     RevBayesCore::DagNode* p = psi->getRevObject().getDagNode();
    //     // taxon sampling fraction
    //     RevBayesCore::DagNode* r = rho->getRevObject().getDagNode();
    //
    //     // rate change times
    //     RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ht = NULL;
    //     RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* lt = NULL;
    //     RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* mt = NULL;
    //     RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* pt = NULL;
    //     RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    //
    //     if ( lambda_timeline->getRevObject() != RevNullObject::getInstance() )
    //     {
    //         lt = static_cast<const ModelVector<RealPos> &>( lambda_timeline->getRevObject() ).getDagNode();
    //     }
    //
    //     if ( mu_timeline->getRevObject() != RevNullObject::getInstance() )
    //     {
    //         mt = static_cast<const ModelVector<RealPos> &>( mu_timeline->getRevObject() ).getDagNode();
    //     }
    //
    //     if ( psi_timeline->getRevObject() != RevNullObject::getInstance() )
    //     {
    //         pt = static_cast<const ModelVector<RealPos> &>( psi_timeline->getRevObject() ).getDagNode();
    //     }
    //
    //     if ( rho_timeline->getRevObject() != RevNullObject::getInstance() )
    //     {
    //         rt = static_cast<const ModelVector<RealPos> &>( rho_timeline->getRevObject() ).getDagNode();
    //     }
    //
    //     if ( timeline->getRevObject() != RevNullObject::getInstance() )
    //     {
    //         ht = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    //     }
    //
    //     d = new RevBayesCore::PiecewiseConstantSerialSampledBirthDeathProcess(sa, l, m, p, r, ht, lt, mt, pt, rt, cond, t, uo, init);
    // }

    RevBayesCore::AbstractBirthDeathProcess* d;

    // Get the parameters :

    // the start age
    RevBayesCore::TypedDagNode<double>* sa       = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // speciation rate
    RevBayesCore::TypedDagNode<double>* l       = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();

    // extinction rate
    RevBayesCore::TypedDagNode<double>* m       = static_cast<const RealPos &>( mu->getRevObject() ).getDagNode();

    // sampling rate
    RevBayesCore::TypedDagNode<double>* p       = static_cast<const RealPos &>( psi->getRevObject() ).getDagNode();

    // occurrence rate
    RevBayesCore::TypedDagNode<double>* o       = static_cast<const RealPos &>( omega->getRevObject() ).getDagNode();

    // sampling probability at present
    RevBayesCore::TypedDagNode<double>* rh      = static_cast<const RealPos &>( rho->getRevObject() ).getDagNode();

    // removal probability
    RevBayesCore::TypedDagNode<double>* r       = static_cast<const RealPos &>( removalPr->getRevObject() ).getDagNode();

    // maximum number of hidden lineages
    RevBayesCore::TypedDagNode< long >* n       = static_cast<const Natural &>( maxHiddenLin->getRevObject() ).getDagNode();

    // sampling condition
    const std::string&                  cond    = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon>    tn      = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // occurrence ages
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*  occAges = NULL;
    if ( occurrence_ages->getRevObject() != RevNullObject::getInstance() )
    {
        occAges                                  = static_cast<const ModelVector<Real> &>( occurrence_ages->getRevObject() ).getDagNode();
    }
    // RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*  occAges = static_cast<const ModelVector<Real> &>( occurrence_ages->getRevObject() ).getDagNode();


    // the start condition
    bool                                uo      = ( start_condition == "originAge" ? true : false );

    // boolean : use Mt, otherwise use Lt
    bool                                mt      = static_cast<const RlBoolean &>( useMt->getRevObject() ).getValue();

    // tree for initialization
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tr = NULL;
    if ( initialTree->getRevObject() != RevNullObject::getInstance() )
    {
       tr                                       = static_cast<const TimeTree &>( initialTree->getRevObject() ).getDagNode();
    }


    d = new RevBayesCore::OccurrenceBirthDeathProcess(sa, l, m, p, o, rh, r, n, cond, tn, occAges, uo, mt, tr);

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
    std::string d_name = "OccurrenceBirthDeath";

    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the constant-rate birth-death process are:
 *
 * (0) all member rules specified by BirthDeathProcess.
 *
 * (1) the start time of the process can be either the root age or the origin age
 *
 * (2) the speciation rate lambda must be a positive real.
 * (3) the extinction rate mu must be a positive real.
 * (4) the fossil sampling rate psi must be a positive real.
 * (5) the occurrence rate omega must be a positive real.
 * (6) the vector of occurrence ages must contain reals.
 *
 * (7) the sampling fraction at present rho must be a real between 0 and 1.
 * (8) the removal probability removalPr must be a real between 0 and 1.

 * (9) the process can be conditioned either on the time or on the survival.
 *
 * (10) the initial tree must be a time tree.
 *
 * \return The member rules.
 */
const MemberRules& Dist_occurrenceBirthDeathProcess::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        dist_member_rules.push_back( new ArgumentRule( aliases,             RealPos::getClassTypeSpec()    , "The start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // std::vector<TypeSpec> paramTypes;
        // paramTypes.push_back( RealPos::getClassTypeSpec() );
        // paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "lambda",            RealPos::getClassTypeSpec(), "The speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu",                RealPos::getClassTypeSpec(), "The extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "psi",               RealPos::getClassTypeSpec(), "The fossil sampling rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "omega",             RealPos::getClassTypeSpec(), "The occurrence rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );


        // std::vector<TypeSpec> rho_paramTypes;
        // rho_paramTypes.push_back( Probability::getClassTypeSpec() );
        // rho_paramTypes.push_back( ModelVector<Probability>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "rho",               Probability::getClassTypeSpec(), "The sampling fraction at present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "removalPr",         Probability::getClassTypeSpec(), "The removal probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0) ) );

        dist_member_rules.push_back( new ArgumentRule( "maxHiddenLin",      Natural::getClassTypeSpec(), "The maximum number of hidden lineages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(30) ) );

        // dist_member_rules.push_back( new ArgumentRule( "timeline",          ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        // dist_member_rules.push_back( new ArgumentRule( "lambdaTimes",       ModelVector<RealPos>::getClassTypeSpec(), "The speciation rate change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        // dist_member_rules.push_back( new ArgumentRule( "muTimes",           ModelVector<RealPos>::getClassTypeSpec(), "The extinction rate change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        // dist_member_rules.push_back( new ArgumentRule( "psiTimes",          ModelVector<RealPos>::getClassTypeSpec(), "The serial sampling rate change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        // dist_member_rules.push_back( new ArgumentRule( "rhoTimes",          ModelVector<RealPos>::getClassTypeSpec(), "The episodic taxon sampling fraction change times.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        dist_member_rules.push_back( new OptionRule( "condition",           new RlString("time"), optionsCondition, "The condition of the process." ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  ,            ModelVector<Taxon>::getClassTypeSpec(), "The taxa used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        dist_member_rules.push_back( new ArgumentRule( "occurrence_ages" ,  ModelVector<Real>::getClassTypeSpec() , "The fixed occurrence ages", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "useMt",             RlBoolean::getClassTypeSpec(), "If true computes densities with the Mt forward traversal algorithm otherwise uses Lt backward one.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( true ) ) );

        dist_member_rules.push_back( new ArgumentRule( "initialTree" ,      TimeTree::getClassTypeSpec() , "Instead of drawing a tree from the distribution, initialize distribution with this tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

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
    if ( name == "rootAge" || name == "originAge" )
    {
        start_age = var;
        start_condition = name;
    }
    else if ( name == "lambda" )
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
    else if ( name == "omega" )
    {
       omega = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "removalPr" )
    {
        removalPr = var;
    }
    else if ( name == "maxHiddenLin" )
    {
        maxHiddenLin = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "occurrence_ages" )
    {
        occurrence_ages = var;
    }
    else if ( name == "useMt" )
    {
        useMt = var;
    }
    // else if ( name == "lambdaTimes" )
    // {
    //     lambda_timeline = var;
    // }
    // else if ( name == "muTimes" )
    // {
    //     mu_timeline = var;
    // }
    // else if ( name == "psiTimes" )
    // {
    //     psi_timeline = var;
    // }
    // else if ( name == "rhoTimes" )
    // {
    //     rho_timeline = var;
    // }

    else if ( name == "initialTree" )
    {
        initialTree = var;
    }

    else
    {
        BirthDeathProcess::setConstParameter(name, var);
    }

}

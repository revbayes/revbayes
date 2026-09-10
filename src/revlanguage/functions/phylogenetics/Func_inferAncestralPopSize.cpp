#include "Func_inferAncestralPopSize.h"
#include "InferAncestralPopSizeFunction.h"
#include "InferAncestralPopSizeFunctionPiecewise.h"


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "Natural.h"
#include "RlMatrixReal.h"
#include "MatrixReal.h"
#include "RlTypedFunction.h"
#include "RlString.h"
#include "RlTree.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { template <class valueType> class RbVector; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

/** default constructor */
Func_inferAncestralPopSize::Func_inferAncestralPopSize (void) : TypedFunction<MatrixReal>( )
{}

/** clone function */
Func_inferAncestralPopSize* Func_inferAncestralPopSize::clone( void ) const
{
    return new Func_inferAncestralPopSize( *this );
}

/** you probably have to change the first line of this function */
RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_inferAncestralPopSize::createFunction( void ) const
{
    /**first we get parameters that are shared by the piecewise and constant rate models */

    RevBayesCore::TypedDagNode< double >*                           sa              = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // sampling condition
    const std::string&                                              cdt             = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // occurrence ages
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> > *  occAges         = static_cast<const ModelVector<Real> &>( occurrence_ages->getRevObject() ).getDagNode();

    std::vector<double>                                             tau             = static_cast<const ModelVector<Real> &>( time_points->getRevObject() ).getValue();

    // verbose
    bool                                                            vb              = static_cast<const RlBoolean &>( verbose->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< double >*                           rh              = static_cast<const RealPos &>( rho->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< std::int64_t >*                             n               = static_cast<const Natural &>( maxHiddenLin->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< RevBayesCore::Tree >*               tr              = static_cast<const TimeTree &>( timeTree->getRevObject() ).getDagNode();

    /**If a timeline is provided, go to the piecwise version*/

    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >*time            = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();

        RevBayesCore::DagNode*                                      l               = ( lambda->getRevObject() ).getDagNode();

        RevBayesCore::DagNode*                                      m               = ( mu->getRevObject() ).getDagNode();

        RevBayesCore::DagNode*                                      p               = ( psi->getRevObject() ).getDagNode();

        RevBayesCore::DagNode*                                      o               = ( omega->getRevObject() ).getDagNode();

        RevBayesCore::DagNode*                                      r               = ( removalPr->getRevObject() ).getDagNode();

        RevBayesCore::InferAncestralPopSizeFunctionPiecewise* fxn = new RevBayesCore::InferAncestralPopSizeFunctionPiecewise( sa, l, m, p, o, rh, r, n, cdt, occAges, tau, tr, time );

        return fxn;
    }

    RevBayesCore::TypedDagNode< double >*                           l               = static_cast<const RealPos &>( lambda->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< double >*                           m               = static_cast<const RealPos &>( mu->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< double >*                           p               = static_cast<const RealPos &>( psi->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< double >*                           o               = static_cast<const RealPos &>( omega->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< double >*                           r               = static_cast<const RealPos &>( removalPr->getRevObject() ).getDagNode();

    RevBayesCore::InferAncestralPopSizeFunction* fxn = new RevBayesCore::InferAncestralPopSizeFunction( sa, l, m, p, o, rh, r, n, cdt, occAges, tau, vb, tr );

    return fxn;
}


/* Get argument rules */
const ArgumentRules& Func_inferAncestralPopSize::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        argumentRules.push_back( new ArgumentRule( aliases,             RealPos::getClassTypeSpec(), "Start age of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( RealPos::getClassTypeSpec() );
        paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "lambda",            paramTypes, "Speciation/birth rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu",                paramTypes, "Extinction/death rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "psi",               paramTypes, "Serial sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "omega",             paramTypes, "Occurrence sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "rho",               Probability::getClassTypeSpec(), "Sampling probability at present time.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );
        argumentRules.push_back( new ArgumentRule( "removalPr",         Probability::getClassTypeSpec(), "Probabilit(y|ies) of death upon sampling (treatment).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "maxHiddenLin",      Natural::getClassTypeSpec(), "Maximum number of hidden lineages (algorithm accuracy).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(30) ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "survival" );
        optionsCondition.push_back( "survival2" );
        argumentRules.push_back( new OptionRule( "condition",           new RlString("none"), optionsCondition, "Condition of the process on the survival of either 1 (survival) or 2 lineages (survival2) to the present." ) );

        argumentRules.push_back( new ArgumentRule( "occurrence_ages",   ModelVector<Real>::getClassTypeSpec(), "Vector of occurrence ages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "time_points",       ModelVector<Real>::getClassTypeSpec(), "Time points at which we compute the density.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "timeTree" ,         TimeTree::getClassTypeSpec(), "Tree for which ancestral pop. size has to be computed.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "timeline",          ModelVector<RealPos>::getClassTypeSpec(), "Rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        argumentRules.push_back( new ArgumentRule( "verbose",           RlBoolean::getClassTypeSpec(), "If true displays warnings and information messages.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        rules_set = true;
    }

    return argumentRules;
}

const std::string& Func_inferAncestralPopSize::getClassType(void)
{
    static std::string rev_type = "Func_inferAncestralPopSize";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_inferAncestralPopSize::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_inferAncestralPopSize::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnInferAncestralPopSize";

    return f_name;
}

const TypeSpec& Func_inferAncestralPopSize::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
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
void Func_inferAncestralPopSize::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "time_points" )
    {
        time_points = var;
    }
    else if ( name == "timeTree" )
    {
        timeTree = var;
    }
    else if ( name == "verbose" )
    {
        verbose = var;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else {
      Function::setConstParameter(name, var);
    }
}

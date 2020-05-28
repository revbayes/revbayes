#include "Func_inferAncestralPopSize.h"
#include "InferAncestralPopSizeFunction.h"
#include "InferAncestralPopSizeFunctionPiecewise.h"


// you can probably get rid of some of these
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
{

}

/** clone function */
Func_inferAncestralPopSize* Func_inferAncestralPopSize::clone( void ) const
{
    return new Func_inferAncestralPopSize( *this );
}

/** you probably have to change the first line of this function */
RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_inferAncestralPopSize::createFunction( void ) const
{

/**first we get parameters that are shared by the piecewise and constant rate models */

  RevBayesCore::TypedDagNode< double >*                           sa              = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

  // sampling condition
  const std::string&                                              cdt             = static_cast<const RlString &>( this->args[8].getVariable()->getRevObject() ).getValue();

  // occurrence ages
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> > *  O               = static_cast<const ModelVector<Real> &>( this->args[9].getVariable()->getRevObject() ).getDagNode();

  std::vector<double>                                             tau             = static_cast<const ModelVector<Real> &>( this->args[10].getVariable()->getRevObject() ).getValue();

  bool                                                            uo              = ( start_condition == "originAge" ? true : false );

  // verbose
  bool                                                            vb              = static_cast<const RlBoolean &>( this->args[13].getVariable()->getRevObject() ).getValue();

  RevBayesCore::TypedDagNode< double >*                           rh              = static_cast<const RealPos &>( this->args[5].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< long >*                             n               = static_cast<const Natural &>( this->args[7].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< RevBayesCore::Tree >*               tr              = static_cast<const TimeTree &>( this->args[11].getVariable()->getRevObject() ).getDagNode();

  /**If a timeline is provided, go to the piecwise version*/

  if ( this->args[12].getVariable()->getRevObject() != RevNullObject::getInstance() )
  {
      RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >*time = static_cast<const ModelVector<RealPos> &>( this->args[12].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::DagNode* l = ( this->args[1].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::DagNode* m = ( this->args[2].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::DagNode* p = ( this->args[3].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::DagNode* o = ( this->args[4].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::DagNode* r = ( this->args[6].getVariable()->getRevObject() ).getDagNode();

      RevBayesCore::InferAncestralPopSizeFunctionPiecewise* fxn = new RevBayesCore::InferAncestralPopSizeFunctionPiecewise( sa, l, m, p, o, rh, r, n, cdt, O, tau, uo, tr, time );

      return fxn;

  }


  RevBayesCore::TypedDagNode< double >*                           l               = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< double >*                           m               = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< double >*                           p               = static_cast<const RealPos &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< double >*                           o               = static_cast<const RealPos &>( this->args[4].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< double >*                           r               = static_cast<const RealPos &>( this->args[6].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::InferAncestralPopSizeFunction* fxn = new RevBayesCore::InferAncestralPopSizeFunction( sa, l, m, p, o, rh, r, n, cdt, O, tau, uo, vb, tr );

  return fxn;


}


/* Get argument rules */
/* you probably have to edit the argument rules here */
/* but after this I don't think you need to change anything else in this file */
const ArgumentRules& Func_inferAncestralPopSize::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        argumentRules.push_back( new ArgumentRule( aliases,             RealPos::getClassTypeSpec(), "The start age of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( RealPos::getClassTypeSpec() );
        paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "lambda",            paramTypes, "The speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu",                paramTypes, "The extinction rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "psi",               paramTypes, "The fossil sampling rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "omega",             paramTypes, "The occurrence rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "rho",               Probability::getClassTypeSpec(), "The sampling fraction at present.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );
        argumentRules.push_back( new ArgumentRule( "removalPr",         Probability::getClassTypeSpec(), "The removal probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "maxHiddenLin",      Natural::getClassTypeSpec(), "The number of hidden lineages (algorithm accuracy).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(30) ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        argumentRules.push_back( new OptionRule( "condition",           new RlString("time"), optionsCondition, "The condition of the process." ) );

        argumentRules.push_back( new ArgumentRule( "occurrence_ages",   ModelVector<Real>::getClassTypeSpec(), "Occurrence ages for incomplete fossils.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "time_points",       ModelVector<Real>::getClassTypeSpec(), "Time points for which we compute density.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "timeTree" ,         TimeTree::getClassTypeSpec(), "Tree for which ancestral pop. size has to be computed.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "timeline",    ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

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

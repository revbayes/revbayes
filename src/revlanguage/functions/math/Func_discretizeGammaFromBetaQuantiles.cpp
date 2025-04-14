#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_discretizeGammaFromBetaQuantiles.h"
#include "Integer.h"
#include "NormalizeVectorFunction.h"
#include "ModelVector.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

#include "DiscretizeGammaFromBetaQuantilesFunction.h"


using namespace RevLanguage;

Func_discretizeGammaFromBetaQuantiles::Func_discretizeGammaFromBetaQuantiles() : TypedFunction< ModelVector< RealPos > >() {

}

/* Clone object */
Func_discretizeGammaFromBetaQuantiles* Func_discretizeGammaFromBetaQuantiles::clone( void ) const
{

    return new Func_discretizeGammaFromBetaQuantiles( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizeGammaFromBetaQuantiles::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<double>* alpha = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* beta = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* shape = static_cast<const Real &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* rate = static_cast<const Real &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* numCats = static_cast<const Integer &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
    bool median = static_cast<const RlBoolean &>( args[5].getVariable()->getRevObject() ).getValue();

    RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction *func = new RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction( alpha, beta, shape, rate, numCats, median );

    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizeGammaFromBetaQuantiles::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "alpha"   , RealPos::getClassTypeSpec()  , "The alpha parameter of the beta distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"    , RealPos::getClassTypeSpec()  , "The beta parameter of the beta distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "shape"   , RealPos::getClassTypeSpec()  , "The shape parameter of the gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rate"    , RealPos::getClassTypeSpec()  , "The rate parameter (rate = 1/scale) of the gamma distribution", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCats", Integer::getClassTypeSpec()   , "The number of categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "median"  , RlBoolean::getClassTypeSpec(), "Should we use the median or mean?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }

    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizeGammaFromBetaQuantiles::getClassName(void)
{

    static std::string rbClassName = "Func_discretizeGammaFromBetaQuantiles";

	return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizeGammaFromBetaQuantiles::getClassTypeSpec(void)
{

    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rbClass;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_discretizeGammaFromBetaQuantiles::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDiscretizeGammaFromBetaQuantiles";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_discretizeGammaFromBetaQuantiles::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}

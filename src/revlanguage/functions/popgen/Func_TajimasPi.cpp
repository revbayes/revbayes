#include "TajimasPiFunction.h"
#include "Func_TajimasPi.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_TajimasPi::Func_TajimasPi( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_TajimasPi* Func_TajimasPi::clone( void ) const
{
    
    return new Func_TajimasPi( *this );
}


RevBayesCore::TypedFunction< double >* Func_TajimasPi::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData >* d = static_cast<const AbstractHomologousDiscreteCharacterData &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Boolean >* ps = static_cast<const RlBoolean &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Boolean >* ex = static_cast<const RlBoolean &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TajimasPiFunction* f = new RevBayesCore::TajimasPiFunction( d, ps->getValue(), ex->getValue() );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_TajimasPi::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "data",    AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The character data matrix for which to compute the summary.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "perSite", RlBoolean::getClassTypeSpec(), "Is the statistic normalized per site?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "excludeAmbiguous", RlBoolean::getClassTypeSpec(), "Should we exclude ambiguous or missing characters?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_TajimasPi::getClassType(void)
{
    
    static std::string rev_type = "Func_TajimasPi";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_TajimasPi::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<RealPos>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_TajimasPi::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTajimasPi";
    
    return f_name;
}


const TypeSpec& Func_TajimasPi::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

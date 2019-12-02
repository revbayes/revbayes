#include "Func_distanceRateModifier.h"
#include "DistanceRateModifier.h"
#include "DistanceRateModifierFunction.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlAtlas.h"
#include "RlBoolean.h"
#include "RlCharacterHistoryRateModifier.h"
#include "RlDeterministicNode.h"
#include "RlRateGeneratorSequence.h"
#include "RlSimplex.h"
#include "TimeAtlas.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_distanceRateModifier::Func_distanceRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_distanceRateModifier* Func_distanceRateModifier::clone( void ) const {
    
    return new Func_distanceRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_distanceRateModifier::createFunction( void ) const
{
    
    const RevBayesCore::TimeAtlas* atlas = &( static_cast<const RlAtlas&>( this->args[0].getVariable()->getRevObject() ).getValue() );
    
    RevBayesCore::TypedDagNode<double>* dp = static_cast<const Real&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    bool udd    = static_cast<const RlBoolean &>( this->args[2].getVariable()->getRevObject() ).getValue();
    bool uav    = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getValue();
    bool uadj   = false; // static_cast<const RlBoolean &>( this->args[4].getVariable()->getRevObject() ).getValue();


    RevBayesCore::DistanceRateModifierFunction* f = new RevBayesCore::DistanceRateModifierFunction(dp, atlas, uadj, uav, udd);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_distanceRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "distances"    , ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "Matrix of distances", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//        argumentRules.push_back( new ArgumentRule( "atlas"        , RlAtlas::getClassTypeSpec()   , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "distancePower", Real::getClassTypeSpec()      , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(1e-5) ) );
        argumentRules.push_back( new ArgumentRule( "useDistances" , RlBoolean::getClassTypeSpec() , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );
        argumentRules.push_back( new ArgumentRule( "useAvailable" , RlBoolean::getClassTypeSpec() , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );
//        argumentRules.push_back( new ArgumentRule( "useAdjacency" , RlBoolean::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_distanceRateModifier::getClassType(void)
{
    
    static std::string rev_type = "Func_biogeo_grm";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_distanceRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_distanceRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDistanceRateModifier";
    
    return f_name;
}


const TypeSpec& Func_distanceRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "GeoSSECladogeneticBirthDeathFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "Func_GeoSSECladoEvents.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_GeoSSECladoEvents::Func_GeoSSECladoEvents( void ) : TypedFunction<CladogeneticSpeciationRateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_GeoSSECladoEvents* Func_GeoSSECladoEvents::clone( void ) const {
    
    return new Func_GeoSSECladoEvents( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix >* Func_GeoSSECladoEvents::createFunction( void ) const
{
    
    // sympatry rate vector
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // allopatry rate vector
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ar = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        
    RevBayesCore::GeoSSECladogeneticBirthDeathFunction* f = new RevBayesCore::GeoSSECladogeneticBirthDeathFunction( sr, ar );

    // dispersal-rate matrix (if it exists)
    if (this->args[2].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* jr = static_cast<const ModelVector<ModelVector<RealPos>> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
        f->setJumpRates(jr);
    }

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_GeoSSECladoEvents::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    { 
        argumentRules.push_back( new ArgumentRule( "sympatryRates",  ModelVector<RealPos>::getClassTypeSpec() , "The vector of sympatric speciation rates (per geographic region).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "allopatryRates", ModelVector<RealPos>::getClassTypeSpec() , "The vector of allopatric speciation rates (per geographic region).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "jumpRates",      ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "Jump-dispersal rate between pairs of areas.", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT, NULL ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_GeoSSECladoEvents::getClassType(void)
{
    
    static std::string rev_type = "Func_GeoSSECladoEvents";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_GeoSSECladoEvents::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_GeoSSECladoEvents::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGeoSSECladoEvents";
    
    return f_name;
}


const TypeSpec& Func_GeoSSECladoEvents::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

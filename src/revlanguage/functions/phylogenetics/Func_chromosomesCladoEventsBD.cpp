//
//  Func_chromosomesCladoEventsBD.cpp
//
//  Created by Will Freyman on 6/22/16.
//

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ChromosomesCladogenicBirthDeathFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "Func_chromosomesCladoEventsBD.h"
#include "RealPos.h"
#include "ModelVector.h"
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
Func_chromosomesCladoEventsBD::Func_chromosomesCladoEventsBD( void ) : TypedFunction<CladogeneticSpeciationRateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chromosomesCladoEventsBD* Func_chromosomesCladoEventsBD::clone( void ) const {
    
    return new Func_chromosomesCladoEventsBD( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix >* Func_chromosomesCladoEventsBD::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    int ns = (int)static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::ChromosomesCladogenicBirthDeathFunction* f = new RevBayesCore::ChromosomesCladogenicBirthDeathFunction( sr, ns );
    
    if (this->args[1].getVariable()->getRevObject() != RevNullObject::getInstance())
    {    
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rm = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
        if (rm->getValue().size() > 1)
        {
            throw RbException("Only one hidden rate category currently implemented.");
        }
        f->setRateMultipliers(rm);
    }
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_chromosomesCladoEventsBD::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "speciation_rates", ModelVector<RealPos>::getClassTypeSpec() , "The speciation rates for different cladogenetic event types.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rate_multipliers", ModelVector<RealPos>::getClassTypeSpec() , "The rate multipliers for hidden rate classes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "max_chromosomes", Natural::getClassTypeSpec(), "The maximum number of chromosomes.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chromosomesCladoEventsBD::getClassType(void)
{
    
    static std::string rev_type = "Func_chromosomesCladoEventsBD";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_chromosomesCladoEventsBD::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_chromosomesCladoEventsBD::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChromosomesCladoEventsBD";
    
    return f_name;
}


const TypeSpec& Func_chromosomesCladoEventsBD::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

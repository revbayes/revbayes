//
//  Func_chromosomesPloidyCladoEventsBD.cpp
//
//  Created by Will Freyman on 5/22/17.
//

#include <iosfwd>
#include <string>
#include <vector>

#include "ChromosomesPloidyCladogenicBirthDeathFunction.h"
#include "Func_chromosomesPloidyCladoEventsBD.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "Natural.h"
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

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_chromosomesPloidyCladoEventsBD::Func_chromosomesPloidyCladoEventsBD( void ) : TypedFunction<CladogeneticSpeciationRateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chromosomesPloidyCladoEventsBD* Func_chromosomesPloidyCladoEventsBD::clone( void ) const {
    
    return new Func_chromosomesPloidyCladoEventsBD( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix >* Func_chromosomesPloidyCladoEventsBD::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    int ns = (int)static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::ChromosomesPloidyCladogenicBirthDeathFunction* f = new RevBayesCore::ChromosomesPloidyCladogenicBirthDeathFunction( sr, ns );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_chromosomesPloidyCladoEventsBD::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "speciation_rates", ModelVector<RealPos>::getClassTypeSpec() , "The speciation rates for different cladogenetic event types.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "max_chromosomes", Natural::getClassTypeSpec(), "The maximum number of chromosomes.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chromosomesPloidyCladoEventsBD::getClassType(void)
{
    
    static std::string rev_type = "Func_chromosomesPloidyCladoEventsBD";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_chromosomesPloidyCladoEventsBD::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_chromosomesPloidyCladoEventsBD::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChromosomesPloidyCladoEventsBD";
    
    return f_name;
}


const TypeSpec& Func_chromosomesPloidyCladoEventsBD::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

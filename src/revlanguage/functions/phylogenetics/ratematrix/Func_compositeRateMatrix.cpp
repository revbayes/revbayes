//
//  Func_compositeRateMatrix.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/13/24.
//

#include "Func_compositeRateMatrix.h"

#include <cstddef>

#include "CompositeRateMatrixFunction.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "TypedDagNode.h"
#include "AbstractRateMatrix.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_compositeRateMatrix::Func_compositeRateMatrix( void ) : TypedFunction<RateGenerator>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_compositeRateMatrix* Func_compositeRateMatrix::clone( void ) const
{
    
    return new Func_compositeRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_compositeRateMatrix::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm1 = static_cast<const ModelVector<RateGenerator> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm2 = static_cast<const ModelVector<RateGenerator> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
//    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* sr = static_cast<const RateMatrix&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
//    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* cr = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
//    bool rescale = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getDagNode()->getValue();
    
    // sanity check
    size_t num_rm1 = rm1->getValue().size();
    size_t num_rm2 = rm2->getValue().size();
    size_t size_rm1 = rm1->getValue()[0].size();
    size_t size_rm2 = rm2->getValue()[0].size();
    
    if (num_rm1 != size_rm2) {
        std::stringstream ss;
        ss << "Rate matrix vector 1 size (" << num_rm1 << ") and rate matrix vector 2 element size (" << size_rm2 << ") do not match";
        throw RbException( ss.str() );
    }
    if (num_rm2 != size_rm1) {
        std::stringstream ss;
        ss << "Rate matrix vector 2 size (" << num_rm2 << ") and rate matrix vector 1 element size (" << size_rm1 << ") do not match";
        throw RbException( ss.str() );
    }
    
//    if ( sr->getValue().size() != cr->getValue().size() )
//    {
//        throw RbException( "switch_rates and clock_rates have different numbers of classes." );
//    }
    
//    rm = rm1->getValue()
    
    
//    RevBayesCore::JcRateMatrixFunction* f = new RevBayesCore::JcRateMatrixFunction(4);
    RevBayesCore::CompositeRateMatrixFunction* f = new RevBayesCore::CompositeRateMatrixFunction( rm1, rm2 );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_compositeRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "Q1"                    , ModelVector<RateGenerator>::getClassTypeSpec(), "The first vector of M NxN rate matrices", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "Q2"                    , ModelVector<RateGenerator>::getClassTypeSpec(), "The second vector for N MxM rate matrices", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//        argumentRules.push_back( new ArgumentRule( "switch_rates"         , RateMatrix::getClassTypeSpec(), "The class-switching rate matrix", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//        argumentRules.push_back( new ArgumentRule( "clock_rates"          , ModelVector<RealPos>::getClassTypeSpec(), "The rate multipliers per class", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//        argumentRules.push_back( new ArgumentRule( "rescaled", RlBoolean::getClassTypeSpec(), "Should the matrix be normalized?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_compositeRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_compositeRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_compositeRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_compositeRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCompositeRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_compositeRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

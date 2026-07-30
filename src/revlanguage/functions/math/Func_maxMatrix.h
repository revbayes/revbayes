#ifndef Func_maxMatrix_H
#define Func_maxMatrix_H

#include "RlTypedFunction.h"
#include "Real.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the maximum value function for MatrixReal.
     * This non-templated version handles MatrixReal objects, returning Real.
     */
    class Func_maxMatrix : public TypedFunction<Real> {
        
    public:
        Func_maxMatrix( void );
        
        // Basic utility functions
        Func_maxMatrix*                                 clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };

}


#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "GenericFunction.h"
#include "RbConstants.h"
#include "MatrixReal.h"


namespace RevBayesCore {
    double* max_matrix(const MatrixReal& m)
    {
        double result = RbConstants::Double::neginf;
        for (size_t row = 0; row < m.size(); row++)
        {
            for (size_t col = 0; col < m[row].size(); col++)
            {
                if (m[row][col] > result)
                {
                    result = m[row][col];
                }
            }
        }
        return new double(result);
    }
}

using namespace RevLanguage;

 
/** default constructor */
Func_maxMatrix::Func_maxMatrix( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_maxMatrix* Func_maxMatrix::clone( void ) const
{
    
    return new Func_maxMatrix( *this );
}


RevBayesCore::TypedFunction<double>* Func_maxMatrix::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >* arg_matrix =
        static_cast< const MatrixReal& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    return RevBayesCore::generic_function_ptr<double>(RevBayesCore::max_matrix, arg_matrix);
}


/* Get argument rules */
const ArgumentRules& Func_maxMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", MatrixReal::getClassTypeSpec(), "A vector/matrix of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_maxMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_maxMatrix";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_maxMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_maxMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "max";
    
    return f_name;
}


const TypeSpec& Func_maxMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif

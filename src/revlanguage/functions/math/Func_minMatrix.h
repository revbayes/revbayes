#ifndef Func_minMatrix_H
#define Func_minMatrix_H

#include "RlTypedFunction.h"
#include "Real.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the minimum value function for MatrixReal.
     * This non-templated version handles MatrixReal and ModelVector<Real>, returning Real.
     *
     * @copybrief RevBayesCore::MinFunction
     * @see RevBayesCore::MinFunction for the internal object
     */
    class Func_minMatrix : public TypedFunction<Real> {
        
    public:
        Func_minMatrix( void );
        
        // Basic utility functions
        Func_minMatrix*                                 clone(void) const;                                          //!< Clone the object
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
#include "MinFunction.h"
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


namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;


/** default constructor */
Func_minMatrix::Func_minMatrix( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_minMatrix* Func_minMatrix::clone( void ) const
{
    
    return new Func_minMatrix( *this );
}


RevBayesCore::TypedFunction<double>* Func_minMatrix::createFunction( void ) const
{
    RevBayesCore::DagNode* node = this->args[0].getVariable()->getRevObject().getDagNode();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg_vector = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*>(node);
    RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >* arg_matrix = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >*>(node);
    
    RevBayesCore::MinFunction* f = arg_vector != NULL ? new RevBayesCore::MinFunction( arg_vector ) : new RevBayesCore::MinFunction( arg_matrix );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_minMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        paramTypes.push_back( MatrixReal::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "x", paramTypes, "A vector/matrix of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_minMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_minMatrix";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_minMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_minMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "min";
    
    return f_name;
}


const TypeSpec& Func_minMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif

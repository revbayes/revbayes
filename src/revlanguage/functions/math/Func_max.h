#ifndef Func_max_H
#define Func_max_H

#include "RlTypedFunction.h"
#include "Real.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the maximum value function.
     * This non-templated version handles MatrixReal and ModelVector<Real>, returning Real.
     *
     * @copybrief RevBayesCore::MaxFunction
     * @see RevBayesCore::MaxFunction for the internal object
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

    /**
     * The RevLanguage wrapper of the maximum value function for ModelVector types.
     *
     * @copybrief RevBayesCore::MaxFunction
     * @see RevBayesCore::MaxFunction for the internal object
     */
    template <typename valType, typename retType>
    class Func_max : public TypedFunction<retType> {
    
    public:
        Func_max( void );
    
        // Basic utility functions
        Func_max*                                                      clone(void) const;                           //!< Clone the object
        static const std::string&                                      getClassType(void);                          //!< Get Rev type
        static const TypeSpec&                                         getClassTypeSpec(void);                      //!< Get class type spec
        std::string                                                    getFunctionName(void) const;                 //!< Get the primary name of the function in Rev
        const TypeSpec&                                                getTypeSpec(void) const;                     //!< Get the type spec of the instance
    
        // Function functions you have to override
        RevBayesCore::TypedFunction< typename retType::valueType >*    createFunction(void) const;                  //!< Create internal function object
        const ArgumentRules&                                           getArgumentRules(void) const;                //!< Get argument rules
    
    };

}


#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "GenericFunction.h"
#include "Integer.h"
#include "IntegerPos.h"
#include "MaxFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RealPos.h"
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

 
/**
 * ===========================================
 * Non-templated (MatrixReal) Func_max version
 * ===========================================
 */


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
    RevBayesCore::DagNode* node = this->args[0].getVariable()->getRevObject().getDagNode();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg_vector = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*>(node);
    RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >* arg_matrix = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >*>(node);
    
    RevBayesCore::MaxFunction* f = arg_vector != NULL ? new RevBayesCore::MaxFunction( arg_vector ) : new RevBayesCore::MaxFunction( arg_matrix );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_maxMatrix::getArgumentRules( void ) const
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


/**
 * ========================================
 * Templated (ModelVector) Func_max version
 * ========================================
 */


namespace RevBayesCore {
    template <typename T>
    T* max_value(const RbVector<T>& v)
    {
        T m = v[0];
        for (size_t i = 1; i < v.size(); ++i)
        {
            if (v[i] > m)
            {
                m = v[i];
            }
        }
        return new T(m);
    }
}


template <typename valType, typename retType>
RevLanguage::Func_max<valType, retType>::Func_max( void ) : TypedFunction<retType>( )
{
    
}


template <typename valType, typename retType>
RevLanguage::Func_max<valType, retType>* RevLanguage::Func_max<valType, retType>::clone( void ) const
{
    
    return new Func_max<valType, retType>( *this );
}


template <typename valType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::Func_max<valType, retType>::createFunction( void ) const
{
    // Use generic_function_ptr to handle arbitrary numeric types. valType is the element type (e.g., Integer, RealPos),
    // and ModelVector<valType> wraps RbVector<valType::valueType>. So the element type is valType::valueType.
    typedef typename valType::valueType elementType;
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<elementType> >* arg_vector =
        static_cast< const ModelVector<valType>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // Use generic_function_ptr with the return type, which matches elementType
    return RevBayesCore::generic_function_ptr<typename retType::valueType>(RevBayesCore::max_value<elementType>, arg_vector);
}


template <typename valType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func_max<valType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        // Only accept ModelVector<valType>: MatrixReal is handled by the non-templated Func_max
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<valType>::getClassTypeSpec(), "A vector/matrix of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType, typename retType>
const std::string& RevLanguage::Func_max<valType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func_max<" + valType::getClassType() + "," + retType::getClassType() + ">";
    
    return rev_type;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_max<valType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


template <typename valType, typename retType>
std::string RevLanguage::Func_max<valType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "max";
    
    return f_name;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_max<valType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif

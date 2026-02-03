#ifndef Func_min_H
#define Func_min_H

#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the minimum value function for ModelVector types.
     */
    template <typename valType, typename retType>
    class Func_min : public TypedFunction<retType> {
        
    public:
        Func_min( void );
        
        // Basic utility functions
        Func_min*                                                      clone(void) const;                           //!< Clone the object
        static const std::string&                                      getClassType(void);                          //!< Get Rev type
        static const TypeSpec&                                         getClassTypeSpec(void);                      //!< Get class type spec
        std::string                                                    getFunctionName(void) const;                 //!< Get the primary name of the function in Rev
        const TypeSpec&                                                getTypeSpec(void) const;                     //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< typename retType::valueType >*    createFunction(void) const;                  //!< Create internal function object
        const ArgumentRules&                                           getArgumentRules(void) const;                //!< Get argument rules
        
    };
    
}


#include "GenericFunction.h"
#include "Integer.h"
#include "IntegerPos.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"


namespace RevBayesCore {
    template <typename T>
    T* min_value(const RbVector<T>& v)
    {
        T m = v[0];
        for (size_t i = 1; i < v.size(); ++i)
        {
            if (v[i] < m)
            {
                m = v[i];
            }
        }
        return new T(m);
    }
}


template <typename valType, typename retType>
RevLanguage::Func_min<valType, retType>::Func_min( void ) : TypedFunction<retType>( )
{
    
}


template <typename valType, typename retType>
RevLanguage::Func_min<valType, retType>* RevLanguage::Func_min<valType, retType>::clone( void ) const
{
    
    return new Func_min<valType, retType>( *this );
}


template <typename valType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::Func_min<valType, retType>::createFunction( void ) const
{
    // Use generic_function_ptr to handle arbitrary numeric types. valType is the element type (e.g., Integer, RealPos),
    // and ModelVector<valType> wraps RbVector<valType::valueType>. So the element type is valType::valueType.
    typedef typename valType::valueType elementType;
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<elementType> >* arg_vector = 
        static_cast< const ModelVector<valType>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // Use generic_function_ptr with the return type, which matches elementType
    return RevBayesCore::generic_function_ptr<typename retType::valueType>(RevBayesCore::min_value<elementType>, arg_vector);
}


template <typename valType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func_min<valType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        // Only accept ModelVector<valType>: MatrixReal is handled by the non-templated Func_min
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<valType>::getClassTypeSpec(), "A vector/matrix of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType, typename retType>
const std::string& RevLanguage::Func_min<valType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func_min<" + valType::getClassType() + "," + retType::getClassType() + ">";
    
    return rev_type;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_min<valType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


template <typename valType, typename retType>
std::string RevLanguage::Func_min<valType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "min";
    
    return f_name;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_min<valType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif

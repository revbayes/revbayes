#ifndef Func_median_H
#define Func_median_H

#include "RlTypedFunction.h"

#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the arithmetic median function.
     *
     * This function computes the arithmetic median of a vector of real number:
     *   y = x[ (x.size()-1)/2 ]                             if x.size() is even
     *   y = (x[ (x.size()-1)/2 ]+x[ x.size()/2 ]) / 2       otherwise
     */
    template <typename valType, typename retType>
    class Func_median : public TypedFunction<retType> {
        
    public:
        Func_median( void );
        
        // Basic utility functions
        Func_median*                                                   clone(void) const;                           //!< Clone the object
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
#include "RbConstants.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"


namespace RevBayesCore {
    // Compute median as double first to preserve fractional results, then convert to ReturnType
    template <typename ReturnType, typename InputType>
    ReturnType* median_value(const RbVector<InputType>& v)
    {
        double m = 0.0;
        
        // get a copy of the original vector
        RbVector<InputType> w = v;
        
        if ( w.size() == 0 )
        {
            m = RbConstants::Double::nan;
        }
        else if ( w.size() == 1 )
        {
            m = double(w[0]);
        }
        else
        {
            // we need to sort the vector first
            w.sort();

            // do we have an even number of elements
            if (w.size() % 2 == 0)
            {
                size_t index = w.size() / 2;
                // Convert to double before division to preserve fractional results
                m = (double(w[index-1]) + double(w[index])) / 2.0;
            }
            else
            {
                size_t index = floor(w.size() / 2);
                m = double(w[index]);
            }
        }
        
        // Convert the double result to ReturnType
        return new ReturnType(ReturnType(m));
    }
}


template <typename valType, typename retType>
RevLanguage::Func_median<valType, retType>::Func_median( void ) : TypedFunction<retType>( )
{
    
}


template <typename valType, typename retType>
RevLanguage::Func_median<valType, retType>* RevLanguage::Func_median<valType, retType>::clone( void ) const
{
    
    return new Func_median<valType, retType>( *this );
}


template <typename valType, typename retType>
RevBayesCore::TypedFunction< typename retType::valueType >* RevLanguage::Func_median<valType, retType>::createFunction( void ) const
{
    // Use generic_function_ptr to handle arbitrary numeric types. valType is the element type (e.g., Integer, RealPos),
    // and ModelVector<valType> wraps RbVector<valType::valueType>. So the element type is valType::valueType.
    typedef typename valType::valueType elementType;
    typedef typename retType::valueType returnType;
    
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<elementType> >* arg_vector =
        static_cast< const ModelVector<valType>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // Use a wrapper function that converts from elementType to returnType
    return RevBayesCore::generic_function_ptr<returnType>(RevBayesCore::median_value<returnType, elementType>, arg_vector);
}


template <typename valType, typename retType>
const RevLanguage::ArgumentRules& RevLanguage::Func_median<valType, retType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<valType>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType, typename retType>
const std::string& RevLanguage::Func_median<valType, retType>::getClassType(void)
{
    
    static std::string rev_type = "Func_median<" + valType::getClassType() + "," + retType::getClassType() + ">";
    
    return rev_type;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_median<valType, retType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


template <typename valType, typename retType>
std::string RevLanguage::Func_median<valType, retType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "median";
    
    return f_name;
}


template <typename valType, typename retType>
const RevLanguage::TypeSpec& RevLanguage::Func_median<valType, retType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif


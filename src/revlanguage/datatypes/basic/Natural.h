#ifndef Natural_H
#define Natural_H

#include "RbBoolean.h"
#include "Integer.h"

#include <ostream>
#include <string>
#include <type_traits>  // For std::enable_if and std::is_integral



namespace RevLanguage {

    class RealPos;

    /**
     * Primitive type for Natural numbers (positive integers including 0).
     *
     * Note that we derive this from Integer. To make
     * sure inheritance is safe, we restrict the range
     * of natural numbers from 0 to to INT_MAX
     */
    class Natural : public Integer {

        public:
        Natural(void);                                                                                      //!< Default constructor (value is 0)
        Natural(RevBayesCore::TypedDagNode<long> *v);   

        template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
            explicit Natural(T x) : Integer(x) {
            if (x < 0) {
                throw RbException("Negative value for " + getClassType());
            }
        }

        explicit Natural(const RevBayesCore::Boolean& b) : Integer(b ? 1 : 0) {}
      
        // Basic operator functions
        RevObject*                  add(const RevObject &rhs) const;                                        //!< Addition operator used for example in '+=' statements
        Natural*                    add(const Natural &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        RealPos*                    add(const RealPos &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        RevObject*                  divide(const RevObject &rhs) const;                                     //!< Division operator used for example in '/=' statements
        RealPos*                    divide(const Natural &rhs) const;                                       //!< Division operator used for example in '/=' statements
        RealPos*                    divide(const RealPos &rhs) const;                                       //!< Division operator used for example in '/=' statements
        RevObject*                  multiply(const RevObject &rhs) const;                                   //!< Multiplication operator used for example in '*=' statements
        Natural*                    multiply(const Natural &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements
        RealPos*                    multiply(const RealPos &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements

        // Basic utility functions
        Natural*                    clone(void) const;                                                      //!< Clone object
        RevObject*                  convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&   getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&      getClassTypeSpec(void);                                                 //!< Get class type spec
        const TypeSpec&             getTypeSpec(void) const;                                                //!< Get language type of the object
        double                      isConvertibleTo(const TypeSpec& type, bool once) const;                 //!< Is convertible to type?
        
        std::string                 getGuiName(void) { return "Natural"; }
        std::string                 getGuiUnicodeSymbol(void) { return "N"; }
        std::string                 getGuiInfo(void) { return ""; }
    };
    
}

#endif


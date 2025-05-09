#ifndef Integer_H
#define Integer_H

#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>
#include <cstdint>

namespace RevLanguage {
    
    class Real;

    class Integer : public ModelObject<std::int64_t> {

    public:
        Integer(void);                                                                                          //!< Default constructor
        Integer(RevBayesCore::TypedDagNode<std::int64_t> *v);                                                   //!< Constructor from DAG node
        Integer(std::int64_t v);                                                                                //!< Constructor from std::int64_t
        template <typename T, typename U=std::enable_if_t<std::is_integral_v<T>,T>> Integer(T v):               //!< Constructor from T
           Integer(std::int64_t(v)) {}

        // Basic operator functions
        virtual RevObject*              add(const RevObject &rhs) const;                                        //!< Addition operator used for example in '+=' statements
        Integer*                        add(const Integer &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        Real*                           add(const Real &rhs) const;                                             //!< Addition operator used for example in '+=' statements
        void                            decrement(void);                                                        //!< Decrement operator used for example in 'a--' statements
        virtual RevObject*              divide(const RevObject &rhs) const;                                     //!< Division operator used for example in '/=' statements
        Real*                           divide(const Integer &rhs) const;                                       //!< Division operator used for example in '/=' statements
        Real*                           divide(const Real &rhs) const;                                          //!< Division operator used for example in '/=' statements
        void                            increment(void);                                                        //!< Increment operator used for example in 'a++' statements
        virtual RevObject*              multiply(const RevObject &rhs) const;                                   //!< Multiplication operator used for example in '*=' statements
        Integer*                        multiply(const Integer &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements
        Real*                           multiply(const Real &rhs) const;                                        //!< Multiplication operator used for example in '*=' statements
        virtual RevObject*              subtract(const RevObject &rhs) const;                                   //!< Subtraction operator used for example in '-=' statements
        Integer*                        subtract(const Integer &rhs) const;                                     //!< Subtraction operator used for example in '-=' statements
        Real*                           subtract(const Real &rhs) const;                                        //!< Subtraction operator used for example in '-=' statements

        // Basic utility functions
        virtual Integer*                clone(void) const;                                                      //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        virtual const TypeSpec&         getTypeSpec(void) const;                                                //!< Get language type of the object
        virtual double                  isConvertibleTo(const TypeSpec& type, bool convert_by_value) const;                 //!< Is convertible to type?
        
        std::string                     getGuiName(void) { return "Integer"; }
        std::string                     getGuiUnicodeSymbol(void) { return "Z"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
    
}

#endif


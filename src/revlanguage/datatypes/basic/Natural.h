#ifndef Natural_H
#define Natural_H

#include "Integer.h"

#include <ostream>
#include <string>


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
        Natural(RevBayesCore::TypedDagNode<std::int64_t> *v);                                                       //!< Constructor with DAG node
        Natural(std::int64_t x);                                                                                    //!< Constructor from int
//        Natural(std::uint64_t x);                                                                           //!< Constructor from size_t

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


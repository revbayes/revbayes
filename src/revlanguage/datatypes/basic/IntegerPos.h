#ifndef IntegerPos_H
#define IntegerPos_H

#include "Natural.h"

#include <ostream>
#include <string>



namespace RevLanguage {

    class RealPos;

    /**
     * Primitive type for strictly positive integer numbers (same as Natural without 0).
     *
     * Note that we derive this from Natural. To make
     * sure inheritance is safe, we restrict the range
     * of positive integers numbers from 1 to INT_MAX
     */
    class IntegerPos : public Natural {

        public:
        IntegerPos(void);                                                                                      //!< Default constructor (value is 0)
        IntegerPos(RevBayesCore::TypedDagNode<std::int64_t> *v);                                                       //!< Constructor with DAG node
        IntegerPos(std::int64_t x);                                                                                    //!< Constructor from int
//        Natural(std::uint64_t x);                                                                           //!< Constructor from size_t

        // Basic operator functions
        RevObject*                  add(const RevObject &rhs) const;                                        //!< Addition operator used for example in '+=' statements
        IntegerPos*                 add(const IntegerPos &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        RealPos*                    add(const RealPos &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        RevObject*                  divide(const RevObject &rhs) const;                                     //!< Division operator used for example in '/=' statements
        RealPos*                    divide(const IntegerPos &rhs) const;                                       //!< Division operator used for example in '/=' statements
        RealPos*                    divide(const RealPos &rhs) const;                                       //!< Division operator used for example in '/=' statements
        RevObject*                  multiply(const RevObject &rhs) const;                                   //!< Multiplication operator used for example in '*=' statements
        IntegerPos*                 multiply(const IntegerPos &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements
        RealPos*                    multiply(const RealPos &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements

        // Basic utility functions
        IntegerPos*                 clone(void) const;                                                      //!< Clone object
        RevObject*                  convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&   getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&      getClassTypeSpec(void);                                                 //!< Get class type spec
        const TypeSpec&             getTypeSpec(void) const;                                                //!< Get language type of the object
        double                      isConvertibleTo(const TypeSpec& type, bool once) const;                 //!< Is convertible to type?
        
        std::string                 getGuiName(void) { return "IntegerPos"; }
        std::string                 getGuiUnicodeSymbol(void) { return "+Z"; }
        std::string                 getGuiInfo(void) { return ""; }
    };
    
}

#endif


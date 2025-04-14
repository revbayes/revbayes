#ifndef Real_H
#define Real_H

#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

/**
 * Real is the class used to hold a real value. Internally, the real is represented by
 * a double.
 *
 * @note Some functions are virtual because RealPos is derived from Real
 */

namespace RevLanguage {
    
    class Integer;

    class Real : public ModelObject<double> {

        public:

        Real(void);                                                                                             //!< Default constructor (0.0)
        Real(RevBayesCore::TypedDagNode<double> *v);                                                            //!< Construct from DAG node
        Real(double v);                                                                                         //!< Construct from double
        Real(std::int64_t v);                                                                                           //!< Construct from int
        Real(bool v);                                                                                           //!< Construct from bool
        Real(const Real& x);                                                                                    //!< Copy constructor
        
        // Basic operator functions
        virtual RevObject*              add(const RevObject &rhs) const;                                        //!< Addition operator used for example in '+=' statements
        Real*                           add(const Real &rhs) const;                                             //!< Addition operator used for example in '+=' statements
        Real*                           add(const Integer &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        void                            decrement(void);                                                        //!< Decrement operator used for example in 'a--' statements
        virtual RevObject*              divide(const RevObject &rhs) const;                                     //!< Division operator used for example in '/=' statements
        Real*                           divide(const Real &rhs) const;                                          //!< Division operator used for example in '/=' statements
        Real*                           divide(const Integer &rhs) const;                                       //!< Division operator used for example in '/=' statements
        void                            increment(void);                                                        //!< Decrement operator used for example in 'a++' statements
        virtual RevObject*              multiply(const RevObject &rhs) const;                                   //!< Multiplication operator used for example in '*=' statements
        Real*                           multiply(const Real &rhs) const;                                        //!< Multiplication operator used for example in '*=' statements
        Real*                           multiply(const Integer &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements
        virtual RevObject*              subtract(const RevObject &rhs) const;                                   //!< Subtraction operator used for example in '-=' statements
        Real*                           subtract(const Real &rhs) const;                                        //!< Subtraction operator used for example in '-=' statements
        Real*                           subtract(const Integer &rhs) const;                                     //!< Subtraction operator used for example in '-=' statements

        // Basic utility functions
        virtual Real*                   clone(void) const;                                                      //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        virtual const TypeSpec&         getTypeSpec(void) const;                                                //!< Get language type of the object
        virtual double                  isConvertibleTo(const TypeSpec& type, bool once) const;                 //!< Is convertible to type?
    
        std::string                     getGuiName(void) { return "Real"; }
        std::string                     getGuiUnicodeSymbol(void) { return "R"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
    
}

#endif


/**
 * @file
 * This file contains the declaration of RealPos, which
 * is the primitive RevBayes type for positive real numbers.
 *
 * @brief Declaration of RealPos
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef RealPos_H
#define RealPos_H

#include "Real.h"

#include <ostream>
#include <string>

namespace RevLanguage {

    class Natural;
    
    class RealPos : public Real {

        public:
        RealPos(void);                                                              //!< Default constructor
        RealPos(RevBayesCore::TypedDagNode<double> *x);                             //!< Construct from double
        RealPos(double x);                                                          //!< Construct from double
        RealPos(long x);                                                             //!< Construct from int
        RealPos(bool x);                                                             //!< Construct from bool

        // Basic operator functions
        virtual RevObject*              add(const RevObject &rhs) const;            //!< Addition operator used for example in '+=' statements
        RealPos*                        add(const Natural &rhs) const;              //!< Addition operator used for example in '+=' statements
        RealPos*                        add(const RealPos &rhs) const;              //!< Addition operator used for example in '+=' statements
        virtual RevObject*              divide(const RevObject &rhs) const;         //!< Division operator used for example in '/=' statements
        RealPos*                        divide(const Natural &rhs) const;           //!< Division operator used for example in '/=' statements
        RealPos*                        divide(const RealPos &rhs) const;           //!< Division operator used for example in '/=' statements
        virtual RevObject*              multiply(const RevObject &rhs) const;       //!< Multiplication operator used for example in '*=' statements
        RealPos*                        multiply(const Natural &rhs) const;         //!< Multiplication operator used for example in '*=' statements
        RealPos*                        multiply(const RealPos &rhs) const;         //!< Multiplication operator used for example in '*=' statements

        // Basic utility functions
        virtual RealPos*                clone(void) const;                          //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                         //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                     //!< Get class type spec
        virtual const TypeSpec&         getTypeSpec(void) const;                    //!< Get language type of the object
        virtual double                  isConvertibleTo(const TypeSpec& type, bool once) const;
    
        std::string                     getGuiName(void) { return "Positive Real"; }
        std::string                     getGuiUnicodeSymbol(void) { return "+R"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
    
}

#endif


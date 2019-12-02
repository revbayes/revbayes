/**
 * @file
 * This file contains the declaration of Probability, which
 * is the primitive RevBayes type for positive real numbers between 0 and 1.
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

#ifndef Probability_H
#define Probability_H

#include "RealPos.h"

#include <ostream>
#include <string>

namespace RevLanguage {

    class Probability : public RealPos {
    
    public:
        Probability(void);                                                                                      //!< Default constructor
        Probability(RevBayesCore::TypedDagNode<double> *x);                                                     //!< Construct from double
        Probability(double x);                                                                                  //!< Construct from double
        
        // Basic utility functions
        Probability*                    clone(void) const;                                                      //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                                                //!< Get language type of the object
        virtual double                  isConvertibleTo(const TypeSpec& type, bool once) const;                 //!< Is convertible to type?

        std::string                     getGuiName(void) { return "Probability"; }
        std::string                     getGuiUnicodeSymbol(void) { return "P"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
    
}

#endif


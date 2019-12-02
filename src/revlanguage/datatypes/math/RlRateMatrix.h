/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the RateMatrix class.
 *
 * @brief Declaration of RlRateMatrix
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-06 20:14:22 +0200 (Mon, 06 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: Real.h 1746 2012-08-06 18:14:22Z hoehna $
 */

#ifndef RlRateMatrix_H
#define RlRateMatrix_H

#include "RlRateGenerator.h"
#include "RateMatrix.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class RateMatrix : public RateGenerator {
        
    public:
        
        RateMatrix(void);                                                                                                           //!< Default constructor
        RateMatrix(const RevBayesCore::RateMatrix& m);                                                                              //!< Default constructor
        RateMatrix(RevBayesCore::RateMatrix *m);                                                                                    //!< Default constructor
        RateMatrix(RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator> *d);                                                                                                        //!< Default constructor
        
        // Basic utility functions
        RateMatrix*                         clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
    protected:
        void                                initMethods(void);
    };
    
}

#endif

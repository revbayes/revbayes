/**
 * @file
 * This file contains the declaration of a RootedTripletDistribution, which is
 * class that holds a rooted triplet distribution with their frequencies.
 *
 * @brief Declaration of RootedTripletDistribution
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef RlRootedTripletDistribution_H
#define RlRootedTripletDistribution_H

#include "ModelObject.h"
#include "RootedTripletDistribution.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    class RootedTripletDistribution : public ModelObject<RevBayesCore::RootedTripletDistribution> {
        
    public:
        RootedTripletDistribution(void);                                                                                                     //!< Constructor requires character type
        RootedTripletDistribution(RevBayesCore::RootedTripletDistribution *v);                                                                                //!< Constructor requires character type
        RootedTripletDistribution(const RevBayesCore::RootedTripletDistribution &v);                                                                          //!< Constructor requires character type
        RootedTripletDistribution(RevBayesCore::TypedDagNode<RevBayesCore::RootedTripletDistribution> *n);                                                    //!< Constructor requires character type
        
        typedef RevBayesCore::RootedTripletDistribution valueType;
        
        // Basic utility functions
        RootedTripletDistribution*          clone(void) const;                                                              //!< Clone object
        static const std::string&           getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                         //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found); //!< Map member methods to internal functions
        
        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }        
        std::string                         getGuiInfo(void) { return ""; }
    };
    
}

#endif

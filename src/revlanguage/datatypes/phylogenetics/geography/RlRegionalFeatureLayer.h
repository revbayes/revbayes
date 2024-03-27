//
//  RlRegionalFeatureLayer.hpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#ifndef RlRegionalFeatureLayer_hpp
#define RlRegionalFeatureLayer_hpp

#include "ModelObject.h"
#include "RegionalFeatureLayer.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class RlRegionalFeatureLayer : public ModelObject<RevBayesCore::RegionalFeatureLayer> {
        
    public:
        RlRegionalFeatureLayer(void);                                                                          //!< Default constructor
        RlRegionalFeatureLayer(RevBayesCore::RegionalFeatureLayer *m);                                                    //!< Default constructor
        RlRegionalFeatureLayer(const RevBayesCore::RegionalFeatureLayer &m);                                       //!< Const ref constructor
        RlRegionalFeatureLayer(RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatureLayer> *d);                        //!< Default constructor
        typedef RevBayesCore::RegionalFeatureLayer valueType;
        
        
        // Basic utility functions
        RlRegionalFeatureLayer*             clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions

        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
        
    private:
        void initMethods(void);
//        RevBayesCore::RegionalFeatureLayer*            layer;
    };
    
}

#endif /* RlRegionalFeatureLayer_hpp */

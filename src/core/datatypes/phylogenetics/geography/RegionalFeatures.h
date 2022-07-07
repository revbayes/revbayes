//
//  RegionalFeatures.h
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__RegionalFeatures__
#define __rb_mlandis__RegionalFeatures__

#include <stddef.h>
#include <iosfwd>
#include <vector>

#include "Cloneable.h"

namespace RevBayesCore {
// class GeographicArea;
// class RegionalFeaturesDataReader;

    class RegionalFeatures : public Cloneable {
        
    public:
        RegionalFeatures(std::vector<std::string> fn, std::vector<std::string> fr, std::vector<std::string> ft);
        RegionalFeatures(const RegionalFeatures& a);
        RegionalFeatures&                                      operator=(const RegionalFeatures& a);
        virtual RegionalFeatures*                              clone(void) const;
        
        size_t                                                  getNumLayers(void) const;
        
    
    private:
        std::vector<std::string> featureNames;
        std::vector<std::string> featureRelationships;
        std::vector<std::string> featureTypes;
        // TODO: replace with std::vector< RegionalFeatureLayer >
//        std::vector<RegionalFeatureLayer> featureLayers;
        
        unsigned                                        numLayers;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RegionalFeatures& x);
}


#endif /* defined(__rb_mlandis__RegionalFeatures__) */

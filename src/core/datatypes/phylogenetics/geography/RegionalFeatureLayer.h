//
//  RegionalFeatureLayer.hpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#ifndef RegionalFeatureLayer_hpp
#define RegionalFeatureLayer_hpp

#include "Cloneable.h"
#include <vector>

namespace RevBayesCore {

    class RegionalFeatureLayer : public Cloneable
    {
    public:
        
        RegionalFeatureLayer();
        RegionalFeatureLayer(size_t ti, size_t fi, std::string ft, std::string fr);
        RegionalFeatureLayer(const RegionalFeatureLayer& g);
        RegionalFeatureLayer&               operator=(const RegionalFeatureLayer& g);
        virtual RegionalFeatureLayer*       clone(void) const;
        void                                setFeatures(std::vector<int> f);
        void                                setFeatures(std::vector<double> f);
        void                                setFeatures(std::vector<std::vector<int> > f);
        void                                setFeatures(std::vector<std::vector<double> > f);
        std::vector<int>                    getWithinCategoricalFeatures(void);
        std::vector<double>                 getWithinQuantitativeFeatures(void);
        std::vector<std::vector<int> >      getBetweenCategoricalFeatures(void);
        std::vector<std::vector<double> >   getBetweenQuantitativeFeatures(void);
    
    private:
        size_t                              time_index;
        size_t                              feature_index;
        std::string                         feature_type;
        std::string                         feature_relationship;
        
        // this could be a template
        std::vector<int>                    within_categorical;
        std::vector<double>                 within_quantitative;
        std::vector<std::vector<int> >      between_categorical;
        std::vector<std::vector<double> >   between_quantitative;
    };

    std::ostream& operator<<(std::ostream& o, const RegionalFeatureLayer& x);
}

#endif /* RegionalFeatureLayer_hpp */

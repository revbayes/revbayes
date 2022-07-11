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
#include <map>

#include "Cloneable.h"
#include "RegionalFeatureLayer.h"



namespace RevBayesCore {

    class RegionalFeatures : public Cloneable {
        
    public:
        RegionalFeatures(void);
        RegionalFeatures(std::map<size_t, std::map<size_t, std::vector<long> > > wc,
                         std::map<size_t, std::map<size_t, std::vector<double> > > wq,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<long> > > > bc,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > bq);
        RegionalFeatures(const RegionalFeatures& a);
        RegionalFeatures&                                       operator=(const RegionalFeatures& a);
        virtual RegionalFeatures*                               clone(void) const;
        
        const std::vector<std::vector<RegionalFeatureLayer> >&  getLayers(std::string feature_relationship, std::string feature_type);
        const std::vector<RegionalFeatureLayer>&                getLayers(std::string feature_relationship, std::string feature_type, size_t time_index);
        const RegionalFeatureLayer&                             getLayers(std::string feature_relationship, std::string feature_type, size_t time_index, size_t feature_index);
        
        size_t                                          getNumLayers(void) const;
        
    
    private:
        void initializeFeatures();
        void normalizeWithinQuantitative();
        void normalizeBetweenQuantitative();
        
        std::map<size_t, std::map<size_t, std::vector<long> > > withinCategorical;
        std::map<size_t, std::map<size_t, std::vector<double> > > withinQuantitative;
        std::map<size_t, std::map<size_t, std::vector<std::vector<long> > > > betweenCategorical;
        std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > betweenQuantitative;
        
        // relationship, type, time_index, feature_index
        std::map<std::string, std::map<std::string, std::vector< std::vector< RegionalFeatureLayer > > > > feature_layers;
        
        unsigned                                        numLayers;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RegionalFeatures& x);
}


#endif /* defined(__rb_mlandis__RegionalFeatures__) */

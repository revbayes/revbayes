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



namespace RevBayesCore {
// class GeographicArea;
// class RegionalFeaturesDataReader;


    struct RegionalFeatureLayer {
        RegionalFeatureLayer() :
            time_index(0), feature_index(0), feature_type(""), feature_relationship("")
        { }
        RegionalFeatureLayer(size_t ti, size_t fi, std::string ft, std::string fr) :
            time_index(ti), feature_index(fi), feature_type(ft), feature_relationship(fr)
        { }
        void setFeature(std::vector<int> f) { within_categorical = f; }
        void setFeature(std::vector<double> f) { within_quantitative = f; }
        void setFeature(std::vector<std::vector<int> > f) { between_categorical = f; }
        void setFeature(std::vector<std::vector<double> > f) { between_quantitative = f; }
        
        size_t                              time_index;
        size_t                              feature_index;
        std::string                         feature_type;
        std::string                         feature_relationship;
        std::vector<int>                    within_categorical;
        std::vector<double>                 within_quantitative;
        std::vector<std::vector<int> >      between_categorical;
        std::vector<std::vector<double> >   between_quantitative;
    };

    class RegionalFeatures : public Cloneable {
        
    public:
        RegionalFeatures(std::map<size_t, std::map<size_t, std::vector<int> > > wc,
                         std::map<size_t, std::map<size_t, std::vector<double> > > wq,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<int> > > > bc,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > bq);
        RegionalFeatures(const RegionalFeatures& a);
        RegionalFeatures&                               operator=(const RegionalFeatures& a);
        virtual RegionalFeatures*                       clone(void) const;
        
        const std::vector<std::vector<RegionalFeatureLayer> >&  getLayers(std::string feature_relationship, std::string feature_type);
        const std::vector<RegionalFeatureLayer>&                getLayers(std::string feature_relationship, std::string feature_type, size_t time_index);
        const RegionalFeatureLayer&                             getLayers(std::string feature_relationship, std::string feature_type, size_t time_index, size_t feature_index);
        
        size_t                                          getNumLayers(void) const;
        
    
    private:
        std::map<size_t, std::map<size_t, std::vector<int> > > withinCategorical;
        std::map<size_t, std::map<size_t, std::vector<double> > > withinQuantitative;
        std::map<size_t, std::map<size_t, std::vector<std::vector<int> > > > betweenCategorical;
        std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > betweenQuantitative;
        
        // relationship, type, time_index, feature_index
        std::map<std::string, std::map<std::string, std::vector< std::vector< RegionalFeatureLayer > > > > feature_layers;
        
        // TODO: replace with std::vector< RegionalFeatureLayer >
        // std::vector<RegionalFeatureLayer> featureLayers;
        
        unsigned                                        numLayers;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RegionalFeatures& x);
}


#endif /* defined(__rb_mlandis__RegionalFeatures__) */

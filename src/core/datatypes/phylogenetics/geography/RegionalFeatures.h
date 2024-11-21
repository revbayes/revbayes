//
//  RegionalFeatures.h
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__RegionalFeatures__
#define __rb_mlandis__RegionalFeatures__

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <map>
#include <cstdint>

#include "Cloneable.h"
#include "RegionalFeatureLayer.h"



namespace RevBayesCore {

    class RegionalFeatures : public Cloneable {
        
    public:
        RegionalFeatures(void);
        RegionalFeatures(std::map<size_t, std::map<size_t, std::vector<std::int64_t> > > wc,
                         std::map<size_t, std::map<size_t, std::vector<double> > > wq,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<std::int64_t> > > > bc,
                         std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > bq);
        RegionalFeatures(const RegionalFeatures& a);
        RegionalFeatures&                                       operator=(const RegionalFeatures& a);
        virtual RegionalFeatures*                               clone(void) const;
        
        const std::vector<std::vector<RegionalFeatureLayer> >&  getLayers(std::string feature_relationship, std::string feature_type);
        const std::vector<RegionalFeatureLayer>&                getLayers(std::string feature_relationship, std::string feature_type, size_t time_index);
        const RegionalFeatureLayer&                             getLayers(std::string feature_relationship, std::string feature_type, size_t time_index, size_t feature_index);
        
        
        void            normalizeWithinQuantitative();
        void            normalizeBetweenQuantitative();
        std::map<std::string, std::map<std::string, size_t> >  getNumLayers(void) const;
        size_t          getNumTimeslices(void) const;
        
    
    private:
        void initializeFeatures();
        
        // [time_index][feature_index][this_region]
        std::map<size_t, std::map<size_t, std::vector<std::int64_t> > > withinCategorical;
        std::map<size_t, std::map<size_t, std::vector<double> > > withinQuantitative;
        
        // example: [cat_or_quant][time_index][feature_index][this_region]
        // std::vector< std::vector< std::vector< std::vector<double> > > > withinFeatures;
        
        // [time_index][feature_index][from_region][to_region]
        std::map<size_t, std::map<size_t, std::vector<std::vector<std::int64_t> > > > betweenCategorical;
        std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > betweenQuantitative;
        
        // make everything look like this:
        // std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > >
        
        
        
        // relationship, type, time_index, feature_index
        std::map<std::string, std::map<std::string, std::vector< std::vector< RegionalFeatureLayer > > > > feature_layers;
        
        size_t                                                  numTimeslices;
        std::map<std::string, std::map<std::string, size_t> >   numLayers;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RegionalFeatures& x);
}


#endif /* defined(__rb_mlandis__RegionalFeatures__) */

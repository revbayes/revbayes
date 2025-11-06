//
//  TimeAtlasDataReader.h
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__TimeAtlasDataReader__
#define __rb_mlandis__TimeAtlasDataReader__

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "DelimitedDataReader.h"
#include "boost/property_tree/ptree_fwd.hpp"

namespace RevBayesCore {
class GeographicArea;

    class TimeAtlasDataReader : public DelimitedDataReader
    {
    public:
        
        TimeAtlasDataReader(const path& fn, std::string d="");
        TimeAtlasDataReader(const TimeAtlasDataReader& tadr);
        
        std::vector<double> getEpochs(void);
        std::vector<std::vector<GeographicArea*> > getAreas(void);
        void readJson(void);
        void printJson(boost::property_tree::ptree const& pt);
        void fillData(boost::property_tree::ptree const& pt);
        
        
    protected:
        
        void setEpochs(void);
        void setAreas(void);
        void sortEpochs(void);
        
        std::vector<double> epochs;
        std::vector<std::vector<GeographicArea*> > areas;
        
        struct AgeIndexPair
        {
            double age;
            size_t index;
            AgeIndexPair(double a, size_t idx) : age(a), index(idx) {};
        };
        
        struct AgeCompare
        {
            bool operator()(const AgeIndexPair first, const AgeIndexPair second)
            {
                return first.age > second.age;
            }
        };
        
    };
}

#endif /* defined(__rb_mlandis__TimeAtlasDataReader__) */

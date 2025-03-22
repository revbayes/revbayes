//
//  RegionalFeatures.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//


#include "RegionalFeatures.h"
#include "RegionalFeatureLayer.h"
#include "RbConstants.h"
#include "RbException.h"

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <string>

using namespace RevBayesCore;


RegionalFeatures::RegionalFeatures(void) {
    ;
}


RegionalFeatures::RegionalFeatures(std::map<size_t, std::map<size_t, std::vector<std::int64_t> > > wc,
                                   std::map<size_t, std::map<size_t, std::vector<double> > > wq,
                                   std::map<size_t, std::map<size_t, std::vector<std::vector<std::int64_t> > > > bc,
                                   std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > bq) :
    withinCategorical(wc),
    withinQuantitative(wq),
    betweenCategorical(bc),
    betweenQuantitative(bq)
{
    
    initializeFeatures();
}

RegionalFeatures::RegionalFeatures(const RegionalFeatures& a)
{
    if (this != &a) {
        numLayers = a.numLayers;
        numTimeslices = a.numTimeslices;
        withinCategorical = a.withinCategorical;
        withinQuantitative = a.withinQuantitative;
        betweenCategorical = a.betweenCategorical;
        betweenQuantitative = a.betweenQuantitative;
        feature_layers = a.feature_layers;
        *this = a;
    }
}

RegionalFeatures& RegionalFeatures::operator=(const RegionalFeatures& a)
{
    if (this != &a) {
        numLayers = a.numLayers;
        numTimeslices = a.numTimeslices;
        withinCategorical = a.withinCategorical;
        withinQuantitative = a.withinQuantitative;
        betweenCategorical = a.betweenCategorical;
        betweenQuantitative = a.betweenQuantitative;
        feature_layers = a.feature_layers;
    }
    return *this;
}

RegionalFeatures* RegionalFeatures::clone(void) const
{
    return new RegionalFeatures(*this);
}

void RegionalFeatures::initializeFeatures(void) {
    
    // how many timeslices
    numTimeslices = withinCategorical.size();
    
    // how many layers for each relationship/type
    numLayers["within"]["categorical"] = withinCategorical.begin()->second.size();
    numLayers["within"]["quantitative"] = withinQuantitative.begin()->second.size();
    numLayers["between"]["categorical"] = betweenCategorical.begin()->second.size();
    numLayers["between"]["quantitative"] = betweenQuantitative.begin()->second.size();
    
    for (auto it = withinCategorical.begin(); it != withinCategorical.end(); it++) {
        feature_layers["within"]["categorical"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t time_index = it->first;
        // loop over time-feature maps
        for (auto jt = withinCategorical[time_index].begin(); jt != withinCategorical[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            RegionalFeatureLayer tmp( time_index, feature_index, "within", "categorical" );
            std::vector<std::vector<double> > tmp_val;
            tmp_val.push_back(std::vector<double>());
            for (size_t i = 0; i < withinCategorical[time_index][feature_index].size(); i++) {
                tmp_val[0].push_back((double)withinCategorical[time_index][feature_index][i]);
            }
            tmp.setFeatures(tmp_val);
            feature_layers["within"]["categorical"][time_index-1].push_back(tmp);
        }
    }
    
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        feature_layers["within"]["quantitative"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t time_index = it->first;
        // loop over time-feature maps
        for (auto jt = withinQuantitative[time_index].begin(); jt != withinQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            RegionalFeatureLayer tmp( it->first, jt->first, "within", "quantitative" );
            std::vector<std::vector<double> > tmp_val;
            tmp_val.push_back( withinQuantitative[time_index][feature_index] );
            tmp.setFeatures(tmp_val);
            feature_layers["within"]["quantitative"][time_index-1].push_back(tmp);
        }
    }
    
    for (auto it = betweenCategorical.begin(); it != betweenCategorical.end(); it++) {
        feature_layers["between"]["categorical"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t time_index = it->first;
        // loop over time-feature maps
        for (auto jt = betweenCategorical[time_index].begin(); jt != betweenCategorical[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            RegionalFeatureLayer tmp(time_index, feature_index, "between", "categorical" );
            std::vector<std::vector<double> > tmp_val;
            for (size_t i = 0; i < betweenCategorical[time_index][feature_index].size(); i++) {
                tmp_val.push_back(std::vector<double>());
                for (size_t j = 0; j < betweenCategorical[time_index][feature_index][i].size(); j++) {
                    tmp_val[i].push_back((double)betweenCategorical[time_index][feature_index][i][j]);
                }
            }
            tmp.setFeatures(tmp_val);
            feature_layers["between"]["categorical"][time_index-1].push_back(tmp);
        }
    }
    
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        feature_layers["between"]["quantitative"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t time_index = it->first;
        // loop over time-feature maps
        for (auto jt = betweenQuantitative[time_index].begin(); jt != betweenQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            RegionalFeatureLayer tmp(time_index, feature_index, "between", "quantitative" );
            std::vector<std::vector<double> > tmp_val = betweenQuantitative[time_index][feature_index];
            tmp.setFeatures(tmp_val);
            feature_layers["between"]["quantitative"][time_index-1].push_back(tmp);
        }
    }
    
    // normalize all quantitative features across time slices
//    normalizeWithinQuantitative();
//    normalizeBetweenQuantitative();
    
}

void RegionalFeatures::normalizeWithinQuantitative(void) {
    
    size_t n_elem = 0;
    
    // compute sample mean
    double m = 0.0;
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = withinQuantitative[time_index].begin(); jt != withinQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["within"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v[0].size(); a++) {
                if (std::isnan( v[0][a] )) {
                    ; // do nothing
                } else {
                    m += v[0][a];
                    n_elem += 1;
                }
            }
        }
    }
    double mean = m / n_elem;
    
    // compute sample standard deviation
    double stddev = 0;
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = withinQuantitative[time_index].begin(); jt != withinQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["within"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v[0].size(); a++) {
                if (std::isnan( v[0][a] )) {
                    ; // do nothing
                } else {
                    stddev += std::pow(v[0][a] - m, 2);
                }
            }
        }
    }
    
    // do not use stddev if it equals zero
    stddev = std::sqrt( stddev / n_elem );
    if (stddev == 0.0) {
        throw RbException("RegionalFeatures::normalizeWithinQuantitative can only standardize data if stddev != 0 (i.e. features must contain variation.)");
    }
    
    // standardize all values
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = withinQuantitative[time_index].begin(); jt != withinQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["within"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v[0].size(); a++) {
                if (std::isnan( v[0][a] )) {
                    ; // do nothing
                } else {
                    v[0][a] = (v[0][a] - mean) / stddev;
                }
            }
            feature_layers["within"]["quantitative"][i][j].setFeatures(v);
        }
    }
    return;
}

void RegionalFeatures::normalizeBetweenQuantitative(void) {
    
    size_t n_elem = 0;
    
    // compute sample mean
    double m = 0.0;
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = betweenQuantitative[time_index].begin(); jt != betweenQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["between"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v.size(); a++) {
                for (size_t b = 0; b < v[a].size(); b++) {
                    if (std::isnan(v[a][b]) || a == b) {
                        ; // do nothing
                    } else {
                        m += v[a][b];
                        n_elem += 1;
                    }
                }
            }
        }
    }
    double mean = m / n_elem;
    
    // compute sample standard deviation
    double stddev = 0;
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = betweenQuantitative[time_index].begin(); jt != betweenQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["between"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v.size(); a++) {
                for (size_t b = 0; b < v[a].size(); b++) {
                    if (std::isnan(v[a][b]) || a == b) {
                        ; // do nothing
                    } else {
                        stddev += std::pow(v[a][b] - mean, 2);
                    }
                }
            }
        }
    }
    
    stddev = std::sqrt( stddev / n_elem );
    if (stddev == 0.0) {
        throw RbException("RegionalFeatures::normalizeBetweenQuantitative can only standardize data if stddev != 0 (i.e. features must contain variation.)");
    }
    
    
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        size_t time_index = it->first;
        size_t i = time_index - 1;
        for (auto jt = betweenQuantitative[time_index].begin(); jt != betweenQuantitative[time_index].end(); jt++) {
            size_t feature_index = jt->first;
            size_t j = feature_index - 1;
            std::vector<std::vector<double> > v = feature_layers["between"]["quantitative"][i][j].getFeatureValues();
            for (size_t a = 0; a < v.size(); a++) {
                for (size_t b = 0; b < v[a].size(); b++) {
                    if (std::isnan( v[a][b] ) || a == b) {
                        ; // do nothing
                    } else {
                        v[a][b] = (v[a][b] - mean) / stddev;
                    }
                }
            }
            feature_layers["between"]["quantitative"][i][j].setFeatures(v);
        }
    }
    
    return;
}

const std::vector<std::vector<RegionalFeatureLayer> >& RegionalFeatures::getLayers(std::string feature_relationship, std::string feature_type)
{
    return feature_layers[feature_relationship][feature_type];
}
const std::vector<RegionalFeatureLayer>& RegionalFeatures::getLayers(std::string feature_relationship, std::string feature_type, size_t time_index)
{
    
    return feature_layers[feature_relationship][feature_type][time_index];
}
const RegionalFeatureLayer& RegionalFeatures::getLayers(std::string feature_relationship, std::string feature_type, size_t time_index, size_t feature_index)
{
    return feature_layers[feature_relationship][feature_type][time_index][feature_index];
}

std::map<std::string, std::map<std::string, size_t> > RegionalFeatures::getNumLayers(void) const {
    return numLayers;
}

size_t RegionalFeatures::getNumTimeslices(void) const {
    return numTimeslices;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RegionalFeatures& x) {

    std::stringstream s;

    // Generate nice header
    o << std::endl;
    s << "Regional Features" << std::endl;
//    s << "RegionalFeatures with " << x.getAreas().size() << " epochs and " << x.getAreas().size() << " areas" << std::endl;
    o << s.str();

    for ( size_t i = 0; i < s.str().length() - 1; ++i )
        o << "=";
    o << std::endl;
    
    return o;
}


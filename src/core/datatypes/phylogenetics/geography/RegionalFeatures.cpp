//
//  RegionalFeatures.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "RegionalFeatures.h"
#include "RbConstants.h"

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <string>

using namespace RevBayesCore;


RegionalFeatures::RegionalFeatures(std::map<size_t, std::map<size_t, std::vector<int> > > wc,
                                   std::map<size_t, std::map<size_t, std::vector<double> > > wq,
                                   std::map<size_t, std::map<size_t, std::vector<std::vector<int> > > > bc,
                                   std::map<size_t, std::map<size_t, std::vector<std::vector<double> > > > bq) :
    withinCategorical(wc),
    withinQuantitative(wq),
    betweenCategorical(bc),
    betweenQuantitative(bq)
{
    
    
    // relationship, type, time_index, feature_index
    // feature_layers[ relationship ][ type ][ time_index ][ feature_index ]
    
    // loop over time maps
    for (auto it = withinCategorical.begin(); it != withinCategorical.end(); it++) {
        feature_layers["within"]["categorical"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t i = it->first - 1;
        // loop over time-feature maps
        for (auto jt = withinCategorical[it->first].begin(); jt != withinCategorical[it->first].end(); jt++) {
            RegionalFeatureLayer tmp( it->first, jt->first, "within", "categorical" );
            tmp.setFeature( withinCategorical[it->first][jt->first] );
            feature_layers["within"]["categorical"][i].push_back(tmp);
        }
    }
    
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        feature_layers["within"]["quantitative"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t i = it->first - 1;
        // loop over time-feature maps
        for (auto jt = withinQuantitative[it->first].begin(); jt != withinQuantitative[it->first].end(); jt++) {
            RegionalFeatureLayer tmp( it->first, jt->first, "within", "quantitative" );
            tmp.setFeature( withinQuantitative[it->first][jt->first] );
            feature_layers["within"]["quantitative"][i].push_back(tmp);
        }
    }
    
    normalizeWithinQuantitative();
    
    for (auto it = betweenCategorical.begin(); it != betweenCategorical.end(); it++) {
        feature_layers["between"]["categorical"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t i = it->first - 1;
        // loop over time-feature maps
        for (auto jt = betweenCategorical[it->first].begin(); jt != betweenCategorical[it->first].end(); jt++) {
            RegionalFeatureLayer tmp( it->first, jt->first, "between", "categorical" );
            tmp.setFeature( betweenCategorical[it->first][jt->first] );
            feature_layers["between"]["categorical"][i].push_back(tmp);
        }
    }
    
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        feature_layers["between"]["quantitative"].push_back( std::vector<RegionalFeatureLayer>() );
        size_t i = it->first - 1;
        // loop over time-feature maps
        for (auto jt = betweenQuantitative[it->first].begin(); jt != betweenQuantitative[it->first].end(); jt++) {
            RegionalFeatureLayer tmp( it->first, jt->first, "between", "quantitative" );
            tmp.setFeature( betweenQuantitative[it->first][jt->first] );
            feature_layers["between"]["quantitative"][i].push_back(tmp);
        }
    }
    
    std::cout << "initialized!\n";
}

void RegionalFeatures::normalizeWithinQuantitative(void) {
    double m = 0.0;
    size_t n_elem = 0;
    
    // get product for geometric mean
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        size_t i = it->first - 1;
        for (auto jt = withinQuantitative[it->first].begin(); jt != withinQuantitative[it->first].end(); jt++) {
            size_t j = jt->first - 1;
            std::vector<double> v = feature_layers["within"]["quantitative"][i][j].within_quantitative;
            for (size_t a = 0; a < v.size(); a++) {
                if (v[a] > 0) {
                    m += v[a];
                    n_elem += 1;
                }
            }
        }
    }
    double z = m / n_elem;
    for (auto it = withinQuantitative.begin(); it != withinQuantitative.end(); it++) {
        size_t i = it->first - 1;
        for (auto jt = withinQuantitative[it->first].begin(); jt != withinQuantitative[it->first].end(); jt++) {
            size_t j = jt->first - 1;
            std::vector<double> v = feature_layers["within"]["quantitative"][i][j].within_quantitative;
            for (size_t a = 0; a < v.size(); a++) {
                v[a] = v[a] / z;;
            }
            feature_layers["within"]["quantitative"][i][j].within_quantitative = v;
        }
    }
    return;
}

void RegionalFeatures::normalizeBetweenQuantitative(void) {
    double m = 0.0;
    size_t n_elem = 0;
    
    // get product for geometric mean
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        size_t i = it->first - 1;
        for (auto jt = betweenQuantitative[it->first].begin(); jt != betweenQuantitative[it->first].end(); jt++) {
            size_t j = jt->first - 1;
            std::vector<std::vector<double> > v = feature_layers["between"]["quantitative"][i][j].between_quantitative;
            for (size_t a = 0; a < v.size(); a++) {
                for (size_t b = 0; b < v[a].size(); b++) {
                    if (a != b && v[a][b] > 0) {
                        m += v[a][b];
                        n_elem += 1;
                    }
                }
            }
        }
    }
    double z = m / n_elem;
    for (auto it = betweenQuantitative.begin(); it != betweenQuantitative.end(); it++) {
        size_t i = it->first - 1;
        for (auto jt = betweenQuantitative[it->first].begin(); jt != betweenQuantitative[it->first].end(); jt++) {
            size_t j = jt->first - 1;
            std::vector<std::vector<double> > v = feature_layers["between"]["quantitative"][i][j].between_quantitative;
            for (size_t a = 0; a < v.size(); a++) {
                for (size_t b = 0; b < v[a].size(); b++) {
                    if (a != b) {
                        v[a][b] = v[a][b] / z;
                    }
                }
            }
            feature_layers["between"]["quantitative"][i][j].between_quantitative = v;
        }
    }
    return;
}

RegionalFeatures::RegionalFeatures(const RegionalFeatures& a)
{
    *this = a;
}

RegionalFeatures& RegionalFeatures::operator=(const RegionalFeatures& a)
{
    if (this != &a) {
        numLayers = a.numLayers;
    }

    return *this;
}

RegionalFeatures* RegionalFeatures::clone(void) const
{
    return new RegionalFeatures(*this);
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


size_t RegionalFeatures::getNumLayers(void) const {
    return numLayers;
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

//    o << "Origination:                   " << x.getFilename() << std::endl;
//    o << "Number of epochs:              " << x.getEpochs().size() << std::endl;
//    o << "Number of areas:               " << x.getAreas()[0].size() << std::endl;
//    o << std::endl;

    return o;
}


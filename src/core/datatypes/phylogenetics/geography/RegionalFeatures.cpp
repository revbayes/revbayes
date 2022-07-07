//
//  RegionalFeatures.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "RegionalFeatures.h"

#include <vector>
#include <sstream>
#include <string>

// #include "RegionalFeaturesDataReader.h"
#include "GeographicArea.h"

using namespace RevBayesCore;


RegionalFeatures::RegionalFeatures(std::vector<std::string> fn, std::vector<std::string> fr, std::vector<std::string> ft) :
    featureNames(fn),
    featureRelationships(fr),
    featureTypes(ft)
{
    for (size_t i = 0; i < featureNames.size(); i++) {
//        featureLayers.push_back( RegionalFeatureLayer( featureNames[i], featureRelationships[i], featureTypes[i] );
    }
    
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

size_t RegionalFeatures::getNumLayers(void) const {
    return numLayers;
}

//std::vector<double> RegionalFeatures::getEpochs(void) const
//{
//    return epochs;
//}

//std::vector<std::vector<GeographicArea*> > RegionalFeatures::getAreas(void) const
//{
//    return areas;
//}

//std::string RegionalFeatures::getFilename(void) const
//{
//    return filename;
//}

//size_t RegionalFeatures::getNumEpochs(void) const
//{
//    return epochs.size();
//}

//size_t RegionalFeatures::getNumAreas(void) const
//{
//    return areas[0].size();
//}

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


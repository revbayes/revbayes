//
//  RegionalFeatureLayer.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "RegionalFeatureLayer.h"
#include <iostream>
#include <sstream>

using namespace RevBayesCore;

RegionalFeatureLayer::RegionalFeatureLayer() :
    time_index(0),
    feature_index(0),
    feature_type(""),
    feature_relationship("")
{
    ;
}

RegionalFeatureLayer::RegionalFeatureLayer(size_t ti, size_t fi, std::string ft, std::string fr) :
    time_index(ti),
    feature_index(fi),
    feature_type(ft),
    feature_relationship(fr)
{
    ;
}

RegionalFeatureLayer::RegionalFeatureLayer(const RegionalFeatureLayer& g)
{
    if (this != &g) {
        feature_type = g.feature_type;
        feature_relationship = g.feature_relationship;
        feature_index = g.feature_index;
        time_index = g.time_index;
        within_categorical = g.within_categorical;
        within_quantitative = g.within_quantitative;
        between_categorical = g.between_categorical;
        between_quantitative = g.between_quantitative;
    }
}

RegionalFeatureLayer& RegionalFeatureLayer::operator=(const RegionalFeatureLayer &g)
{
    if (this != &g) {
        feature_type = g.feature_type;
        feature_relationship = g.feature_relationship;
        feature_index = g.feature_index;
        time_index = g.time_index;
        within_categorical = g.within_categorical;
        within_quantitative = g.within_quantitative;
        between_categorical = g.between_categorical;
        between_quantitative = g.between_quantitative;
    }
    
    return *this;
}

RegionalFeatureLayer* RegionalFeatureLayer::clone(void) const
{
    return new RegionalFeatureLayer( *this );
}

std::vector<int> RegionalFeatureLayer::getWithinCategoricalFeatures(void)
{
    return within_categorical;
}

std::vector<double> RegionalFeatureLayer::getWithinQuantitativeFeatures(void)
{
    return within_quantitative;
}

std::vector<std::vector<int> > RegionalFeatureLayer::getBetweenCategoricalFeatures(void)
{
    return between_categorical;
}

std::vector<std::vector<double> > RegionalFeatureLayer::getBetweenQuantitativeFeatures(void)
{
    return between_quantitative;
}

void RegionalFeatureLayer::setFeatures(std::vector<int> f)
{
    within_categorical = f;
}

void RegionalFeatureLayer::setFeatures(std::vector<double> f)
{
    within_quantitative = f;
}

void RegionalFeatureLayer::setFeatures(std::vector<std::vector<int> > f)
{
    between_categorical = f;
}

void RegionalFeatureLayer::setFeatures(std::vector<std::vector<double> > f)
{
    between_quantitative = f;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RegionalFeatureLayer& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << "RegionalFeatureLayer" << std::endl;
    o << s.str();
    
    return o;
}

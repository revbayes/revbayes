//
//  RegionalFeatureLayer.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "RegionalFeatureLayer.h"
#include <iostream>
#include <cstdint>
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

RegionalFeatureLayer::RegionalFeatureLayer(size_t ti, size_t fi, std::string fr, std::string ft) :
    time_index(ti),
    feature_index(fi),
    feature_relationship(fr),
    feature_type(ft)
{
    ;
}

RegionalFeatureLayer::RegionalFeatureLayer(const RegionalFeatureLayer& g)
{
    if (this != &g) {
        time_index = g.time_index;
        feature_index = g.feature_index;
        feature_relationship = g.feature_relationship;
        feature_type = g.feature_type;
        feature_values = g.feature_values;
        
        // delete...
//        within_categorical = g.within_categorical;
//        within_quantitative = g.within_quantitative;
//        between_categorical = g.between_categorical;
//        between_quantitative = g.between_quantitative;
        
    }
}

RegionalFeatureLayer& RegionalFeatureLayer::operator=(const RegionalFeatureLayer &g)
{
    if (this != &g) {
        time_index = g.time_index;
        feature_index = g.feature_index;
        feature_relationship = g.feature_relationship;
        feature_type = g.feature_type;
        feature_values = g.feature_values;
//        within_categorical = g.within_categorical;
//        within_quantitative = g.within_quantitative;
//        between_categorical = g.between_categorical;
//        between_quantitative = g.between_quantitative;
    }
    
    return *this;
}

RegionalFeatureLayer* RegionalFeatureLayer::clone(void) const
{
    return new RegionalFeatureLayer( *this );
}

/**
 * Equals operator.
 * We check the species name and the individuals name.
 */
bool RegionalFeatureLayer::operator==(const RevBayesCore::RegionalFeatureLayer &t) const
{
    
    if ( feature_index != t.feature_index )
    {
        return false;
    }
    if ( time_index != t.time_index )
    {
        return false;
    }
    if ( feature_relationship != t.feature_relationship ) {
        return false;
    }
    if ( feature_type != t.feature_type ) {
        return false;
    }

    return true;
}


/**
 * Not equals operator. We simply invert the result of the equals operation.
 */
bool RegionalFeatureLayer::operator!=(const RevBayesCore::RegionalFeatureLayer &t) const
{
    
    return !operator==(t);
}


/**
 * Less-than operator.
 * We check first the species name and then the indidivuals name.
 */
bool RegionalFeatureLayer::operator<(const RevBayesCore::RegionalFeatureLayer &t) const
{
    
    if ( time_index < t.time_index)
    {
        return true;
    }
    else if ( time_index > t.time_index )
    {
        return false;
    }
    
    if ( feature_index < t.feature_index )
    {
        return true;
    }
    else if ( feature_index > t.feature_index )
    {
        return false;
    }
    
    // by default return false.
    return false;
}



/**
 * Less-than or equals operator.
 */
bool RegionalFeatureLayer::operator<=(const RevBayesCore::RegionalFeatureLayer &t) const
{
    return operator<(t) || operator==(t);
}


/**
 * Greater-than operator.
 * We check first the species name and then the indidivuals name.
 */
bool RegionalFeatureLayer::operator>(const RevBayesCore::RegionalFeatureLayer &t) const
{
    return operator<=(t) == false;
}

/**
 * Greater-than or equals operator.
 */
bool RegionalFeatureLayer::operator>=(const RevBayesCore::RegionalFeatureLayer &t) const
{
    return operator>(t) || operator==(t);
}

//std::vector<std::int64_t> RegionalFeatureLayer::getWithinCategoricalFeatures(void) const
//{
//    return within_categorical;
//}
//
//std::vector<double> RegionalFeatureLayer::getWithinQuantitativeFeatures(void) const
//{
//    return within_quantitative;
//}
//
//std::vector<std::vector<std::int64_t> > RegionalFeatureLayer::getBetweenCategoricalFeatures(void) const
//{
//    return between_categorical;
//}
//
//std::vector<std::vector<double> > RegionalFeatureLayer::getBetweenQuantitativeFeatures(void) const
//{
//    return between_quantitative;
//}

std::vector<std::vector<double> > RegionalFeatureLayer::getFeatureValues(void) const
{
    return feature_values;
}

size_t RegionalFeatureLayer::getTimeIndex(void) const
{
    return time_index;
}

size_t RegionalFeatureLayer::getFeatureIndex(void) const
{
    return feature_index;
}

std::string RegionalFeatureLayer::getFeatureType(void) const
{
    return feature_type;
}

std::string RegionalFeatureLayer::getFeatureRelationship(void) const
{
    return feature_relationship;
}


//void RegionalFeatureLayer::setFeatures(std::vector<std::int64_t> f)
//{
//    within_categorical = f;
//}
//
//void RegionalFeatureLayer::setFeatures(std::vector<double> f)
//{
//    within_quantitative = f;
//}
//
//void RegionalFeatureLayer::setFeatures(std::vector<std::vector<std::int64_t> > f)
//{
//    between_categorical = f;
//}

void RegionalFeatureLayer::setFeatures(std::vector<std::vector<double> > f)
{
    feature_values = f;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const RegionalFeatureLayer& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << "RegionalFeatureLayer" << std::endl;
    s << "  time_index           = " << x.getTimeIndex() << std::endl;
    s << "  feature_index        = " << x.getFeatureIndex() << std::endl;
    s << "  feature_relationship = " << x.getFeatureRelationship() << std::endl;
    s << "  feature_type         = " << x.getFeatureType() << std::endl;
    
    std::string type = x.getFeatureType();
    std::string relationship = x.getFeatureRelationship();
//    if (relationship == "within" && type == "categorical") {
//        s << "  within_categorical   = " << x.getWithinCategoricalFeatures() << std::endl;
//    }
    
    o << s.str();
    
    return o;
}

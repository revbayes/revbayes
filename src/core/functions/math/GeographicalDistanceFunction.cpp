//
//  GeographicalDistanceFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "GeographicalDistanceFunction.h"
#include "RbConstants.h"
#include <cmath>

using namespace RevBayesCore;
/**
 * Geographical distance function of TypedDagNode of type double
 *
 * @param lat1 a TypedDagNode of type double
 * @param lon1 a TypedDagNode of type double
 * @param lat2 a TypedDagNode of type double
 * @param lon2 a TypedDagNode of type double
*/

GeographicalDistanceFunction::GeographicalDistanceFunction(const TypedDagNode<double> *lat1,
                                                           const TypedDagNode<double> *lon1,
                                                           const TypedDagNode<double> *lat2,
                                                           const TypedDagNode<double> *lon2) : ContinuousFunction( new double(0.0) ),
latitude1( lat1 ),
longitude1( lon1 ),
latitude2( lat2 ),
longitude2( lon2 )

{
    addParameter( latitude1 );
    addParameter( longitude1 );
    addParameter( latitude2 );
    addParameter( longitude2 );
    
}


GeographicalDistanceFunction* GeographicalDistanceFunction::clone( void ) const
{
    return new GeographicalDistanceFunction(*this);
}


void GeographicalDistanceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == latitude1)
    {
        latitude1 = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == longitude1)
    {
        longitude1 = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == latitude2)
    {
        latitude2 = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == longitude2)
    {
        longitude2 = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void GeographicalDistanceFunction::update( void )
{
    
    double to_rad = RbConstants::PI / 180.0;
    double lat1 = latitude1->getValue() * to_rad;  // phi1
    double lon1 = longitude1->getValue() * to_rad; // lambda1
    double lat2 = latitude2->getValue() * to_rad;  // phi2
    double lon2 = longitude2->getValue() * to_rad; // lambda2
   
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    
    double r = 6371;
//    double h1 = 1.0 - std::cos(dlat) / 2.0;
//    double h2 = 1.0 - std::cos(dlon) / 2.0;
    double x1 = std::sin(dlat / 2.0) * std::sin(dlat / 2.0);
    double x2 = std::cos(lat1) * std::cos(lat2) * std::sin(dlon / 2.0) * std::sin(dlon / 2.0);;
    double x = std::sqrt(x1 + x2);
    double y = 2 * r * std::asin(x);
    /*
   
    double x1 = std::cos(lat2) * std::sin(dlon);
    double x2 = std::cos(lat1) * std::sin(lat2) - std::sin(lat1) * std::cos(lat2) * std::cos(dlon);
    double x3 = std::sin(lat1) * std::sin(lat2) + std::cos(lat1) * std::cos(lat2) * std::cos(dlon);
    
    double y = std::atan2( std::sqrt( x1*x1 + x2*x2 ) / x3 );
    */
    
    
    *value = y; // std::atan( std::sqrt( x1*x1 + x2*x2 ) / x3 );
}



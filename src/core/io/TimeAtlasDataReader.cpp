//
//  TimeAtlasDataReader.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "TimeAtlasDataReader.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <utility>

#include "RbFileManager.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RlUserInterface.h"
#include "GeographicArea.h"

using namespace RevBayesCore;

TimeAtlasDataReader::TimeAtlasDataReader(const path& fn, std::string d) : DelimitedDataReader(fn, d)
{
    readJson();
}

TimeAtlasDataReader::TimeAtlasDataReader(const TimeAtlasDataReader& tadr) : DelimitedDataReader(tadr)
{
    
    areas = tadr.areas;
    epochs = tadr.epochs;
    
}

void TimeAtlasDataReader::readJson(void)
{
    
    std::ifstream readStream( filename.string() );
    if ( not readStream )
    {
        throw RbException()<<"Could not open file "<<filename;
    }
    
    try
    {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(readStream, pt);
        
        unsigned epochIndex = 0;
        BOOST_FOREACH(boost::property_tree::ptree::value_type &node, pt.get_child("epochs"))
        {
            std::vector<GeographicArea*> tmpAreas;
            
            std::string epochName = "";
            if (node.second.find("name") != node.second.not_found())
                epochName = node.second.get<std::string>("name");
            
            double epochStartAge = RbConstants::Double::neginf;
            if (node.second.find("start_age") != node.second.not_found())
                epochStartAge = node.second.get<double>("start_age");
            
            double epochEndAge = 0.0;
            if (node.second.find("end_age") != node.second.not_found())
                epochEndAge = node.second.get<double>("end_age");
            
            unsigned areaIndex = 0;
            boost::property_tree::ptree pt_areas = node.second;
            BOOST_FOREACH(boost::property_tree::ptree::value_type &node_areas, pt_areas.get_child("areas"))
            {
                std::string areaName = "";
                if (node_areas.second.find("name") != node_areas.second.not_found())
                    areaName = node_areas.second.get<std::string>("name");
                
                double areaLatitude = 0.0;
                if (node_areas.second.find("latitude") != node_areas.second.not_found())
                    areaLatitude = node_areas.second.get<double>("latitude");
                
                double areaLongitude = 0.0;
                if (node_areas.second.find("longitude") != node_areas.second.not_found())
                    areaLongitude =  node_areas.second.get<double>("longitude");
                
                double areaAltitude = 0.0;
                if (node_areas.second.find("altitude") != node_areas.second.not_found())
                    areaAltitude = node_areas.second.get<double>("altitude");
                
                double areaSize = 0.0;
                if (node_areas.second.find("size") != node_areas.second.not_found())
                    areaSize = node_areas.second.get<double>("size");
                
                std::vector<double> dispersalValues;
                if (node_areas.second.find("dispersalValues") != node_areas.second.not_found())
                {
                    boost::property_tree::ptree pt_rates = node_areas.second;
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node_rates, pt_rates.get_child("dispersalValues"))
                    {
                        double v = node_rates.second.get_value<double>();
                        dispersalValues.push_back( v );
                    }
                }
         
                std::vector<double> extinctionValues;
                if (node_areas.second.find("extinctionValues") != node_areas.second.not_found())
                {
                    boost::property_tree::ptree pt_rates = node_areas.second;
                    BOOST_FOREACH(boost::property_tree::ptree::value_type &node_rates, pt_rates.get_child("extinctionValues"))
                    {
                        double v = node_rates.second.get_value<double>();
                        extinctionValues.push_back( v );
                    }
                }
                
                
                unsigned areaState = 0;
                if (node_areas.second.find("state") != node_areas.second.not_found())
                    areaState = node_areas.second.get<int>("state");
              
                
                
                GeographicArea* g  = new GeographicArea(areaIndex, areaLatitude, areaLongitude, areaName, areaState, areaAltitude, areaSize, dispersalValues, extinctionValues);
                
                //std::cout << areaIndex << " " << areaName << " " << areaLatitude << " " << areaLongitude << " " << areaState << "\n";
                
                areaIndex++;
                tmpAreas.push_back(g);
            }
            
            epochIndex++;
            areas.push_back(tmpAreas);
            epochs.push_back(epochEndAge);
        }
        
        
        sortEpochs();
    }
    catch (std::exception const& e)
    {
        throw RbException( e.what() );
    }
}

void TimeAtlasDataReader::sortEpochs(void)
{
    std::vector<AgeIndexPair> v;
    for (size_t i = 0; i < epochs.size(); i++)
        v.push_back(AgeIndexPair(epochs[i],i));
    
    std::sort(v.begin(), v.end(), AgeCompare());
    
    std::vector<double> tmpEpochs = epochs;
    std::vector<std::vector<GeographicArea*> > tmpAreas = areas;
    for (size_t i = 0; i < v.size(); i++)
    {
        areas[i] = tmpAreas[v[i].index];
        epochs[i] = tmpEpochs[v[i].index];
        //std::cout << epochs[i] << "\n";
    }

}


void TimeAtlasDataReader::printJson(boost::property_tree::ptree const& pt)
{
    using boost::property_tree::ptree;
    ptree::const_iterator end = pt.end();
    for (ptree::const_iterator it = pt.begin(); it != end; ++it)
    {
        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        printJson(it->second);
    }
}

std::vector<double> TimeAtlasDataReader::getEpochs(void)
{
    return epochs;
}

std::vector<std::vector<GeographicArea*> > TimeAtlasDataReader::getAreas(void)
{
    return areas;
}

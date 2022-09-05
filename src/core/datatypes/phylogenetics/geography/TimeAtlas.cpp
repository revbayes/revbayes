//
//  TimeAtlas.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "TimeAtlas.h"

#include <sstream>
#include <string>

#include "TimeAtlasDataReader.h"
#include "GeographicArea.h"

using namespace RevBayesCore;

TimeAtlas::TimeAtlas(TimeAtlasDataReader* tadr) : areas(tadr->getAreas()), epochs(tadr->getEpochs()), filename(tadr->getFilename())
{
    numEpochs =(unsigned int)areas.size();
    numAreas = 0;
    if (numEpochs > 0)
        numAreas = (unsigned int)areas[0].size();
    
}

TimeAtlas::TimeAtlas(const TimeAtlas& a)
{
    *this = a;
}

TimeAtlas& TimeAtlas::operator=(const TimeAtlas& a)
{
    if (this != &a)
    {
        epochs = a.epochs;
        for (size_t i = 0; i < a.areas.size(); i++)
        {
            std::vector<GeographicArea*> tmpAreas;
            for (size_t j = 0; j < a.areas[i].size(); j++)
            {
                tmpAreas.push_back(new GeographicArea(*a.areas[i][j]));
            }
            areas.push_back(tmpAreas);
        }
        
        numAreas = a.numAreas;
        numEpochs = a.numEpochs;
        filename = a.filename;
    }
    
    return *this;
}

TimeAtlas* TimeAtlas::clone(void) const
{
    return new TimeAtlas(*this);
}

std::vector<double> TimeAtlas::getEpochs(void) const
{
    return epochs;
}

std::vector<std::vector<GeographicArea*> > TimeAtlas::getAreas(void) const
{
    return areas;
}

const path& TimeAtlas::getFilename(void) const
{
    return filename;
}

size_t TimeAtlas::getNumEpochs(void) const
{
    return epochs.size();
}

size_t TimeAtlas::getNumAreas(void) const
{
    return areas[0].size();
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const TimeAtlas& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << "TimeAtlas with " << x.getAreas().size() << " epochs and " << x.getAreas().size() << " areas" << std::endl;
    o << s.str();
    
    for ( size_t i = 0; i < s.str().length() - 1; ++i )
        o << "=";
    o << std::endl;
    
    o << "Origination:                   " << x.getFilename() << std::endl;
    o << "Number of epochs:              " << x.getEpochs().size() << std::endl;
    o << "Number of areas:               " << x.getAreas()[0].size() << std::endl;
    o << std::endl;
    
    return o;
}

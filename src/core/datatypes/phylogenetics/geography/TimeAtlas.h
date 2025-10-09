//
//  TimeAtlas.h
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__TimeAtlas__
#define __rb_mlandis__TimeAtlas__

#include <cstddef>
#include <iosfwd>
#include <vector>
#include "RbFileManager.h"

#include "Cloneable.h"

namespace RevBayesCore {
class GeographicArea;
class TimeAtlasDataReader;

    class TimeAtlas : public Cloneable
    {
        
    public:
        TimeAtlas(TimeAtlasDataReader* tadr);
        TimeAtlas(const TimeAtlas& a);
        TimeAtlas&                                      operator=(const TimeAtlas& a);
        virtual TimeAtlas*                              clone(void) const;
        
        std::vector<double>                             getEpochs(void) const;
        std::vector<std::vector<GeographicArea*> >      getAreas(void) const;
        const path&                                     getFilename(void) const;
        std::string                                     getDataType(void) const;  
        size_t                                          getNumEpochs(void) const;
        size_t                                          getNumAreas(void) const;
        
    protected:
        std::vector<std::vector<GeographicArea*> >      areas;
        std::vector<double>                             epochs;
        
    private:
        unsigned                                        numAreas;
        unsigned                                        numEpochs;
        path                                            filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const TimeAtlas& x);
}


#endif /* defined(__rb_mlandis__TimeAtlas__) */

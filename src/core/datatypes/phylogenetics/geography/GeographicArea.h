//
//  GeographicArea.h
//  rb_mlandis
//
//  Created by Michael Landis on 11/25/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__GeographicArea__
#define __rb_mlandis__GeographicArea__

#include <cstddef>
#include <set>
#include <vector>
#include <iosfwd>

#include "Cloneable.h"

namespace RevBayesCore {

    class GeographicArea : public Cloneable
    {
    public:
        GeographicArea(size_t idx, double lat, double lon, std::string n="", int state=0, double alt=1.0, double size=1.0, std::vector<double> dv=std::vector<double>(), std::vector<double> ev=std::vector<double>(), std::set<GeographicArea*> aa=std::set<GeographicArea*>(), std::set<GeographicArea*> ra=std::set<GeographicArea*>());
//        areaIndex, areaLatitude, areaLongitude, areaName, areaState, areaAltitude, areaSize, areaRates
        GeographicArea(const GeographicArea& g);
        GeographicArea&                                                         operator=(const GeographicArea& g);
        virtual GeographicArea*                                                 clone(void) const;
        std::string                                                             getName(void);
        size_t                                                                  getIndex(void);
        double                                                                  getLatitude(void);
        double                                                                  getLongitude(void);
        double                                                                  getAltitude(void);
        const std::vector<double>&                                              getLatlon(void);
        double                                                                  getSize(void);
        int                                                                     getState(void);
        const std::vector<double>&                                              getDispersalValues(void);
        const std::vector<double>&                                              getExtinctionValues(void);
        const std::set<GeographicArea*>&                                        getAdjacentAreas(void);
        const std::set<GeographicArea*>&                                        getReachableAreas(void);
        void                                                                    addAdjacentArea(GeographicArea* a);
        void                                                                    removeAdjacentArea(GeographicArea* a);
        void                                                                    swapAdjacentAreas(GeographicArea* oldArea, GeographicArea* newArea);
        void                                                                    swapReachableAreas(GeographicArea* oldArea, GeographicArea* newArea);
        std::string                                                             getStr(void);
        
    protected:
        size_t                      index;
        std::string                 name;
        int                         state;
        double                      latitude;
        double                      longitude;
        std::vector<double>         latlon;
        double                      altitude;
        double                      size;
        std::vector<double>         dispersalValues;
        std::vector<double>         extinctionValues;
        std::set<GeographicArea*>   adjacentAreas;
        std::set<GeographicArea*>   reachableAreas;
    };

    std::ostream&                                       operator<<(std::ostream& o, const GeographicArea& x);
}

#endif /* defined(__rb_mlandis__GeographicArea__) */

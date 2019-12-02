//
//  RelativeNodeAgeWeightedConstraints.h
//
//  Created by Bastien Boussau on 4/8/15.
//  Copyright (c) 2015 Bastien Boussau. All rights reserved.
//

#ifndef __RelativeNodeAgeWeightedConstraints__
#define __RelativeNodeAgeWeightedConstraints__

#include "Cloneable.h"
#include "RelativeNodeAgeWeightedConstraintsReader.h"
#include <string>

namespace RevBayesCore {
    class RelativeNodeAgeWeightedConstraints : public Cloneable
    {
        
    public:
        RelativeNodeAgeWeightedConstraints();
        RelativeNodeAgeWeightedConstraints(RelativeNodeAgeWeightedConstraintsReader* tadr);
        RelativeNodeAgeWeightedConstraints(const RelativeNodeAgeWeightedConstraints& a);
        RelativeNodeAgeWeightedConstraints&                     operator=(const RelativeNodeAgeWeightedConstraints& a);
        virtual RelativeNodeAgeWeightedConstraints*             clone(void) const;
        
        const std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double >& getConstraint(size_t i) const;
        const std::vector <std::pair<std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > >& getConstraints( ) const;

        size_t                                    getNumberOfConstraints(void) const;
        std::string                                     getFilename(void) const;
       // std::string                                     getDatatype(void) const;
    protected:
        std::vector < std::pair <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > >								        olderYoungerConstraints;
        
    private:
        std::string                                     filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RelativeNodeAgeWeightedConstraints& x);
}


#endif /* defined(__RelativeNodeAgeWeightedConstraints__) */

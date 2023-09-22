//
//  RelativeNodeAgeWeightedConstraints.h
//
//  Created by Bastien Boussau on 4/8/15.
//  Copyright (c) 2015 Bastien Boussau. All rights reserved.
//

#ifndef __RelativeNodeAgeWeightedConstraints__
#define __RelativeNodeAgeWeightedConstraints__

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>

#include "RbFileManager.h"
#include "Cloneable.h"

namespace RevBayesCore {
class RelativeNodeAgeWeightedConstraintsReader;

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
        const path&                                       getFilename(void) const;
       // std::string                                     getDatatype(void) const;
    protected:
        std::vector < std::pair <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > >								        olderYoungerConstraints;
        
    private:
        path                                            filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RelativeNodeAgeWeightedConstraints& x);
}


#endif /* defined(__RelativeNodeAgeWeightedConstraints__) */

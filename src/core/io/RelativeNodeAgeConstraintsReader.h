#ifndef RelativeNodeAgeConstraintsReader_H
#define RelativeNodeAgeConstraintsReader_H

#include <cstddef>
#include <vector>
#include <iosfwd>
#include <utility>

#include "DelimitedDataReader.h"

namespace RevBayesCore {
    
    
    /**
     * Reader for relative constraints on node ages.
     * Each line corresponds to two nodes, each defined using a pair of tip names separated by spaces (if a tip needs to be defined, it can be defined by repeating the tip name).
     * It is understood that the first node is older (closer to the root) than the second node.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Bastien Boussau)
     * @since 2015-03-03, version 1.0
     *
     */
    class RelativeNodeAgeConstraintsReader : public DelimitedDataReader {
        
    public:
        
        RelativeNodeAgeConstraintsReader(const std::string &fn, std::string d="", size_t ns=0);
        
        const std::vector< std::pair < std::pair<std::string, std::string >, std::pair<std::string, std::string > > >& getConstraints(void);
        const size_t                                                                                     getNumberOfConstraints(void);
        
        
    protected:
        
        std::vector< std::pair < std::pair < std::string, std::string >, std::pair < std::string, std::string > > > OlderYoungerConstraints;
    };
    
}

#endif

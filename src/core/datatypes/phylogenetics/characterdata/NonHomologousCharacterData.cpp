//
//  NonHomologousCharacterData.cpp
//  RevBayes
//
//  Created by Sebastian Hoehna on 3/31/15.
//  Copyright (c) 2015 Sebastian Hoehna. All rights reserved.
//

#include "NonHomologousCharacterData.h"

#include <sstream>
#include <string>

std::ostream& RevBayesCore::operator<<(std::ostream& o, const NonHomologousCharacterData& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << x.getDataType() << " character matrix with " << x.getNumberOfTaxa() << " taxa" << std::endl;
    o << s.str();
    
    for ( size_t i = 0; i < s.str().length() - 1; ++i )
    {
        o << "=";
    }
    o << std::endl;
    
    o << "Origination:                   " << x.getFilename().filename() << std::endl;
    o << "Number of taxa:                " << x.getNumberOfTaxa() << std::endl;
    o << "Number of included taxa:       " << x.getNumberOfIncludedTaxa() << std::endl;
    o << "Datatype:                      " << x.getDataType() << std::endl;
    o << std::endl;
    
    return o;
}

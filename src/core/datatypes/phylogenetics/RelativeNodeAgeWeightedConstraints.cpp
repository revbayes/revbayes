//
//  RelativeNodeAgeWeightedConstraints.cpp
//
//  Created by Bastien Boussau on 4/8/15.
//  Copyright (c) 2015 Bastien Boussau. All rights reserved.
//

#include "RelativeNodeAgeWeightedConstraints.h"

#include <sstream>

#include "RelativeNodeAgeWeightedConstraintsReader.h"

using namespace RevBayesCore;

RelativeNodeAgeWeightedConstraints::RelativeNodeAgeWeightedConstraints(  )
{
    
}


RelativeNodeAgeWeightedConstraints::RelativeNodeAgeWeightedConstraints(RelativeNodeAgeWeightedConstraintsReader* tadr) : filename(tadr->getFilename())
{
    olderYoungerConstraints = tadr->getConstraints();
}

RelativeNodeAgeWeightedConstraints::RelativeNodeAgeWeightedConstraints(const RelativeNodeAgeWeightedConstraints& a)
{
    *this = a;
}


RelativeNodeAgeWeightedConstraints& RelativeNodeAgeWeightedConstraints::operator=(const RelativeNodeAgeWeightedConstraints& a)
{
    if (this != &a)
    {
        olderYoungerConstraints = a.olderYoungerConstraints;
    }
    
    return *this;
}

RelativeNodeAgeWeightedConstraints* RelativeNodeAgeWeightedConstraints::clone(void) const
{
    return new RelativeNodeAgeWeightedConstraints(*this);
}

const std::vector < std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > >& RelativeNodeAgeWeightedConstraints::getConstraints(void) const
{
    return olderYoungerConstraints;
}

const std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double >& RelativeNodeAgeWeightedConstraints::getConstraint(size_t i) const
{
    return olderYoungerConstraints[i];
}

size_t RelativeNodeAgeWeightedConstraints::getNumberOfConstraints(void) const
{
    return olderYoungerConstraints.size();
}


const path& RelativeNodeAgeWeightedConstraints::getFilename(void) const
{
    return filename;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const RelativeNodeAgeWeightedConstraints& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    
    s << "RelativeNodeAgeWeightedConstraints with " << x.getNumberOfConstraints() << " constraints. " << std::endl;
    
    o << s.str();
    
    std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > constraint;
    for ( size_t i = 0; i < x.getNumberOfConstraints(); ++i ) {
        constraint = x.getConstraint(i) ;
            o << "\t" << constraint.first.first.first << "\t" << constraint.first.first.second << "\t" << constraint.first.second.first << "\t" << constraint.first.second.second << "\t" << constraint.second << "\n" ;
        }

    o << std::endl;
    
    return o;
}

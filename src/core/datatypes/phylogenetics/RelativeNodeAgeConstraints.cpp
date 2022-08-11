#include "RelativeNodeAgeConstraints.h"

#include <sstream>

#include "RelativeNodeAgeConstraintsReader.h"

using namespace RevBayesCore;

RelativeNodeAgeConstraints::RelativeNodeAgeConstraints(  )
{
    
}


RelativeNodeAgeConstraints::RelativeNodeAgeConstraints(RelativeNodeAgeConstraintsReader* tadr) : filename(tadr->getFilename())
{
    olderYoungerConstraints = tadr->getConstraints();
}

RelativeNodeAgeConstraints::RelativeNodeAgeConstraints(const RelativeNodeAgeConstraints& a)
{
    *this = a;
}


RelativeNodeAgeConstraints& RelativeNodeAgeConstraints::operator=(const RelativeNodeAgeConstraints& a)
{
    if (this != &a)
    {
        olderYoungerConstraints = a.olderYoungerConstraints;
    }
    
    return *this;
}

RelativeNodeAgeConstraints* RelativeNodeAgeConstraints::clone(void) const
{
    return new RelativeNodeAgeConstraints(*this);
}

const std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > >& RelativeNodeAgeConstraints::getConstraints(void) const
{
    return olderYoungerConstraints;
}

const std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >& RelativeNodeAgeConstraints::getConstraint(size_t i) const
{
    return olderYoungerConstraints[i];
}

size_t RelativeNodeAgeConstraints::getNumberOfConstraints(void) const
{
    return olderYoungerConstraints.size();
}


const path& RelativeNodeAgeConstraints::getFilename(void) const
{
    return filename;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const RelativeNodeAgeConstraints& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    
    s << "RelativeNodeAgeConstraints with " << x.getNumberOfConstraints() << " constraints. " << std::endl;
    
    o << s.str();
    
    std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > constraint;
    for ( size_t i = 0; i < x.getNumberOfConstraints(); ++i ) {
        constraint = x.getConstraint(i) ;
            o << "\t" << constraint.first.first << "\t" << constraint.first.second << "\t" << constraint.second.first << "\t" << constraint.second.second << "\n" ;
        }

    o << std::endl;
    
    return o;
}

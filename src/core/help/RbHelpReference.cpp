#include "RbHelpReference.h"

using namespace RevBayesCore;

const std::string& RbHelpReference::getCitation( void ) const
{
    // return a const reference to the internal value
    return citation;
}


const std::string& RbHelpReference::getDoi( void ) const
{
    // return a const reference to the internal value
    return doi;
}


const std::string& RbHelpReference::getUrl( void ) const
{
    // return a const reference to the internal value
    return url;
}


void RbHelpReference::setCitation(const std::string &c)
{
    // replace the internal value with this new value
    citation = c;
}


void RbHelpReference::setDoi(const std::string &d)
{
    // replace the internal value with this new value
    doi = d;
}


void RbHelpReference::setUrl(const std::string &u)
{
    // replace the internal value with this new value
    url = u;
}

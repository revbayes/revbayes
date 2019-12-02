#include "RbHelpEntry.h"

using namespace RevBayesCore;




const std::vector<std::string>& RbHelpEntry::getAliases( void ) const
{
    // return a const reference to the internal value
    return aliases;
}


const std::vector<std::string>& RbHelpEntry::getAuthor( void ) const
{
    // return a const reference to the internal value
    return author;
}

const std::vector<std::string>& RbHelpEntry::getTypeSpec( void ) const
{
    // return a const reference to the internal value
    return typeSpec;
}

const std::string& RbHelpEntry::getDescription( void ) const
{
    // return a const reference to the internal value
    return description;
}


const std::string& RbHelpEntry::getDetails( void ) const
{
    // return a const reference to the internal value
    return details;
}


const std::string& RbHelpEntry::getExample( void ) const
{
    // return a const reference to the internal value
    return example;
}


const std::string& RbHelpEntry::getName( void ) const
{
    // return a const reference to the internal value
    return name;
}


const std::vector<RbHelpReference>& RbHelpEntry::getReferences( void ) const
{
    // return a const reference to the internal value
    return references;
}


const std::vector<std::string>& RbHelpEntry::getSeeAlso( void ) const
{
    // return a const reference to the internal value
    return seeAlso;
}


const std::string& RbHelpEntry::getTitle( void ) const
{
    // return a const reference to the internal value
    return title;
}


void RbHelpEntry::setAliases(const std::vector<std::string> &a)
{
    // overwrite the internal value with the new value
    aliases = a;
}


void RbHelpEntry::setAuthor(const std::vector<std::string> &a)
{
    // overwrite the internal value with the new value
    author = a;
}


void RbHelpEntry::setTypeSpec(const std::vector<std::string> &a)
{
    // overwrite the internal value with the new value
    typeSpec = a;
}

void RbHelpEntry::setDescription(const std::string &d)
{
    // overwrite the internal value with the new value
    description = d;
}


void RbHelpEntry::setDetails(const std::string &d)
{
    // overwrite the internal value with the new value
    details = d;
}


void RbHelpEntry::setExample(const std::string &e)
{
    // overwrite the internal value with the new value
    example = e;
}


void RbHelpEntry::setName(const std::string &n)
{
    // overwrite the internal value with the new value
    name = n;
}


void RbHelpEntry::setReferences(const std::vector<RbHelpReference> &r)
{
    // overwrite the internal value with the new value
    references = r;
}


void RbHelpEntry::setSeeAlso(const std::vector<std::string> &s)
{
    // overwrite the internal value with the new value
    seeAlso = s;
}


void RbHelpEntry::setTitle(const std::string &t)
{
    // overwrite the internal value with the new value
    title = t;
}



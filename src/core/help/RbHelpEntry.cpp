#include "RbHelpEntry.h"
#include "RbHelpDatabase.h"

using namespace RevBayesCore;




const std::vector<std::string>& RbHelpEntry::getAliases( void ) const
{
    // return a const reference to the internal value
    return aliases;
}


const std::vector<std::string>& RbHelpEntry::getAuthor( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpStringVector(name,"author");
}

const std::vector<std::string>& RbHelpEntry::getTypeSpec( void ) const
{
    // return a const reference to the internal value
    return typeSpec;
}

const std::string& RbHelpEntry::getDescription( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpString(name,"description");
}


const std::string& RbHelpEntry::getDetails( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpString(name,"details");
}


const std::string& RbHelpEntry::getExample( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpString(name,"example");
}


const std::string& RbHelpEntry::getName( void ) const
{
    // return a const reference to the internal value
    return name;
}


const std::vector<RbHelpReference>& RbHelpEntry::getReferences( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpReferences(name);
}


const std::vector<std::string>& RbHelpEntry::getSeeAlso( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpStringVector(name,"see_also");
}


const std::string& RbHelpEntry::getTitle( void ) const
{
    // return a const reference to the internal value
    return RbHelpDatabase::getHelpDatabase().getHelpString(name,"title");
}


void RbHelpEntry::setAliases(const std::vector<std::string> &a)
{
    // overwrite the internal value with the new value
    aliases = a;
}


void RbHelpEntry::setTypeSpec(const std::vector<std::string> &a)
{
    // overwrite the internal value with the new value
    typeSpec = a;
}


void RbHelpEntry::setName(const std::string &n)
{
    // overwrite the internal value with the new value
    name = n;
}



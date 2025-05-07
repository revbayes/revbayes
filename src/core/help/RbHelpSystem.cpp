#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "RbException.h"
#include "RbFileManager.h"
#include "RbHelpSystem.h"
#include "RbHelpDistribution.h"
#include "RbHelpFunction.h"
#include "RbHelpType.h"

namespace RevBayesCore { class RbHelpEntry; }

using namespace RevBayesCore;


RbHelpSystem::RbHelpSystem()
{
    
}


RbHelpSystem::RbHelpSystem( const RbHelpSystem &hs) :
    helpForFunctions( hs.helpForFunctions ),
    helpForMethods( hs.helpForMethods ),
    helpForTypes(  ),
    helpFunctionNames( hs.helpFunctionNames ),
    helpTypeNames( hs.helpTypeNames )
{
    
    
    
    for ( std::map<std::string, RbHelpType*>::const_iterator it = hs.helpForTypes.begin(); it != hs.helpForTypes.end(); ++it)
    {
        helpForTypes.insert( std::pair<std::string, RbHelpType*>( it->first, it->second->clone() ) );
    }
    
}


RbHelpSystem::~RbHelpSystem( void )
{
    
    // free the pointers to the help types
    for ( std::map<std::string, RbHelpType*>::const_iterator it = helpForTypes.begin(); it != helpForTypes.end(); ++it)
    {
        RbHelpType *ht = it->second;
        delete ht;
    }
    helpForTypes.clear();
    
}




RbHelpSystem& RbHelpSystem::operator=( const RbHelpSystem &hs )
{
    
    if ( this != &hs )
    {
        
        // free the pointers to the help types
        for ( std::map<std::string, RbHelpType*>::const_iterator it = helpForTypes.begin(); it != helpForTypes.end(); ++it)
        {
            delete it->second;
        }
        helpForTypes.clear();
        
        helpFunctionNames   = hs.helpFunctionNames;
        helpForFunctions    = hs.helpForFunctions;
        helpTypeNames       = hs.helpTypeNames;
        helpForMethods      = hs.helpForMethods;
        
        for ( std::map<std::string, RbHelpType*>::const_iterator it = hs.helpForTypes.begin(); it != hs.helpForTypes.end(); ++it)
        {
            helpForTypes.insert( std::pair<std::string, RbHelpType*>( it->first, it->second->clone() ) );
        }
    }
    
    return *this;
}


void RbHelpSystem::addHelpDistribution( RbHelpDistribution *h)
{
    
//    helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( h.getName() , h) );
//    helpFunctionNames.insert( h.getName() );
//    
//    // also add all aliases
//    const std::vector<std::string>& aliases = h.getAliases();
//    for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
//    {
//        helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( *alias , h) );
//    }

    
    if ( h != NULL && helpForTypes.find( h->getName() ) == helpForTypes.end() )
    {
        helpForTypes.insert( std::pair<std::string,RbHelpType*>( h->getName() , h ) );
        helpTypeNames.insert( h->getName() );
    
        // create a map for all methods for this type
        std::map<std::string, RbHelpFunction> methodsHelp;
        const std::vector<RbHelpFunction>& method = h->getMethods();
        for (std::vector<RbHelpFunction>::const_iterator m = method.begin(); m != method.end(); ++m)
        {
            methodsHelp.insert( std::pair<std::string,RbHelpFunction>( m->getName() , *m) );
        }
    
        // add the methods to our global map
        helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(h->getName(),methodsHelp) );
    
    
        // also add all aliases
        const std::vector<std::string>& aliases = h->getAliases();
        for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
        {
            if ( helpForTypes.find(*alias ) == helpForTypes.end() )
            {
                helpForTypes.insert( std::pair<std::string,RbHelpType*>( *alias , h->clone() ) );
                helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(*alias,methodsHelp) );
            }
            
        }

    }
    else if ( h != NULL )
    {
        delete h;
    }

}


void RbHelpSystem::addHelpFunction( RbHelpFunction *h )
{
    
    if ( h != NULL )
    {
        helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( h->getName() , *h) );
        helpFunctionNames.insert( h->getName() );
    
        // also add all aliases
        const std::vector<std::string>& aliases = h->getAliases();
        for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
        {
            helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( *alias , *h) );
        }
    }
    
    delete h;
}




void RbHelpSystem::addHelpType( RbHelpType *h )
{
    
    
    if ( h != NULL && helpForTypes.find( h->getName() ) == helpForTypes.end() )
    {
        
        helpForTypes.insert( std::pair<std::string,RbHelpType*>( h->getName() , h ) );
        helpTypeNames.insert( h->getName() );
    
        // create a map for all methods for this type
        std::map<std::string, RbHelpFunction> methodsHelp;
        const std::vector<RbHelpFunction>& method = h->getMethods();
        for (std::vector<RbHelpFunction>::const_iterator m = method.begin(); m != method.end(); ++m)
        {
            methodsHelp.insert( std::pair<std::string,RbHelpFunction>( m->getName() , *m) );
        }
    
        // add the methods to our global map
        helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(h->getName(),methodsHelp) );
    
    
        // also add all aliases
        const std::vector<std::string>& aliases = h->getAliases();
        for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
        {
            helpForTypes.insert( std::pair<std::string,RbHelpType*>( *alias , h->clone() ) );
            helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(*alias,methodsHelp) );
        }
    
    }
    else if ( h != NULL )
    {
        delete h;
    }
    
}


const std::set<std::string>& RbHelpSystem::getFunctionEntries( void ) const
{
    // return a constant reference to the internal value
    return helpFunctionNames;
}


const std::set<std::string>& RbHelpSystem::getTypeEntries( void ) const
{
    // return a constant reference to the internal value
    return helpTypeNames;
}


/** Retrieve the help entry */
const RbHelpEntry& RbHelpSystem::getHelp(const std::string &qs)
{
    
    std::map<std::string, RbHelpFunction>::iterator itFunction = helpForFunctions.find( qs );
    std::map<std::string, RbHelpType*>::iterator itType = helpForTypes.find( qs );
    if ( itFunction != helpForFunctions.end() )
    {
        return itFunction->second;
    }
    else if ( itType != helpForTypes.end() )
    {
        return *(itType->second);
    }
    else
    {
        throw RbException("Could not find help for '" + qs + "'.");
    }
    
}


/** Retrieve the help entry */
const RbHelpEntry& RbHelpSystem::getHelp(const std::string &baseQuery, const std::string &qs)
{
    
    // find the corresponding base type
    std::map<std::string, std::map<std::string, RbHelpFunction> >::iterator itMethods = helpForMethods.find( baseQuery );
    if ( itMethods != helpForMethods.end() )
    {
        const std::map<std::string, RbHelpFunction> &methods = itMethods->second;
        
        // find the corresponding method
        std::map<std::string, RbHelpFunction>::const_iterator itFunction = methods.find( qs );
        
        // check if we found it
        if ( itFunction != methods.end() )
        {
            return itFunction->second;
        }
        
    }

    throw RbException("Could not find help for '" + baseQuery + "." + qs + "'.");
    
}




void RbHelpSystem::addHelpEntry( void )
{

//    if ( parser.testHelpEntry( *it ) == RbHelpParser::FUNCTION )
//    {
//            
//        RbHelpFunction h = parser.parseHelpFunction( *it );
//        helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( h.getName() , h) );
//        helpFunctionNames.insert( h.getName() );
//            
//            
//        // also add all aliases
//        const std::vector<std::string>& aliases = h.getAliases();
//        for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
//        {
//            helpForFunctions.insert( std::pair<std::string,RbHelpFunction>( *alias , h) );
//        }
//            
//    }
//        
//        
//        
//    if ( parser.testHelpEntry( *it ) == RbHelpParser::TYPE || parser.testHelpEntry( *it ) == RbHelpParser::DISTRIBUTION || parser.testHelpEntry( *it ) == RbHelpParser::MONITOR || parser.testHelpEntry( *it ) == RbHelpParser::MOVE )
//    {
//            
//        RbHelpType* h = NULL;
//        if ( parser.testHelpEntry( *it ) == RbHelpParser::TYPE )
//        {
//            h = parser.parseHelpType( *it );
//        }
//        else if ( parser.testHelpEntry( *it ) == RbHelpParser::DISTRIBUTION )
//        {
//            h = parser.parseHelpDistribution( *it );
//        }
//        else if ( parser.testHelpEntry( *it ) == RbHelpParser::MONITOR )
//        {
//            h = parser.parseHelpMonitor( *it );
//        }
//        else if ( parser.testHelpEntry( *it ) == RbHelpParser::MOVE )
//        {
//            h = parser.parseHelpMove( *it );
//        }
//        
//        helpForTypes.insert( std::pair<std::string,RbHelpType*>( h->getName() , h) );
//        helpTypeNames.insert( h->getName() );
//            
//            
//        // create a map for all methods for this type
//        std::map<std::string, RbHelpFunction> methodsHelp;
//        const std::vector<RbHelpFunction>& method = h->getMethods();
//        for (std::vector<RbHelpFunction>::const_iterator m = method.begin(); m != method.end(); ++m)
//        {
//            methodsHelp.insert( std::pair<std::string,RbHelpFunction>( m->getName() , *m) );
//        }
//            
//        // add the methods to our global map
//        helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(h->getName(),methodsHelp) );
//            
//            
//        // also add all aliases
//        const std::vector<std::string>& aliases = h->getAliases();
//        for (std::vector<std::string>::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias)
//        {
//            helpForTypes.insert( std::pair<std::string,RbHelpType*>( *alias , h->clone() ) );
//            helpForMethods.insert( std::pair<std::string, std::map<std::string,RbHelpFunction> >(*alias,methodsHelp) );
//        }
//    }
    
    
}


bool RbHelpSystem::isHelpAvailableForQuery(const std::string &query)
{
    // test if we have a help entry for this query string
    return helpForFunctions.find( query ) != helpForFunctions.end() || helpForTypes.find( query ) != helpForTypes.end();
}


bool RbHelpSystem::isHelpAvailableForQuery(const std::string &baseQuery, const std::string &query)
{
    // test if we have a help entry for this query string
    return helpForTypes.find( baseQuery ) != helpForTypes.end();
}

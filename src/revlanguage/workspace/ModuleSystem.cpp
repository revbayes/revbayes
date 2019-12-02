#include "RbException.h"
#include "RbFileManager.h"
#include "ModuleSystem.h"
#include "RbUtil.h"
#include "StringUtilities.h"

#include <iostream>
#include <sstream>

using namespace RevLanguage;


ModuleSystem::ModuleSystem()
{
    
}


ModuleSystem::ModuleSystem( const ModuleSystem &ms) :
    modules( ms.modules )
{
    
}


ModuleSystem::~ModuleSystem( void )
{
    
}




ModuleSystem& ModuleSystem::operator=( const ModuleSystem &ms )
{
    
    if ( this != &ms )
    {
        
        modules = ms.modules;
    }
    
    return *this;
}


/** Retrieve the module */
const Module& ModuleSystem::getModule(const std::string &qs) const
{
    std::map<std::string, Module>::const_iterator it = modules.find( qs );
    if ( it != modules.end() )
    {
        return it->second;
    }
    else
    {
        throw RbException("Could not find module with name '" + qs + "'.");
    }
}



/** Initialize the modules from an text filew */
void ModuleSystem::loadModules(const std::string &dir)
{
    
    // find the path to the directory containing the help files
    RevBayesCore::RbFileManager fMngr = RevBayesCore::RbFileManager();
    
    
    fMngr.setFilePath( dir );
    if (fMngr.testDirectory() == false)
    {
//        throw RbException("Warning: Cannot find directory containing modules. No modules are available. Path = " + dir);
    }
    else
    {
        // get the files contained in the directory
    
        // gather all text files in dir, filtered by '.ext'
        std::string ext = "Rev";
        std::vector<std::string> files;
        std::vector<std::string> fileNames;
        fMngr.setStringWithNamesOfFilesInDirectory( files );
        for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it)
        {
            RevBayesCore::RbFileManager tmpFM = RevBayesCore::RbFileManager( *it );
            if ( tmpFM.getFileExtension() == ext)
            {
                fileNames.push_back( *it );
            }
        }
    
        for (std::vector<std::string>::iterator it = fileNames.begin(); it != fileNames.end(); ++it)
        {
            Module m = Module(*it);
            RevBayesCore::RbFileManager tmp = RevBayesCore::RbFileManager( *it );
            std::string name = tmp.getFileNameWithoutExtension();
            modules.insert( std::pair<std::string, Module>(name,m) );
        }
    
    }
    
}


bool ModuleSystem::isModuleAvailable(const std::string &query)
{
    // test if we have a help entry for this query string
    return modules.find( query ) != modules.end();
}


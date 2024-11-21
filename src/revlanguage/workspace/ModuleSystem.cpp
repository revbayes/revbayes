#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "RbException.h"
#include "RbFileManager.h"
#include "ModuleSystem.h"
#include "Module.h"

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
        throw RbException() << "Could not find module with name '" <<  qs << "'.";
    }
}



/** Initialize the modules from an text filew */
void ModuleSystem::loadModules(const RevBayesCore::path &dir)
{
    if (not RevBayesCore::is_directory( dir ) )
    {
//        throw RbException() << "Cannot find directory containing modules. No modules are available. Path = " <<  dir;
    }
    else
    {
        // get the files contained in the directory
    
        // gather all text files in dir, filtered by '.ext'
        std::vector<RevBayesCore::path> filenames;

        RevBayesCore::setStringWithNamesOfFilesInDirectory( dir, filenames );
        for (auto& filename: filenames)
        {
            if ( filename.extension() == ".Rev")
            {
                Module m = Module( filename.string() );
                std::string name = filename.stem().string();
                modules.insert( std::pair<std::string, Module>(name,m) );
            }
        }
    }
}


bool ModuleSystem::isModuleAvailable(const std::string &query)
{
    // test if we have a help entry for this query string
    return modules.find( query ) != modules.end();
}


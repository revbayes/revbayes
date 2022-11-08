#ifndef ModuleSystem_H
#define ModuleSystem_H

#include "Module.h"
#include "RbFileManager.h"

#include <set>
#include <string>
#include <map>

namespace RevLanguage {
    
    
    /**
     * \brief The class managing our module system.
     *
     * Our module system consists of several modules-files. Here we load in the files.
     * The module system will provide access when a module is called.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2014-12-04
     *
     */
    class ModuleSystem {
        
    public:
        
        virtual                                    ~ModuleSystem(void);
        
        const Module&                               getModule(const std::string &qs) const;     //!< Format the help information for printing to the terminal
        bool                                        isModuleAvailable(const std::string &q);    //!< Is a module available with the query name
        void                                        loadModules(const RevBayesCore::path &dir);        //!< Initialize the help from an XML file
        
        static ModuleSystem&                        getModuleSystem(void)                       //!< Return a reference to the singleton help
        {
            static ModuleSystem globalInstance = ModuleSystem();
            return globalInstance;
        }
        
        
    private:
        ModuleSystem(void);                                                                     //!< Default constructor
        ModuleSystem(const ModuleSystem&);                                                      //!< Copy constructor (hidden away as this is a singleton class)
        ModuleSystem&                               operator=(const ModuleSystem&);             //!< Assignment operator (hidden away as this is a singleton class)
        
        
        
        std::map<std::string, Module>               modules;
        
    };
    
    
    
}

#endif


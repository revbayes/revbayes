#ifndef RbHelpDatabase_H
#define RbHelpDatabase_H

#include "RbHelpReference.h"

#include <map>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The database of strings and references for insertion into help entry objects
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     *
     */
    class RbHelpDatabase {
        
    public:
        RbHelpDatabase();
        
        // getter
        const std::string&                   getHelpString(std::string name, std::string key) { return help_strings[name][key]; }
        const std::vector<std::string>&      getHelpStringVector(std::string name, std::string key) { return help_arrays[name][key]; }
        const std::vector<RbHelpReference>&  getHelpReferences(std::string name) { return help_references[name]; }
        
        static RbHelpDatabase&               getHelpDatabase(void)
        {
            static RbHelpDatabase globalHelpInstance = RbHelpDatabase();
            return globalHelpInstance;
        }

        
    private:
        
        std::map< std::string, std::map< std::string, std::string> >                   help_strings;
        std::map< std::string, std::map< std::string, std::vector<std::string> > >     help_arrays;
        std::map< std::string, std::vector<RbHelpReference> >                          help_references;
        
    };
    
    
}

#endif


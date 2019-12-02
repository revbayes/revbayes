#ifndef RbHelpSystem_H
#define RbHelpSystem_H

#include "RbHelpDistribution.h"
#include "RbHelpEntry.h"
#include "RbHelpFunction.h"
#include "RbHelpType.h"

#include <set>
#include <string>
#include <map>

namespace RevBayesCore {
    
    
    /**
     * \brief The class managing our help system.
     *
     * Our help system consists of several xml-files. Here we load in the files.
     * The help system will provide access to other classes to the help documentation.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     *
     */
    class RbHelpSystem {
        
    public:
        
        virtual                                     ~RbHelpSystem();
        
        void                                        addHelpEntry( void );
        void                                        addHelpDistribution( RbHelpDistribution *h );
        void                                        addHelpFunction( RbHelpFunction *h );
        void                                        addHelpType( RbHelpType *h );
        const std::set<std::string>&                getFunctionEntries(void) const;
        const std::set<std::string>&                getTypeEntries(void) const;
        const RbHelpEntry&                          getHelp(const std::string &qs);                                         //!< Format the help information for printing to the terminal
        const RbHelpEntry&                          getHelp(const std::string &bq, const std::string &q);                   //!< Format the help information for printing to the terminal
        bool                                        isHelpAvailableForQuery(const std::string &q);
        bool                                        isHelpAvailableForQuery(const std::string &bq, const std::string &q);
        
        static RbHelpSystem&                        getHelpSystem(void)                                                     //!< Return a reference to the singleton help
        {
            static RbHelpSystem globalHelpInstance = RbHelpSystem();
            return globalHelpInstance;
        }
        
    private:
        RbHelpSystem(void);                                                                                                 //!< Default constructor
        RbHelpSystem(const RbHelpSystem&);                                                                                  //!< Copy constructor (hidden away as this is a singleton class)
        RbHelpSystem&                               operator=(const RbHelpSystem&);                                         //!< Assignment operator (hidden away as this is a singleton class)
        
        std::map<std::string, RbHelpFunction>                           helpForFunctions;
        std::map<std::string, std::map<std::string, RbHelpFunction> >   helpForMethods;
        std::map<std::string, RbHelpType*>                              helpForTypes;
        std::set<std::string>                                           helpFunctionNames;                                  //!< Set of finction names without aliases
        std::set<std::string>                                           helpTypeNames;                                      //!< Set of finction names without aliases
    
    };
    
    
    
    
}

#endif


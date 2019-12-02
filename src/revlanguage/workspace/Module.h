#ifndef Module_H
#define Module_H

#include <vector>
#include <string>

namespace RevLanguage {
    
    
    /**
     * \brief The module class.
     *
     * Our module system consists of several modules-files. Here we load in the files.
     * The module system will provide access when a module is called.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2014-12-04
     *
     */
    class Module {
        
    public:
                                            Module(const std::string &fn);          //!< Default constructor
        virtual                            ~Module(void);
        
        const std::vector<std::string>&     getCommandLines(void) const;
        
        
    private:

        std::vector<std::string>            commandLines;
        
        
    };
    
    
    
}

#endif


#ifndef RbHelpReference_H
#define RbHelpReference_H

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for a publication reference.
     *
     * A help entry for a publication reference is used within other help entries, for example,
     * a function call, to provide some citable source. This class is a pure data object.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     *
     */
    class RbHelpReference {
        
    public:
        RbHelpReference() {}
        
        // getter
        const std::string&                  getCitation(void) const;
        const std::string&                  getDoi(void) const;
        const std::string&                  getUrl(void) const;
        
        // setters
        void                                setCitation(const std::string& c);
        void                                setDoi(const std::string &d);
        void                                setUrl(const std::string &u);
        
        
    private:
        
        std::string                         citation;
        std::string                         doi;
        std::string                         url;
        
    };
    
    
}

#endif


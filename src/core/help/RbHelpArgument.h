#ifndef RbHelpArgument_H
#define RbHelpArgument_H

#include <string>
#include <vector>

namespace RevBayesCore {
    
    
    /**
     * \brief The help entry for an argument of a function call.
     *
     * A help entry for an argument of a function call provides the specific description for a single argument
     * in a function call. This class is a pure data object.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     *
     */
    class RbHelpArgument {
    public:
        
        RbHelpArgument() {}
        
        // getter
        const std::string&                  getArgumentDagNodeType(void) const;
        const std::string&                  getArgumentPassingMethod(void) const;
        const std::string&                  getDefaultValue(void) const;
        const std::string&                  getDescription(void) const;
        const std::string&                  getLabel(void) const;
        const std::vector<std::string>&     getOptions(void) const;
        const std::string&                  getValueType(void) const;
        
        // setters
        void                                setArgumentDagNodeType(const std::string &dt);
        void                                setArgumentPassingMethod(const std::string &pm);
        void                                setDefaultValue(const std::string &df);
        void                                setDescription(const std::string &d);
        void                                setLabel(const std::string &l);
        void                                setOptions(const std::vector<std::string> &o);
        void                                setValueType(const std::string &vt);
        
        
    private:
        
        std::string                         argumentDagNodeType;
        std::string                         argumentPassingMethod;
        std::string                         defaultValue;
        std::string                         description;
        std::string                         label;
        std::vector<std::string>            options;
        std::string                         valueType;
        
    };
    
    
}

#endif


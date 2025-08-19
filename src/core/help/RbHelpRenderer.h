#ifndef RbHelpRenderer_H
#define RbHelpRenderer_H

#include <cstddef>
#include <iosfwd>
#include <string> // IWYU pragma: keep

namespace RevBayesCore {
class RbHelpArgument;
class RbHelpEntry;
class RbHelpFunction;
class RbHelpType;
    
    /**
     * @brief Class rendering the help to the console output.
     *
     * The Help renderer converts the help entry into a nice format that can be printed directly to the console output.
     *
     * \copyright (c) Copyright 2009-2014 (GPL version 3)
     * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
     * \since Version 1.0, 2014-09-15
     */
    class HelpRenderer {
        
    public:
        
        HelpRenderer(void);
        std::string                 renderHelp(const RbHelpEntry &entryHelp, size_t w);
        std::string                 renderHelp(const RbHelpFunction &functionHelp, size_t w);
        std::string                 renderHelp(const RbHelpType &typeHelp, size_t w);
        
        
        
    private:
        
        // HelpEntry
        std::string                 renderAliases(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderAuthors(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderDescription(const RbHelpEntry &typeHelp, size_t w) const;
        std::string                 renderDetails(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderExample(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderReferences(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderSeeAlso(const RbHelpEntry &helpEntry, size_t w) const;
        std::string                 renderTitle(const RbHelpEntry &helpEntry, size_t w) const;

        // HelpFunction
        std::string                 renderArguments(const RbHelpFunction &functionHelp, size_t w) const;
        std::string                 renderArgument(const RbHelpArgument &a, size_t w, int l, const std::string &s) const;
        std::string                 renderReturnType(const RbHelpFunction &functionHelp, size_t w, bool dist = false) const;
        std::string                 renderUsage(const RbHelpFunction &functionHelp, size_t w) const;

        // HelpType
        std::string                 renderConstructors(const RbHelpType &typeHelp, size_t w) const;
        std::string                 renderMethods(const RbHelpType &typeHelp, size_t w) const;

        std::string                 section_break;
        std::string                 line_break;
        
    };
    
}

#endif


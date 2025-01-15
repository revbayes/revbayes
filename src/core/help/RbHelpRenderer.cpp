#include "RbHelpRenderer.h"

#include <string>
#include <vector>

#include "RbException.h"
#include "RbHelpDistribution.h"
#include "StringUtilities.h"
#include "TerminalFormatter.h"
#include "RbHelpArgument.h"
#include "RbHelpEntry.h"
#include "RbHelpFunction.h"
#include "RbHelpReference.h"
#include "RbHelpType.h"
#include "RbHelpDatabase.h"

using namespace RevBayesCore;

HelpRenderer::HelpRenderer(void)
{
    
    section_break = "\n\n";
    line_break = "\n";

}


std::string HelpRenderer::renderAliases(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::vector<std::string> &aliases = helpEntry.getAliases();

    std::string result = "";
    
    // aliases
    if ( aliases.size() > 0 )
    {
        result.append( TerminalFormatter::makeUnderlined("Aliases") );
        result.append( section_break );
        
        for (std::vector<std::string>::const_iterator it = aliases.begin(); it != aliases.end(); ++it)
        {
            result.append( StringUtilities::formatTabWrap(*it, 1, w) );
            result.append( line_break );
        }
        result.append( section_break );
    }
    
    return result;
}


std::string HelpRenderer::renderArguments(const RbHelpFunction &functionHelp, size_t w) const
{
    const std::vector<RbHelpArgument> &arguments = functionHelp.getArguments();

    std::string result = "";
    
    if ( arguments.size() > 0)
    {
        // check if we have multiple arguments
        if ( arguments.size() == 1 )
        {
            result.append( TerminalFormatter::makeUnderlined("Argument") );
        }
        else
        {
            result.append( TerminalFormatter::makeUnderlined("Arguments") );
        }
        result.append( section_break );
        
        const std::vector<RbHelpArgument>& args = arguments;
        int longestArgLabel = 0;
        for (std::vector<RbHelpArgument>::const_iterator it = args.begin(); it != args.end(); ++it)
        {
            if ( longestArgLabel < it->getLabel().size() )
            {
                longestArgLabel = int(it->getLabel().size());
            }
        }
        
        int indentSize = longestArgLabel + 4 + 3;
        std::string spaces = "";
        for (size_t i = 0; i < indentSize; ++i)
        {
            spaces += " ";
        }
        
        for (std::vector<RbHelpArgument>::const_iterator it = args.begin(); it != args.end(); ++it)
        {
            const RbHelpArgument& arg = *it;
            
            result.append( renderArgument( arg, w, longestArgLabel, spaces ) );
            
            //            result.append( line_break );
            
        }
        
        result.append( section_break );
    }
    
    return result;
}


std::string HelpRenderer::renderArgument(const RevBayesCore::RbHelpArgument &arg, size_t w, int longestArgLabel, const std::string &spaces) const
{
    std::string result = "";
    
    // label
    std::string l = arg.getLabel();
    StringUtilities::fillWithSpaces(l, longestArgLabel+4, false);
    std::string formattedLabel = StringUtilities::formatStringForScreen(l, "", spaces, w);
    formattedLabel = formattedLabel.substr(0, formattedLabel.size()-1);
    formattedLabel.append( " : " );
    
    // argument description
    formattedLabel.append( arg.getDescription() );
    result.append( StringUtilities::formatStringForScreen(formattedLabel, "", spaces, w) );
    //            result.append( line_break );
    
    // value type
    result.append( spaces );
    result.append( "Type:       " );
    result.append( arg.getValueType() );
    result.append( ", " );
    result.append( arg.getArgumentDagNodeType() );
    result.append( ", " );
    result.append( arg.getArgumentPassingMethod() );
    result.append( line_break );
    
    if ( arg.getDefaultValue().size() > 0 )
    {
        result.append( spaces );
        result.append( "Default:    " );
        result.append( arg.getDefaultValue() );
        result.append( line_break );
    }
    
    if (arg.getOptions().size() > 0)
    {
        result.append( spaces );
        result.append( "Options:    " );
        const std::vector<std::string> & options = arg.getOptions();
        for (std::vector<std::string>::const_iterator opt = options.begin(); opt != options.end(); ++opt)
        {
            // append a '|' separator if we had an option element before
            if ( opt != options.begin() )
            {
                result.append("|");
            }
            
            result.append( *opt );
            
        }
        result.append( line_break );
        
    }
    
    return result;
}


std::string HelpRenderer::renderAuthors(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::vector<std::string> &authors = helpEntry.getAuthor();

    std::string result = "";
    
    if ( authors.size() > 0 )
    {
        
        result.append( TerminalFormatter::makeUnderlined("Author") );
        result.append( section_break );
    
        for (std::vector<std::string>::const_iterator it = authors.begin(); it != authors.end(); ++it)
        {
            result.append( StringUtilities::formatTabWrap(*it, 1, w, false) );
            result.append( line_break );
        }
        result.append( section_break );

    }
    
    return result;
}


std::string HelpRenderer::renderConstructors(const RbHelpType &typeHelp, size_t w) const
{
    const std::vector<RbHelpFunction> &ctors = typeHelp.getConstructors();

    std::string result = "";
    
    if ( ctors.size() > 0)
    {
        // check if we have multiple arguments
        if ( ctors.size() == 1 )
        {
            if (ctors[0].getUsage() == "c_name()")
            {
                return result;
            }
            result.append( TerminalFormatter::makeUnderlined("Constructor") );
        }
        else
        {
            result.append( TerminalFormatter::makeUnderlined("Constructors") );
        }
        result.append( section_break );
        
        for (std::vector<RbHelpFunction>::const_iterator it = ctors.begin(); it != ctors.end(); ++it)
        {
            const RbHelpFunction &functionHelp = *it;
            
            result.append( StringUtilities::formatTabWrap(functionHelp.getUsage(),1,w,true) );
            result.append( line_break );
            result.append( section_break );
            
            // arguments
            result.append( renderArguments( functionHelp, w ) );
            
        }
        
    }

    
    return result;
}


std::string HelpRenderer::renderDescription(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::string &description = helpEntry.getDescription();

    std::string result = "";
    
    // descriptions
    if ( description.size() > 0 )
    {
        result.append( TerminalFormatter::makeUnderlined("Description") );
        result.append( section_break );
    
        result.append( StringUtilities::formatTabWrap(description, 1, w, false) );
        result.append( line_break );
        result.append( section_break );
    }
    
    return result;
}



std::string HelpRenderer::renderDetails(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::string &details = helpEntry.getDetails();

    std::string result = "";
    
    // details
    if ( details.size() > 0)
    {
        result.append( TerminalFormatter::makeUnderlined("Details") );
        result.append( section_break );
        
        result.append( StringUtilities::formatTabWrap(details, 1, w, false) );
        result.append( line_break );
        result.append( section_break );
    }
    
    return result;
}


std::string HelpRenderer::renderExample(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::string &example = helpEntry.getExample();

    std::string result = "";

    // example
    if ( example.size() > 0 )
    {
        result.append( TerminalFormatter::makeUnderlined("Example") );
        result.append( section_break );

        result.append( StringUtilities::formatTabWrap(example, 1, w, false) );
        result.append( section_break );
    }

    return result;
}


std::string HelpRenderer::renderHelp(const RbHelpEntry &entryHelp, size_t w)
{
    if ( dynamic_cast< const RbHelpFunction *>( &entryHelp ) != NULL)
    {
        return renderHelp( static_cast<const RbHelpFunction &>( entryHelp ), w);
    }
    else if ( dynamic_cast< const RbHelpType *>( &entryHelp ) != NULL)
    {
        return renderHelp( static_cast<const RbHelpType &>( entryHelp ), w);
    }

    throw RbException("Unexpected help type.");
}


std::string HelpRenderer::renderHelp(const RbHelpFunction &functionHelp, size_t w)
{
    
    std::string result = "";
    
    // we do not print the name at the top
    
    // title
    result.append( renderTitle( functionHelp, w ) );
    
    // aliases
    result.append( renderAliases( functionHelp, w ) );
    
    // description
    result.append( renderDescription( functionHelp, w ) );

    // usage
    result.append( renderUsage( functionHelp, w ) );
    
    // arguments
    result.append( renderArguments( functionHelp, w ) );
    
    // return type
    result.append( renderReturnType( functionHelp, w ) );
    
    // details
    result.append( renderDetails( functionHelp, w ) );
    
    // example
    result.append( renderExample( functionHelp, w ) );
    
    // author
    result.append( renderAuthors( functionHelp, w ) );
    
    // references
    result.append( renderReferences( functionHelp, w ) );
    
    // see also
    result.append( renderSeeAlso( functionHelp, w ) );
    
    
    return result;
    
}


std::string HelpRenderer::renderHelp(const RbHelpType &typeHelp, size_t w)
{
    std::string result = "";
    
    // we do not print the name at the top
    
    // title
    result.append( renderTitle( typeHelp, w ) );
    
    // aliases
    result.append( renderAliases( typeHelp, w ) );
    
    // description
    result.append( renderDescription( typeHelp, w ) );
    
    // constructors
    result.append( renderConstructors( typeHelp, w ) );

    if ( dynamic_cast< const RbHelpDistribution *>( &typeHelp ) != NULL)
    {
        // domain type
        result.append( renderReturnType( typeHelp.getConstructors().front(), w, true ) );
    }
        
    // details
    result.append( renderDetails( typeHelp, w ) );

    // example
    result.append( renderExample( typeHelp, w ) );
    
    // methods
    result.append( renderMethods( typeHelp, w ) );
    
    // authors
    result.append( renderAuthors( typeHelp, w ) );
    
    // references
    result.append( renderReferences( typeHelp, w ) );
    
    // see also
    result.append( renderSeeAlso( typeHelp, w ) );
    
    
    return result;
}


// For rendering general database entries
std::string HelpRenderer::renderHelp(RbHelpDatabase &hd, std::string entry_name, size_t w)
{
    std::string result = "";
    
    const std::string &title_prelim = hd.getHelpString(entry_name, "title");
    const std::string &description_prelim = hd.getHelpString(entry_name, "description");
    const std::string &details_prelim = hd.getHelpString(entry_name, "details");
    const std::string &example_prelim = hd.getHelpString(entry_name, "example");
    const std::vector<std::string> &authors_prelim = hd.getHelpStringVector(entry_name, "authors");
    const std::vector<std::string> &see_also_prelim = hd.getHelpStringVector(entry_name, "see_also");
    
    std::string title = "";
    std::string description = "";
    std::string details = "";
    std::string example = "";
    std::string authors = "";
    std::string see_also = "";

    if ( title_prelim.size() > 0 )
    {
        title.append( title_prelim );
        title.append( line_break );
        title.append( section_break );
    }
    
    if ( description_prelim.size() > 0 )
    {
        description.append( TerminalFormatter::makeUnderlined("Description") );
        description.append( section_break );
    
        description.append( StringUtilities::formatTabWrap(description_prelim, 1, w, false) );
        description.append( line_break );
        description.append( section_break );
    }
    
    if ( details_prelim.size() > 0)
    {
        details.append( TerminalFormatter::makeUnderlined("Details") );
        details.append( section_break );
        
        details.append( StringUtilities::formatTabWrap(details_prelim, 1, w, false) );
        details.append( line_break );
        details.append( section_break );
    }
    
    if ( example_prelim.size() > 0 )
    {
        example.append( TerminalFormatter::makeUnderlined("Example") );
        example.append( section_break );

        example.append( StringUtilities::formatTabWrap(example_prelim, 1, w, false) );
        example.append( section_break );
    }
    
    if ( authors_prelim.size() > 0 )
    {
        authors.append( TerminalFormatter::makeUnderlined("Author") );
        authors.append( section_break );
    
        for (std::vector<std::string>::const_iterator it = authors_prelim.begin(); it != authors_prelim.end(); ++it)
        {
            authors.append( StringUtilities::formatTabWrap(*it, 1, w, false) );
            authors.append( line_break );
        }
        authors.append( section_break );
    }
    
    if ( see_also_prelim.size() > 0 )
    {
        see_also.append( TerminalFormatter::makeUnderlined("See also") );
        see_also.append( section_break );
        
        for (std::vector<std::string>::const_iterator it = see_also_prelim.begin(); it != see_also_prelim.end(); ++it)
        {
            see_also.append( StringUtilities::formatTabWrap(*it, 1, w) );
            see_also.append( line_break );
        }
        see_also.append(line_break);
    }
    
    result.append( title );
    result.append( description );
    result.append( details );
    result.append( example );
    result.append( authors );
    result.append( see_also );
    
    
    return result;
}


std::string HelpRenderer::renderMethods( const RbHelpType &typeHelp, size_t w ) const
{
    const std::vector<RbHelpFunction>& methods = typeHelp.getMethods();

    const std::string &base_name = typeHelp.getName();

    std::string result = "";
    
    // methods
    if ( methods.size() > 1)
    {
        // check if we have multiple arguments
        if ( methods.size() == 1 )
        {
            result.append( TerminalFormatter::makeUnderlined("Method") );
        }
        else
        {
            result.append( TerminalFormatter::makeUnderlined("Methods") );
        }
        result.append( section_break );
        
        for (std::vector<RbHelpFunction>::const_iterator it = methods.begin(); it != methods.end(); ++it)
        {
            if ( it->getUsage() == "methods()" )
            {
                continue;
            }

            result.append( StringUtilities::formatTabWrap(it->getUsage(),1,w,true) );
            result.append( line_break );
            
            if ( it->getDescription().size() > 0 )
            {
                result.append( StringUtilities::formatTabWrap(it->getDescription(), 1, w) );
                result.append( line_break );
                result.append( section_break );
            }
        }
        
        result.append( line_break );
        
        std::string methodCall = "For more detailed information on methods type, for example, '?" + base_name + "." + methods[0].getName() + "'";
        result.append( StringUtilities::formatTabWrap( methodCall, 1, w) );
        
        result.append( section_break );
    }
    
    
    return result;
}



std::string HelpRenderer::renderReferences(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::vector<RbHelpReference> &refs = helpEntry.getReferences();

    std::string result = "";
    
    if ( refs.size() > 0 )
    {
        // check if we have multiple references
        if ( refs.size() == 1 )
        {
            result.append( TerminalFormatter::makeUnderlined("Reference") );
        }
        else
        {
            result.append( TerminalFormatter::makeUnderlined("References") );
        }
        result.append( section_break );
        
        for (std::vector<RbHelpReference>::const_iterator it = refs.begin(); it != refs.end(); ++it)
        {
            const RbHelpReference &ref = *it;
            
            if ( ref.getCitation() != "" )
            {
                result.append( StringUtilities::formatTabWrap(ref.getCitation(), 1, w) );
                result.append( line_break );
            }
            
            if ( ref.getUrl() != "" )
            {
                result.append( StringUtilities::formatTabWrap(ref.getUrl(), 1, w) );
                result.append( line_break );
            }
            
            if ( ref.getDoi() != "" )
            {
                result.append( StringUtilities::formatTabWrap(ref.getDoi(), 1, w) );
                result.append( line_break );
            }
            result.append( line_break );
            
            
        }
        
        result.append( line_break );
        
    }

    return result;
}

std::string HelpRenderer::renderReturnType(const RbHelpFunction &functionHelp, size_t w, bool dist) const
{
    const std::string &return_type = functionHelp.getReturnType();

    std::string result = "";

    // example
    if ( return_type.size() > 0 )
    {
        result.append( TerminalFormatter::makeUnderlined(dist ? "Domain Type" : "Return Type") );
        result.append( section_break );
        result.append( StringUtilities::formatTabWrap(functionHelp.getReturnType(), 1, w) );
        result.append( line_break );
        result.append( section_break );
    }

    return result;
}

std::string HelpRenderer::renderSeeAlso(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::vector<std::string> &others = helpEntry.getSeeAlso();

    std::string result = "";
    
    if ( others.size() > 0 )
    {
        result.append( TerminalFormatter::makeUnderlined("See also") );
        result.append( section_break );
        
        for (std::vector<std::string>::const_iterator it = others.begin(); it != others.end(); ++it)
        {
            result.append( StringUtilities::formatTabWrap(*it, 1, w) );
            result.append( line_break );
        }
        result.append(line_break);
        
    }

    return result;
}

std::string HelpRenderer::renderTitle(const RbHelpEntry &helpEntry, size_t w) const
{
    const std::string &title = helpEntry.getTitle();

    std::string result = "";

    if ( title.size() > 0 )
    {
        result.append( title );
        result.append( line_break );
        result.append( section_break );
    }

    return result;
}

std::string HelpRenderer::renderUsage(const RbHelpFunction &functionHelp, size_t w) const
{
    const std::string &usage = functionHelp.getUsage();

    std::string result = "";
    
    if (usage.size() > 0)
    {
        result.append( TerminalFormatter::makeUnderlined("Usage") );
        result.append( section_break );
        
        result.append( StringUtilities::formatTabWrap( usage, 1, w, true) );
        result.append( line_break );
        result.append( section_break );
    }
    
    return result;
}

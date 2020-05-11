#include "YAMLHelpRenderer.h"

#include "RbHelpArgument.h"
#include "RbHelpDistribution.h"
#include "RbHelpFunction.h"
#include "RbHelpReference.h"
#include "StringUtilities.h"

using namespace YAML;

std::string HelpRenderer::renderHelp(const RevBayesCore::RbHelpType &help)
{
    Node node( getAttributes(help) );

    return node.render();
}

std::string HelpRenderer::renderHelp(const RevBayesCore::RbHelpFunction &help)
{
    Node node( getAttributes(help) );

    return node.render();
}

std::vector<Element*> HelpRenderer::getAttributes(const RevBayesCore::RbHelpType &help)
{
    std::vector<Element*> yaml;

    yaml.push_back( new Attribute("name", help.getName() ) );
    yaml.push_back( new Attribute("title", help.getTitle() ) );
    yaml.push_back( new Attribute("aliases", help.getAliases() ) );
    yaml.push_back( new Attribute("description", help.getDescription() ) );
    yaml.push_back( new Attribute("details", help.getDetails() ) );
    yaml.push_back( new Attribute("example", help.getExample() ) );
    yaml.push_back( new Attribute("authors", help.getAuthor() ) );
    yaml.push_back( new Attribute("see_also", help.getSeeAlso() ) );
    yaml.push_back( new Attribute("type_spec", help.getTypeSpec() ) );

    Attribute* constructors = new Attribute("constructor");
    for( size_t i = 0; i < help.getConstructors().size(); i++)
    {
        if ( help.getConstructors()[i].getUsage() == "c_name()" )
        {
            continue;
        }
        if ( i == 0 && dynamic_cast<const RevBayesCore::RbHelpDistribution*>( &help ) != NULL )
        {
            yaml.push_back( new Attribute("domain", help.getConstructors()[i].getReturnType() ) );
        }
        constructors->children.push_back( new Node( getAttributes(help.getConstructors()[i], CONSTRUCTOR) ) );
    }
    yaml.push_back(constructors);

    Attribute* methods = new Attribute("methods");
    for( size_t i = 0; i < help.getMethods().size(); i++)
    {
        if ( help.getMethods()[i].getName() == "methods" )
        {
            continue;
        }
        methods->children.push_back( new Node( getAttributes(help.getMethods()[i], METHOD) ) );
    }
    yaml.push_back(methods);

    Attribute* refs = new Attribute("references");
    for( size_t i = 0; i < help.getReferences().size(); i++)
    {
        refs->children.push_back( new Node( getAttributes(help.getReferences()[i]) ) );
    }
    yaml.push_back(refs);

    return yaml;
}

std::vector<Element*> HelpRenderer::getAttributes(const RevBayesCore::RbHelpFunction &help, FunctionType type )
{
    std::vector<Element*> yaml;

    if ( type != CONSTRUCTOR )
    {
        yaml.push_back( new Attribute("name", help.getName() ) );
        yaml.push_back( new Attribute("return_type", help.getReturnType() ) );
    }

    yaml.push_back( new Attribute("usage", help.getUsage() ) );

    Attribute* args = new Attribute("arguments");
    for( size_t i = 0; i < help.getArguments().size(); i++)
    {
        args->children.push_back( new Node( getAttributes(help.getArguments()[i]) ) );
    }
    yaml.push_back(args);

    if ( type == FUNCTION )
    {
        yaml.push_back( new Attribute("title", help.getTitle() ) );
        yaml.push_back( new Attribute("aliases", help.getAliases() ) );
        yaml.push_back( new Attribute("description", help.getDescription() ) );
        yaml.push_back( new Attribute("details", help.getDetails() ) );
        yaml.push_back( new Attribute("example", help.getExample() ) );
        yaml.push_back( new Attribute("authors", help.getAuthor() ) );
        yaml.push_back( new Attribute("see_also", help.getSeeAlso() ) );
        yaml.push_back( new Attribute("type_spec", help.getTypeSpec() ) );

        Attribute* refs = new Attribute("references");
        for( size_t i = 0; i < help.getReferences().size(); i++)
        {
            refs->children.push_back( new Node( getAttributes(help.getReferences()[i]) ) );
        }
        yaml.push_back(refs);
    }

    return yaml;
}

std::vector<Element*> HelpRenderer::getAttributes(const RevBayesCore::RbHelpArgument &help)
{
    std::vector<Element*> yaml;

    yaml.push_back( new Attribute("label", help.getLabel() ) );
    yaml.push_back( new Attribute("description", help.getDescription() ) );
    yaml.push_back( new Attribute("dag_type", help.getArgumentDagNodeType() ) );
    yaml.push_back( new Attribute("pass_by", help.getArgumentPassingMethod() ) );
    yaml.push_back( new Attribute("value_type", help.getValueType() ) );
    yaml.push_back( new Attribute("default", help.getDefaultValue() ) );
    yaml.push_back( new Attribute("options", help.getOptions() ) );

    return yaml;
}

std::vector<Element*> HelpRenderer::getAttributes(const RevBayesCore::RbHelpReference &help)
{
    std::vector<Element*> yaml;

    yaml.push_back( new Attribute("citation", help.getCitation() ) );
    yaml.push_back( new Attribute("doi", help.getDoi() ) );
    yaml.push_back( new Attribute("url", help.getUrl() ) );

    return yaml;
}

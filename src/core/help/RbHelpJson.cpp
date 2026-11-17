#include "RbHelpJson.h"

#include <set>
#include <string>
#include <vector>

#include "RbHelpArgument.h"
#include "RbHelpDatabase.h"
#include "RbHelpDistribution.h"
#include "RbHelpFunction.h"
#include "RbHelpReference.h"
#include "RbHelpSystem.h"
#include "RbHelpType.h"
#include "StringUtilities.h"
#include "nlohmann-json.h"

using nlohmann::ordered_json;
using namespace RevBayesCore;

namespace
{
    enum FunctionKind { FUNCTION, METHOD, CONSTRUCTOR };

    std::string trimRight(const std::string& s)
    {
        size_t end = s.find_last_not_of(" \n\r\t");
        return (end == std::string::npos) ? "" : s.substr(0, end + 1);
    }

    // Empty strings and lists are null.
    ordered_json value(const std::string& s)
    {
        return s.empty() ? ordered_json() : ordered_json(trimRight(s));
    }

    ordered_json value(const std::vector<std::string>& v)
    {
        if (v.empty()) return ordered_json();
        ordered_json a = ordered_json::array();
        for (const std::string& s : v) a.push_back(trimRight(s));
        return a;
    }

    ordered_json value(const ordered_json& a)
    {
        return a.empty() ? ordered_json() : a;
    }

    ordered_json toJson(const RbHelpReference& ref)
    {
        ordered_json j;
        j["citation"] = value(ref.getCitation());
        j["doi"]      = value(ref.getDoi());
        j["url"]      = value(ref.getUrl());
        return j;
    }

    ordered_json toJson(const std::vector<RbHelpReference>& refs)
    {
        ordered_json a = ordered_json::array();
        for (const RbHelpReference& ref : refs) a.push_back(toJson(ref));
        return value(a);
    }

    ordered_json toJson(const RbHelpArgument& arg)
    {
        ordered_json j;
        j["label"]       = value(arg.getLabel());
        j["description"] = value(arg.getDescription());
        j["dag_type"]    = value(arg.getArgumentDagNodeType());
        j["pass_by"]     = value(arg.getArgumentPassingMethod());
        j["value_type"]  = value(arg.getValueType());
        j["default"]     = value(arg.getDefaultValue());
        j["options"]     = value(arg.getOptions());
        return j;
    }

    ordered_json toJson(const RbHelpFunction& help, FunctionKind kind)
    {
        ordered_json j;
        if (kind != CONSTRUCTOR)
        {
            j["name"]        = value(help.getName());
            j["return_type"] = value(help.getReturnType());
        }
        j["usage"] = value(help.getUsage());

        ordered_json args = ordered_json::array();
        for (const RbHelpArgument& arg : help.getArguments()) args.push_back(toJson(arg));
        j["arguments"] = value(args);

        if (kind == FUNCTION)
        {
            j["title"]       = value(help.getTitle());
            j["aliases"]     = value(help.getAliases());
            j["description"] = value(help.getDescription());
            j["details"]     = value(help.getDetails());
            j["example"]     = value(help.getExample());
            j["authors"]     = value(help.getAuthor());
            j["see_also"]    = value(help.getSeeAlso());
            j["type_spec"]   = value(help.getTypeSpec());
            j["references"]  = toJson(help.getReferences());
        }
        return j;
    }

    ordered_json toJson(const RbHelpType& help)
    {
        ordered_json j;
        j["name"]        = value(help.getName());
        j["title"]       = value(help.getTitle());
        j["aliases"]     = value(help.getAliases());
        j["description"] = value(help.getDescription());
        j["details"]     = value(help.getDetails());
        j["example"]     = value(help.getExample());
        j["authors"]     = value(help.getAuthor());
        j["see_also"]    = value(help.getSeeAlso());
        j["type_spec"]   = value(help.getTypeSpec());

        ordered_json constructors = ordered_json::array();
        const std::vector<RbHelpFunction>& ctors = help.getConstructors();
        for (size_t i = 0; i < ctors.size(); i++)
        {
            if (ctors[i].getUsage() == "c_name()") continue;
            if (i == 0 && dynamic_cast<const RbHelpDistribution*>(&help) != nullptr)
            {
                j["domain"] = value(ctors[i].getReturnType());
            }
            constructors.push_back(toJson(ctors[i], CONSTRUCTOR));
        }
        j["constructor"] = value(constructors);

        ordered_json methods = ordered_json::array();
        for (const RbHelpFunction& method : help.getMethods())
        {
            if (method.getName() == "methods") continue;
            methods.push_back(toJson(method, METHOD));
        }
        j["methods"] = value(methods);

        j["references"] = toJson(help.getReferences());
        return j;
    }
}


std::string RevBayesCore::renderHelpJson(void)
{
    RbHelpSystem& help = RbHelpSystem::getHelpSystem();
    const std::set<std::string>& function_names = help.getFunctionEntries();
    const std::set<std::string>& type_names = help.getTypeEntries();

    ordered_json entries = ordered_json::array();

    for (const std::string& name : function_names)
    {
        if (name.empty()) continue;

        // d/p/q/r functions are documented with their distribution
        std::string base = name.substr(1);
        if (type_names.count("dn" + base) || type_names.count("dn" + StringUtilities::firstCharToUpper(base))) continue;

        const RbHelpFunction& entry = static_cast<const RbHelpFunction&>(help.getHelp(name));
        if (entry.getUsage() == "c_name()") continue;

        entries.push_back(toJson(entry, FUNCTION));
    }

    for (const std::string& name : type_names)
    {
        // vector types are documented with their element type
        if (name.find("[]") != std::string::npos) continue;
        if (name.empty() || name == "c_name" || function_names.count(name)) continue;

        const RbHelpType& entry = static_cast<const RbHelpType&>(help.getHelp(name));
        entries.push_back(toJson(entry));
    }

    return entries.dump(2) + "\n";
}


std::vector<std::string> RevBayesCore::checkHelpPages(void)
{
    return RbHelpDatabase::getHelpDatabase().getProblems();
}

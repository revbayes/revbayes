#include "RbHelpDatabase.h"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

using namespace RevBayesCore;

namespace
{
    struct HelpPage
    {
        const char* name;
        const char* text;
    };

    // Generated at build time from help/md/*.md by help/embed_help.cmake or help/embed_help.py
    const HelpPage help_pages[] =
    {
#include "RbHelpText.inc"
    };

    bool isSpace(char c)
    {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
    }

    std::string trimRight(const std::string& s)
    {
        size_t end = s.size();
        while (end > 0 && isSpace(s[end-1])) end--;
        return s.substr(0, end);
    }

    std::string trim(const std::string& s)
    {
        size_t begin = 0;
        while (begin < s.size() && isSpace(s[begin])) begin++;
        return trimRight(s.substr(begin));
    }

    // Split on sep, dropping trailing empty fields.
    std::vector<std::string> split(const std::string& s, const std::string& sep)
    {
        std::vector<std::string> fields;
        size_t start = 0;
        while (true)
        {
            size_t pos = s.find(sep, start);
            if (pos == std::string::npos)
            {
                fields.push_back(s.substr(start));
                break;
            }
            fields.push_back(s.substr(start, pos - start));
            start = pos + sep.size();
        }
        while (!fields.empty() && fields.back().empty()) fields.pop_back();
        return fields;
    }

    std::vector<std::string> parseList(const std::string& value)
    {
        std::vector<std::string> items;
        for (const std::string& line : split(value, "\n"))
        {
            items.push_back(trim(line));
        }
        return items;
    }

    // Each reference is a YAML-like item with citation, doi and url keys, whose values may wrap across lines.
    RbHelpReference parseReference(const std::string& item)
    {
        // drop leading whitespace on every line, including blank lines
        std::string a;
        for (size_t i = 0; i < item.size(); i++)
        {
            if (i == 0 || item[i-1] == '\n')
            {
                while (i < item.size() && isSpace(item[i])) i++;
                if (i == item.size()) break;
            }
            a += item[i];
        }

        // join wrapped lines, along with any backslash continuations
        std::string b;
        for (size_t i = 0; i < a.size(); i++)
        {
            size_t j = i;
            if (a[j] == '\\' && j+1 < a.size() && a[j+1] == '\n') j++;
            if (a[j] != '\n')
            {
                b += a[i];
                continue;
            }
            j++;
            while (j < a.size() && isSpace(a[j])) j++;
            if (j < a.size() && a[j] == '\\') j++;
            b += ' ';
            i = j - 1;
        }

        // put each key on its own line
        std::string c;
        for (size_t i = 0; i < b.size(); i++)
        {
            bool before_key = b[i] == ' ' && (b.compare(i+1, 8, "citation") == 0 || b.compare(i+1, 3, "doi") == 0 || b.compare(i+1, 3, "url") == 0);
            c += before_key ? '\n' : b[i];
        }

        std::map<std::string, std::string> fields;
        for (const std::string& line : split(c, "\n"))
        {
            size_t first = std::string::npos;
            std::string key;
            for (const std::string k : {"citation", "doi", "url"})
            {
                size_t pos = line.find(k + ": ");
                if (pos < first)
                {
                    first = pos;
                    key = k;
                }
            }
            if (first == std::string::npos) continue;

            std::string value = line.substr(first + key.size() + 2);
            size_t null = value.find("null");
            if (null != std::string::npos)
            {
                size_t end = null + 4;
                while (end < value.size() && isSpace(value[end])) end++;
                value.erase(null, end - null);
            }
            fields[key] = value;
        }

        return RbHelpReference(fields["citation"], fields["doi"], fields["url"]);
    }

    std::vector<RbHelpReference> parseReferences(std::string value)
    {
        size_t first = 0;
        while (first < value.size() && isSpace(value[first])) first++;
        if (value.compare(first, 2, "- ") == 0) value.erase(0, first + 2);

        std::vector<RbHelpReference> refs;
        for (const std::string& item : split(value, "- "))
        {
            refs.push_back(parseReference(item));
        }
        return refs;
    }

    // Remove the indentation common to all non-blank lines.
    std::string dedent(const std::string& value)
    {
        std::vector<std::string> lines = split(value, "\n");

        size_t indent = std::string::npos;
        for (const std::string& line : lines)
        {
            size_t spaces = line.find_first_not_of(' ');
            if (spaces != std::string::npos && spaces < indent) indent = spaces;
        }
        if (indent == std::string::npos) return value;

        std::string result;
        for (size_t i = 0; i < lines.size(); i++)
        {
            if (i > 0) result += '\n';
            if (indent < lines[i].size()) result += trimRight(lines[i].substr(indent));
        }
        return result;
    }

    std::string parseString(const std::string& key, const std::string& value)
    {
        std::string s;
        for (char c : value)
        {
            if (c == '\t') s += "        ";
            else s += c;
        }
        if (key == "example") s = dedent(s);
        return trimRight(s);
    }

    // Split a page into its "## key" sections.
    std::map<std::string, std::string> parseSections(const std::string& text)
    {
        std::map<std::string, std::string> sections;
        std::string key;
        std::string value;
        bool in_section = false;

        size_t start = 0;
        while (start < text.size())
        {
            size_t end = text.find('\n', start);
            end = (end == std::string::npos) ? text.size() : end + 1;
            std::string line = text.substr(start, end - start);
            start = end;

            if (line.compare(0, 3, "## ") == 0)
            {
                if (in_section && !key.empty()) sections[key] = value;
                key = line.substr(3);
                if (!key.empty() && key.back() == '\n') key.pop_back();
                value = "";
                in_section = true;
            }
            else
            {
                value += line;
            }
        }
        if (in_section) sections[key] = value;

        return sections;
    }
}


RbHelpDatabase::RbHelpDatabase()
{
    for (const HelpPage& page : help_pages)
    {
        for (const auto& [key, value] : parseSections(page.text))
        {
            if (key == "references")
            {
                for (const RbHelpReference& ref : parseReferences(value))
                {
                    help_references[page.name].push_back(ref);
                }
            }
            else if (key == "see_also" || key == "authors")
            {
                for (const std::string& item : parseList(value))
                {
                    help_arrays[page.name][key].push_back(item);
                }
            }
            else
            {
                std::string s = parseString(key, value);
                if (!s.empty()) help_strings[page.name][key] = s;
            }
        }
    }
}

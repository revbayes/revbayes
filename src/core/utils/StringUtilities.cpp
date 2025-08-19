/**
 * @file
 * This file contains the implementation of StringUtilities, which 
 * helper functions for manipulating strings.
 *
 * @brief Declaration of StringUtilities
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since version 1.0 2009-09-02
 *
 * $Id$
 */

#include "StringUtilities.h"

#include <cstdio>
#include <cstdint>
#include <iomanip>

#include <algorithm>
#include <string>
#include <cstdlib>

#include "RbFileManager.h"
#include "RbVector.h"


using std::string;
using std::vector;


/** Convert the string s to a number */
int StringUtilities::asIntegerNumber(const std::string& s)
{
    
    return std::atoi( s.c_str() );
}

/**
 * Fill this string with spaces so that it has the required length.
 * Either fill the spaces on the right if left aligned (true)
 * or on the left if right aligned.
 */
void StringUtilities::fillWithSpaces(std::string &s, int l, bool left)
{
    
    for (int i=int(s.length()); i<l ; ++i)
    {
        // either left algined
        if (left == true)
        {
            s.push_back(' ');
        }
        else
        {
            s.insert(s.begin(), ' ');
        }
        
    }
    
}


/**
 * Find the first occurence of the given character.
 * We return string::npos if it wasn't found.
 */
size_t StringUtilities::findFirstOf(const std::string &s, char c)
{
    size_t pos = std::string::npos;
    
    for (size_t i=0; i<s.length(); ++i)
    {
        if ( s[i] == c )
        {
            pos = i;
            break;
        }
    }
    
    return pos;
}


/**
 * Find the first occurence of the given character.
 * We return string::npos if it wasn't found.
 */
size_t StringUtilities::findFirstOf(const std::string &a, const std::string &b)
{
    size_t pos = a.find(b);
    
    return pos;
}


/**
 * Find the last occurence of the given character.
 * We return string::npos if it wasn't found.
 */
size_t StringUtilities::findLastOf(const std::string &s, char c)
{
    size_t pos = std::string::npos;
    
    for (size_t i=s.length(); i>0; --i)
    {
        if ( s[i-1] == c )
        {
            pos = i-1;
            break;
        }
    }
    
    return pos;
}

/**
 * Fill this string with spaces so that it has the required length.
 * Either fill the spaces on the right if left aligned (true)
 * or on the left if right aligned.
 */
void StringUtilities::formatFixedWidth(std::string &s, int l, bool left)
{
    
    if ( s.length() > l )
    {
        if ( l > 2)
        {
            s = s.substr(0, l-2);
        }
        
        s += "..";
    }
    
    fillWithSpaces(s, l, left);
}


/**
 * Wraps text so that each line doesn't exceeds column width.
 * Prepends tabs to each line.
 *
 * @param s             string to format
 * @param tabs          number of tabs to prepend to each line
 * @param width         column width
 * @param removeFormat  strips newline, tab etc if set to true
 * @return              formatted text
 */
std::string StringUtilities::formatTabWrap(std::string s, size_t tabs, size_t width, bool removeFormat)
{
    
    std::string tabbing = "";
    for (int i = 0; i < tabs * 4; i++)
    {
        tabbing.append(" ");
    }
    
    std::string result = tabbing;
    size_t w = width - tabbing.size();
    int cc = 0; // character count
    char lastChar = '\0';
    
    for (unsigned int i = 0; i < s.size(); i++)
    {
        if ( result.size() > 0 )
        {
            lastChar = result[result.size() - 1];
        }
        
        // strip consecutive spaces
        if ( removeFormat == false )
        {
            result += s[i];
            lastChar = s[i];
            cc++;
        }
        else if ( !StringUtilities::isFormattingChar(s[i]) )
        {
            result += s[i];
            lastChar = s[i];
            cc++;
        }
        
        if( i < s.size() - 1 )
        {
            if (lastChar == '\n')
            {
                cc = 0;
                result.append(tabbing);
            }
            else if (lastChar == ' ')
            {
                // we now have a possible point where to wrap the line.
                // peek ahead and see where next possible wrap point is:
                std::string sub_str = s.substr(i+1);
                size_t next = StringUtilities::findFirstOf(sub_str, ' ');

                // if next wrap point is beyond the width, then wrap line now
                if (cc + next >= w)
                {
                    result.append("\n").append(tabbing);
                    // reset char count for next line
                    cc = 0;
                }

            }
        }
        
    }
    
    return result;
}


/** Format string for printing to screen, with word wrapping, and various indents */
std::string StringUtilities::formatStringForScreen(const std::string &s, const std::string &firstLinePad, const std::string &hangingPad, size_t screenWidth)
{

    std::string outputString = "";

    std::vector<std::string> lineList = std::vector<std::string>();
    std::string del = "\n";
    StringUtilities::stringSplit( s, del, lineList );

    for ( size_t i=0; i<lineList.size(); i++ )
    {
    
        std::vector<std::string> stringList = std::vector<std::string>();
        std::string space = " ";
        std::string line = lineList[i];
        StringUtilities::stringSplit(line, space, stringList);

        if ( stringList.size() > 0 )
        {
            outputString += firstLinePad;

            size_t count = firstLinePad.size();
            for ( std::vector<std::string>::iterator it = stringList.begin(); it != stringList.end(); it++ )
            {

                if ( count + (*it).size() > screenWidth && count != 0 )
                {
                    outputString[outputString.size()-1] = '\n';
                    outputString += hangingPad;
                    count = hangingPad.size();
                }
                count += (*it).size() + 1;
                outputString += (*it) + " ";
            }
            outputString[outputString.size()-1] = '\n';
        }
        else
        {
            outputString += "\n";
        }
        
    }

    return outputString;
}


/**
 * Indicates if a char is affecting text formatting
 * @param c
 * @return
 */
bool StringUtilities::isFormattingChar(char c)
{
    return c == '\n' || c == '\t' || c == '\r';
}


/** Determine if the string s represents a number */
bool StringUtilities::isIntegerNumber(const std::string& s)
{
    
    bool exponent = false;
    bool sign = false;
    bool digit = false;

    for (size_t i=0; i<s.size(); i++)
    {
        if ( isdigit(s[i]) )
        {
            digit = true;
        }
        else if(s[i] == '.')
        {
            return false;
        }
        else if(s[i] == 'e')
        {
            if( exponent || !digit ) return false;

            exponent = true;

            sign = false;
            digit = false;
        }
        else if(s[i] == '+')
        {
            if( sign || digit ) return false;

            sign = true;
        }
        else if(s[i] == '-')
        {
            if( sign || digit || exponent ) return false;

            sign = true;
        }
        else
        {
            return false;
        }

    }
    
    return true;
}


/** Determine if the string s represents a number */
bool StringUtilities::isNumber(const std::string& s)
{
    bool exponent = false;
    bool sign = false;
    bool decimal = false;
    bool digit = false;

    for (size_t i=0; i<s.size(); i++)
    {
        if ( isdigit(s[i]) )
        {
            digit = true;
        }
        else if(s[i] == '.')
        {
            if( decimal ) return false;

            decimal = true;
        }
        else if(s[i] == 'e')
        {
            if( exponent || !digit ) return false;

            exponent = true;

            sign = false;
            decimal = false;
            digit = false;
        }
        else if(s[i] == '+' || s[i] == '-')
        {
            if( sign || digit || decimal ) return false;

            sign = true;
        }
        else
        {
            return false;
        }
        
    }
    
    return true;
}


/**
 * Utility function for getting a one-line summary being max maxLen std::int64_t.
 * We find the first non-empty line in the input. If it is longer than maxLen,
 * we truncate it at maxLen - 3 and add "..." at the end. If it is shorter, we
 * just return the complete line (without line break).
 */
std::string StringUtilities::oneLiner( const std::string& input, size_t maxLen )
{
    std::string oneLiner;
    
    std::ostringstream x;
    size_t begin = 0;
    size_t i = 0;
    for ( size_t i = 0; i < input.size(); ++i )
    {
        if ( input[ i ] == '\n' || input[ i ] == '\r' )
        {
            begin = i + 1;
        }
        
        if ( isgraph( input[ i ] ) )
        {
            break;
        }
        
    }
    size_t firstGraph = i;
    
    // Empty input; return empty oneliner
    if ( firstGraph >= input.size() )
    {
        return oneLiner;
    }

    // Construct the oneliner
    for ( i = begin; i < input.size() && i < maxLen; ++i )
    {
        if ( input[ i ] == '\n' || input[ i ] == '\r' )
        {
            break;
        }
        oneLiner.push_back( input[ i ] );
    }

    // If the oneliner was cut short by a newline, check whether there
    // are more printable graphs in the input. If so, we append dots.
    if ( input[i] == '\n' || input[i] == '\r' )
    {
        for ( ; i < input.size(); ++i )
        {
            if ( isgraph( input[i] ) )
                break;
        }
        
        if ( i < input.size() )
        {
            if ( maxLen - oneLiner.size() < 3 )
            {
                oneLiner[ maxLen - 1 ] = '.';
                oneLiner[ maxLen - 2 ] = '.';
                oneLiner[ maxLen - 3 ] = '.';
            }
            else
            {
                oneLiner += "...";
            }
        }
    }
    // If not, we check if the oneliner was cut short by the maximum length,
    // and, if so, insert dots into the last three character positions
    else if ( oneLiner.size() == maxLen )
    {
        oneLiner[ maxLen - 1 ] = '.';
        oneLiner[ maxLen - 2 ] = '.';
        oneLiner[ maxLen - 3 ] = '.';
    }

    return oneLiner;
}


void StringUtilities::replaceSubstring(std::string& str, const std::string& oldStr, const std::string& newStr)
{
    size_t pos = 0;
    while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
    
}


void StringUtilities::replaceAllOccurrences(std::string& str, char old_ch, char new_ch)
{
    
    for (size_t i=0; i<str.size(); ++i)
    {
        if ( str[i] == old_ch )
        {
            str[i] = new_ch;
        }
    }
    
}

void StringUtilities::join(std::ostream& o, const vector<string>& ss, const string& sep)
{
    for(int i=0;i<ss.size();i++)
    {
        o<<ss[i];
        if (i+1 < ss.size())
            o<<sep;
    }
}

string StringUtilities::join(const vector<string>& ss, const string& sep)
{
    std::ostringstream o;
    join(o, ss, sep);
    return o.str();
}

/*!
 * Utility function for index sorting, inspired by the following Stack Overflow posts:
 * Lukasz Wiklendt, https://stackoverflow.com/a/12399290 and
 * vsoftco, https://stackoverflow.com/a/37732329
 *
 * \brief Index sorting.
 * \param v is a vector that we want to sort and return the indices of sorted elements.
 * \return Vector of indices of the sorted elements of v. We use uint32_t rather than size_t to speed up memory access.
 * \throws Does not throw an error.
 */
std::vector<uint32_t> StringUtilities::stringSortIndices(const std::vector<std::string>& v)
{
    std::vector<uint32_t> result( v.size() );
    std::iota(result.begin(), result.end(), 0);
    std::sort(result.begin(), result.end(),
              [&v](uint32_t i1, uint32_t i2) { return v[i1] < v[i2]; }
             );
    return result;
}

/**
 * Utility function for dividing string into pieces
 * If no delimiter is specified, then the string is split on whitespace.
 */
void StringUtilities::stringSplit(std::string str, std::string delim, std::vector<std::string>& results, bool trim)
{

    std::string::iterator cut;

    // if we are delimiting on whitespace*,
    // then assume the first and last fields are non-empty
    // i.e. remove leading and trailing whitespace
    if ( delim.empty() )
    {
        // erase trailing whitespace
        str.erase(std::find_if (str.rbegin(), str.rend(), [](char c) {return not isspace(c);} ).base(), str.end());
        // erase leading whitespace
        str.erase(str.begin(), std::find_if (str.begin(), str.end(), [](char c) {return not isspace(c);} ));
    }

    while ( true )
    {
        if ( delim.empty() )
        {
            // find first whitespace character
            cut = std::find_if(str.begin(), str.end(), ::isspace);
        }
        else
        {
            // find the first occurrence of the full delimiter string
            cut = std::search(str.begin(), str.end(), delim.begin(), delim.end());
        }

        std::string substr(str.begin(), cut);

        if ( trim )
        {
            // erase trailing whitespace
            substr.erase(std::find_if (substr.rbegin(), substr.rend(), [](char c) {return not isspace(c);}).base(), substr.end());
            // erase leading whitespace
            substr.erase(substr.begin(), std::find_if (substr.begin(), substr.end(), [](char c) {return not isspace(c);}));
        }

        results.push_back(substr);

        if ( cut == str.end() )
        {
            break;
        }

        str = std::string(cut + delim.size(), str.end());

        if ( delim.empty() )
        {
            // erase leading whitespace in remaining string
            str.erase(str.begin(), std::find_if (str.begin(), str.end(), [](char c) {return not isspace(c);}));
        }
    }
}


/** Utility function for converting string to all lower case */
void StringUtilities::toLower(std::string& str)
{

    for (size_t i=0; i<str.size(); i++)
    {
        str[i] = char( tolower(str[i]) );
    }
    
}

std::string StringUtilities::toString(double x, int digits)
{
    
    std::stringstream o;
//    o << std::setw(12) << std::setprecision(digits);
    o << std::setprecision(digits);
    o << x;
    return o.str();

}


/** Utility function for converting string to all lower case */
std::string& StringUtilities::firstCharToUpper(std::string& str)
{
    
    if ( str.size() > 0)
    {
        str[0] = char( toupper(str[0]) );
    }
    
    return str;
}

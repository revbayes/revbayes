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

#include <stdio.h>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "RbFileManager.h"
#include "RbVector.h"





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


std::string StringUtilities::getStringWithDeletedLastPathComponent(const std::string& s)
{
    
#	ifdef RB_WIN
    std::string pathSeparator = "\\";
#	else
    std::string pathSeparator = "/";
#   endif
    
    std::string tempS = s;
	size_t location = StringUtilities::findLastOf(tempS, pathSeparator[0]);
	if ( location == std::string::npos )
    {
		/* There is no path in this string. We
         must have only the file name. */
		return "";
    }
	else if ( location == tempS.length() - 1 )
    {
		/* It looks like the last character is "/", which
         means that no file name has been provided. */
		return tempS;
    }
	else
    {
		/* We can divide the path into the path and the file. */
		tempS.erase( location );
		return tempS;
    }
    
    return "";
}


/** Return file contents as string with '\n' line breaks */
std::string StringUtilities::getFileContentsAsString(const std::string& s)
{

    RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(s);
    
    // open file
	std::ifstream fStrm;
    fStrm.open(fm.getFullFileName().c_str(), std::ios::in);
    if ( !fStrm.is_open() )
        return "";
        
    // read the file
    int ch;
    std::string retStr = "";
    while ( (ch = fStrm.get()) != EOF)
    {
        char c = (char)ch;
        
        if (ch == '\n' || ch == '\r' || ch == EOF)
            retStr += '\n';
        else
            retStr += c;
        } 

    // close file
    fStrm.close();

    return retStr;
}


/** Find the last component of a file path */
std::string StringUtilities::getLastPathComponent(const std::string& s)
{
#	ifdef RB_WIN
    std::string pathSeparator = "\\";
#	else
    std::string pathSeparator = "/";
#   endif

    std::vector<std::string> sVec;
    stringSplit(s, pathSeparator, sVec);
    std::string lastComponent = "";
    if (sVec.size() > 0)
    {
        lastComponent = sVec[sVec.size()-1];
    }
    
    return lastComponent;
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
 * Utility function for getting a one-line summary being max maxLen long.
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


/** Utility function for dividing string into pieces */
void StringUtilities::stringSplit(const std::string &s, const std::string &delim, std::vector<std::string>& results)
{
    
    // create our own copy of the string
    std::string str = s;

    size_t cutAt;
    while ( (cutAt = StringUtilities::findFirstOf(str, delim)) != str.npos )
    {
        if (cutAt > 0)
        {
            results.push_back(str.substr(0, cutAt));
        }
        else
        {
            results.push_back( "" );
        }
        str = str.substr(cutAt+delim.size());
    }
    
    if (str.length() > 0)
    {
        results.push_back(str);
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


/** Global operator for appending double to std::string */
std::string operator+( const std::string& A, double B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( const std::string& A, int B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( const std::string& A, long B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( const std::string& A, size_t B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( double A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( int A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( long A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string operator+( size_t A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( const std::string& A, double B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( const std::string& A, int B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( const std::string& A, long B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( const std::string& A, size_t B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( double A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( int A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( long A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevLanguage::operator+( size_t A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, double B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, int B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, long B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, size_t B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, const RbVector<double> &B )
{
    RbVector<double> C = B;
    std::stringstream o;
    o << A << C;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( const std::string& A, const RbVector<long> &B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( double A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( int A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( long A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}


/** Global operator for appending double to std::string */
std::string RevBayesCore::operator+( size_t A, const std::string& B )
{
    
    std::stringstream o;
    o << A << B;
    return o.str();
}





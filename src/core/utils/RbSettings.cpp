#include "RbSettings.h"

#include <unistd.h>
#include <cstdlib> //includes std::atof
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include <optional>

#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"

#	ifdef _WIN32
#include <windows.h>
#   endif

using namespace RevBayesCore;


/** Default constructor: The default settings are first read, and
 * then potentially overwritten by values contained in a file.  */
RbSettings::RbSettings(void)
{

    readUserSettings();
}


const path& RbSettings::getModuleDir( void ) const
{
    return moduleDir;
}


size_t RbSettings::getLineWidth( void ) const
{
    // return the internal value
    return lineWidth;
}

size_t RbSettings::getScalingDensity( void ) const
{
    // return the internal value
    return scalingDensity;
}

bool RbSettings::getUseScaling( void ) const
{
    // return the internal value
    return useScaling;
}

int RbSettings::getDebugMCMC( void ) const
{
    // return the internal value
    return debugMCMC;
}

int RbSettings::getLogMCMC( void ) const
{
    // return the internal value
    return logMCMC;
}

std::string bool_to_string(bool b)
{
    return b ? "TRUE" : "FALSE";
}

std::string RbSettings::getOption(const std::string &key) const
{
    if ( key == "moduledir" )
    {
        return moduleDir.string();
    }
    else if ( key == "outputPrecision" )
    {
        return StringUtilities::to_string(outputPrecision);
    }
    else if ( key == "printNodeIndex" )
    {
        return printNodeIndex ? "true" : "false";
    }
    else if ( key == "tolerance" )
    {
        return StringUtilities::to_string(tolerance);
    }
    else if ( key == "linewidth" )
    {
        return StringUtilities::to_string(lineWidth);
    }
    else if ( key == "scalingDensity" )
    {
        return StringUtilities::to_string(scalingDensity);
    }
    else if ( key == "useScaling" )
    {
        return useScaling ? "true" : "false";
    }
    else if ( key == "debugMCMC" )
    {
        return std::to_string(debugMCMC);
    }
    else if ( key == "logMCMC" )
    {
        return std::to_string(logMCMC);
    }
    else if ( key == "echo" )
    {
        return bool_to_string(echo);
    }
    else if ( key == "continueOnError" )
    {
        return bool_to_string(continue_on_error);
    }
    else
    {
        std::cout << "Unknown user setting with key '" << key << "'." << std::endl;
    }
    
    return "";
}


size_t RbSettings::getOutputPrecision( void ) const
{
    // return the internal value
    return outputPrecision;
}


bool RbSettings::getPrintNodeIndex( void ) const
{
    // return the internal value
    return printNodeIndex;
}


double RbSettings::getTolerance( void ) const
{
    
    return tolerance;
}

bool RbSettings::getEcho( void ) const
{
    return echo;
}

bool RbSettings::getContinueOnError( void ) const
{
    return continue_on_error;
}


/** Initialize the user settings */
void RbSettings::readUserSettings(void)
{
    path user_dir = RevBayesCore::expandUserDir("~");
    
    // read the ini file, override defaults if applicable
    path settings_file_name = user_dir / ".RevBayes.ini";

    //    bool failed = false; //unused
    if ( is_regular_file( settings_file_name) )
    {
        std::ifstream readStream( settings_file_name.string() );
        std::string readLine = "";
        while ( safeGetline(readStream,readLine) )
        {
            std::vector<std::string> tokens = std::vector<std::string>();
            StringUtilities::stringSplit(readLine, "=", tokens);
            if (tokens.size() > 1)
            {
                setOption(tokens[0], tokens[1], false);
            }
        }
        
        readStream.close();
    }

//    writeUserSettings();
}


void RbSettings::listOptions() const
{
    std::cout << "moduledir = " << moduleDir << std::endl;
    std::cout << "outputPrecision = " << outputPrecision << std::endl;
    std::cout << "printNodeIndex = " << (printNodeIndex ? "true" : "false") << std::endl;
    std::cout << "tolerance = " << tolerance << std::endl;
    std::cout << "linewidth = " << lineWidth << std::endl;
    std::cout << "useScaling = " << (useScaling ? "true" : "false") << std::endl;
    std::cout << "scalingDensity = " << scalingDensity << std::endl;
    std::cout << "debugMCMC = " << debugMCMC << std::endl;
    std::cout << "logMCMC = " << logMCMC << std::endl;
    std::cout << "echo = " << echo << std::endl;
    std::cout << "continueOnError = " << continue_on_error << std::endl;
}


void RbSettings::setModuleDir(const path &md)
{
    if ( not is_directory(md) )
    {
        throw RbException()<<"Cannot set the help directory to "<<md<<".";
    }
    
    moduleDir = md;
    
    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::setLineWidth(size_t w)
{
    // replace the internal value with this new value
    lineWidth = w;
    
    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setUseScaling(bool w)
{
    // replace the internal value with this new value
    useScaling = w;

    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setScalingDensity(size_t w)
{
    // replace the internal value with this new value
    scalingDensity = w;
    
    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::setDebugMCMC(int d)
{
    // replace the internal value with this new value
    debugMCMC = d;
    
    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::setLogMCMC(int d)
{
    // replace the internal value with this new value
    logMCMC = d;

    // save the current settings for the future.
    writeUserSettings();
}

std::optional<bool> string_to_bool(const std::string& option)
{
    if (option == "TRUE" or option == "true")
        return true;
    else if (option == "FALSE" or option == "false")
        return false;
    else
        return {};
}


void RbSettings::setOption(const std::string &key, const std::string &v, bool write)
{

    std::string value = v;
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);

    if ( key == "moduledir" )
    {
        // Read from stream to handle quotes.
        std::istringstream input(value);
        input >> moduleDir;
    }
    else if ( key == "outputPrecision" )
    {
        outputPrecision = boost::lexical_cast<int>(value);
    }
    else if ( key == "printNodeIndex" )
    {
        printNodeIndex = value == "true";
    }
    else if ( key == "tolerance" )
    {
        tolerance = boost::lexical_cast<double>(value);
    }
    else if ( key == "linewidth" )
    {
        lineWidth = boost::lexical_cast<int>(value);
    }
    else if ( key == "useScaling" )
    {
        useScaling = value == "true";
    }
    else if ( key == "scalingDensity" )
    {
        size_t w = boost::lexical_cast<int>(value);
        if (w < 1)
            throw(RbException("scalingDensity must be an integer greater than 0"));
        
        scalingDensity = boost::lexical_cast<int>(value);
    }
    else if ( key == "debugMCMC" )
    {
        debugMCMC = boost::lexical_cast<int>(value);
    }
    else if ( key == "logMCMC" )
    {
        logMCMC = boost::lexical_cast<int>(value);
    }
    else if ( key == "echo" )
    {
        auto b = string_to_bool(value);
        if (b)
            echo = *b;
        else
            throw RbException()<<"setOption: expected a boolean, but got '"<<value<<"'";
    }
    else if ( key == "continueOnError" )
    {
        auto b = string_to_bool(value);
        if (b)
            continue_on_error = *b;
        else
            throw RbException()<<"setOption: expected a boolean, but got '"<<value<<"'";
    }
    else
    {
        std::cout << "Unknown user setting with key '" << key << "'." << std::endl;
    }

    // Maybe we should just write this particular setting instead of all settings?
    if ( write == true )
    {
        writeUserSettings();
    }
    
}


void RbSettings::setOutputPrecision(size_t p)
{
    // replace the internal value with this new value
    outputPrecision = p;

    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::setPrintNodeIndex(bool tf)
{
    // replace the internal value with this new value
    printNodeIndex = tf;
}


void RbSettings::setTolerance(double t)
{
    // replace the internal value with this new value
    tolerance = t;
    
    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setEcho(bool b)
{
    echo = b;

    // This is a per-session setting and so should not persist across sessions.
}

void RbSettings::setContinueOnError(bool b)
{
    continue_on_error = b;

    // This is a per-session setting and so should not persist across sessions.
}


void RbSettings::writeUserSettings( void )
{
    // Does this always work on windows?
    path user_dir = expandUserDir("~");
    
    // open the ini file
    path settings_file_name = user_dir / ".RevBayes.ini";

    std::ofstream writeStream( settings_file_name.string() );
    assert( moduleDir == "modules" or is_directory(moduleDir) );
    writeStream << "moduledir=" << moduleDir << std::endl;
    writeStream << "outputPrecision=" << outputPrecision << std::endl;
    writeStream << "printNodeIndex=" << (printNodeIndex ? "true" : "false") << std::endl;
    writeStream << "tolerance=" << tolerance << std::endl;
    writeStream << "linewidth=" << lineWidth << std::endl;
    writeStream << "useScaling=" << (useScaling ? "true" : "false") << std::endl;
    writeStream << "scalingDensity=" << scalingDensity << std::endl;

    // "echo" and "continueOnError" are per-session settings and so should not persist across sessions.

    writeStream.close();
}

void showDebug(const std::string& s, int level)
{
    if (RbSettings::userSettings().getLogMCMC() >= level)
	std::cerr<<s<<"\n";
}

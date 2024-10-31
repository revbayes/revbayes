#include "RbSettings.h"

#include <unistd.h>
#include <cstdlib> //includes std::atof
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

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
        outputPrecision = atoi(value.c_str());
    }
    else if ( key == "printNodeIndex" )
    {
        printNodeIndex = value == "true";
    }
    else if ( key == "tolerance" )
    {
        //std::string::size_type sz;     // alias of size_t
        //tolerance = std::stod (value,&sz);
        tolerance = (double)atof(value.c_str());
    }
    else if ( key == "linewidth" )
    {
        //std::string::size_type sz;     // alias of size_t
        //lineWidth = std::stoi (value,&sz);
        lineWidth = atoi(value.c_str());
    }
    else if ( key == "useScaling" )
    {
        useScaling = value == "true";
    }
    else if ( key == "scalingDensity" )
    {
        size_t w = atoi(value.c_str());
        if (w < 1)
            throw(RbException("scalingDensity must be an integer greater than 0"));
        
        scalingDensity = atoi(value.c_str());
    }
    else
    {
        std::cout << "Unknown user setting with key '" << key << "'." << std::endl;
    }
    
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
    writeStream.close();

}

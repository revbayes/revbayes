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
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#	ifdef _WIN32
#include <windows.h>
#   endif

namespace fs = boost::filesystem;




/** Default constructor: The default settings are first read, and
 * then potentially overwritten by values contained in a file.  */
RbSettings::RbSettings(void)
{

	initializeUserSettings();
}


const std::string& RbSettings::getModuleDir( void ) const
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

bool RbSettings::getCollapseSampledAncestors( void ) const
{
    // return the internal value
    return collapseSampledAncestors;
}

std::string RbSettings::getOption(const std::string &key) const
{
    if ( key == "moduledir" )
    {
        return moduleDir;
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
    else if ( key == "collapseSampledAncestors" )
    {
        return collapseSampledAncestors ? "true" : "false";
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


std::string RbSettings::getWorkingDirectory( void ) const
{
    return fs::current_path().make_preferred().string();
}


/** Initialize the user settings */
#define	MAX_DIR_PATH	2048
void RbSettings::initializeUserSettings(void)
{
    moduleDir = "modules";      // the default module directory
    useScaling = true;         // the default useScaling
    scalingDensity = 1;         // the default scaling density
    lineWidth = 160;            // the default line width
    tolerance = 10E-10;         // set default value for tolerance comparing doubles
    outputPrecision = 7;
    printNodeIndex = true;      // print node indices of tree nodes as comments
    collapseSampledAncestors = true;
    
    std::string user_dir = RevBayesCore::RbFileManager::expandUserDir("~");
    
    // read the ini file, override defaults if applicable
    std::string settings_file_name = ".RevBayes.ini";
    RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(user_dir, settings_file_name);

 	
    //    bool failed = false; //unused
    if ( fm.isFile() )
    {
        std::ifstream readStream;
        fm.openFile( readStream );
        std::string readLine = "";
        while ( fm.safeGetline(readStream,readLine) )
        {
            std::vector<std::string> tokens = std::vector<std::string>();
            StringUtilities::stringSplit(readLine, "=", tokens);
            if (tokens.size() > 1)
            {
                setOption(tokens[0], tokens[1], false);
            }
        }
        
        fm.closeFile(readStream);
    }

    // initialize the current directory to be the directory the binary is sitting in
#	ifdef _WIN32
    
    char buffer[MAX_DIR_PATH];
    GetModuleFileName( NULL, buffer, MAX_DIR_PATH );
    std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
    fs::current_path( std::string( buffer ).substr( 0, pos) );

#	else

    char cwd[MAX_DIR_PATH+1];
    if ( getcwd(cwd, MAX_DIR_PATH+1) )
    {
        fs::current_path( cwd );
    }
    
#   endif

    
    // save the current settings for the future.
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
    std::cout << "collapseSampledAncestors = " << (collapseSampledAncestors ? "true" : "false") << std::endl;
}


void RbSettings::setModuleDir(const std::string &md)
{
    
    RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(md);
    
    if ( !fm.isDirectory() )
    {
        throw RbException("Cannot set the help directory to '" + md + "'.");
    }
    
    moduleDir = fm.getFullFilePath();
    
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


void RbSettings::setCollapseSampledAncestors(bool w)
{
    // replace the internal value with this new value
    collapseSampledAncestors = w;

    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::setOption(const std::string &key, const std::string &v, bool write)
{

    std::string value = v;
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);

    if ( key == "moduledir" )
    {
        moduleDir = value;
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
    else if ( key == "collapseSampledAncestors" )
    {
        collapseSampledAncestors = value == "true";
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


void RbSettings::setWorkingDirectory(const std::string &wd)
{
    try
    {
        RevBayesCore::RbFileManager fm;
        std::string wd2 = fm.expandUserDir(wd);
        fs::current_path( wd2 );
    }
    catch (...)
    {
        throw RbException()<<"Cannot set the current directory to \""<<wd<<"\"";
    }
    
    // save the current settings for the future.
    writeUserSettings();
}


void RbSettings::writeUserSettings( void )
{
    std::string user_dir = RevBayesCore::RbFileManager::expandUserDir("~");
    
    // open the ini file
    std::string settings_file_name = ".RevBayes.ini";
    RevBayesCore::RbFileManager fm = RevBayesCore::RbFileManager(user_dir, settings_file_name);

    std::ofstream writeStream;
    fm.openFile( writeStream );
    writeStream << "moduledir=" << moduleDir << std::endl;
    writeStream << "outputPrecision=" << outputPrecision << std::endl;
    writeStream << "printNodeIndex=" << (printNodeIndex ? "true" : "false") << std::endl;
    writeStream << "tolerance=" << tolerance << std::endl;
    writeStream << "linewidth=" << lineWidth << std::endl;
    writeStream << "useScaling=" << (useScaling ? "true" : "false") << std::endl;
    writeStream << "scalingDensity=" << scalingDensity << std::endl;
    writeStream << "collapseSampledAncestors=" << (collapseSampledAncestors ? "true" : "false") << std::endl;
    fm.closeFile( writeStream );

}

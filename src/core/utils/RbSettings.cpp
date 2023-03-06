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

	initializeUserSettings();
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

bool RbSettings::getCollapseSampledAncestors( void ) const
{
    // return the internal value
    return collapseSampledAncestors;
}


bool RbSettings::getUseBeagle( void ) const
{
#if defined( RB_BEAGLE )
    // return the internal value
    return useBeagle;
#else
    return false;
#endif
}

#if defined( RB_BEAGLE )
const std::string& RbSettings::getBeagleDevice( void ) const
{
    // return the internal value
    return beagleDevice;
}

size_t RbSettings::getBeagleResource( void ) const
{
    // return the internal value
    return beagleResource;
}

bool RbSettings::getBeagleUseDoublePrecision( void ) const
{
    // return the internal value
    return beagleUseDoublePrecision;
}

size_t RbSettings::getBeagleMaxCPUThreads( void ) const
{
    // return the internal value
    return beagleMaxCPUThreads;
}

const std::string& RbSettings::getBeagleScalingMode( void ) const
{
    // return the internal value
    return beagleScalingMode;
}

size_t RbSettings::getBeagleDynamicScalingFrequency( void ) const
{
    // return the internal value
    return beagleDynamicScalingFrequency;
}
#endif /* RB_BEAGLE */


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
    else if ( key == "collapseSampledAncestors" )
    {
        return collapseSampledAncestors ? "true" : "false";
    }
#if defined( RB_BEAGLE )
    else if ( key == "useBeagle" )
    {
        return useBeagle ? "true" : "false";
    }
    else if ( key == "beagleDevice" )
    {
        return beagleDevice;
    }
    else if ( key == "beagleResource" )
    {
        return StringUtilities::to_string(beagleResource);
    }
    else if ( key == "beagleUseDoublePrecision" )
    {
        return beagleUseDoublePrecision ? "true" : "false";
    }
    else if ( key == "beagleMaxCPUThreads" )
    {
        return StringUtilities::to_string(beagleMaxCPUThreads);
    }
    else if ( key == "beagleScalingMode" )
    {
        return beagleScalingMode;
    }
    else if ( key == "beagleDynamicScalingFrequency" )
    {
        return StringUtilities::to_string(beagleDynamicScalingFrequency);
    }
#endif /* RB_BEAGLE */
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
#define	MAX_DIR_PATH	2048
void RbSettings::initializeUserSettings(void)
{
    moduleDir = "modules";      // the default module directory
    useScaling = true;          // the default useScaling
    scalingDensity = 1;         // the default scaling density
    lineWidth = 160;            // the default line width
    tolerance = 10E-10;         // set default value for tolerance comparing doubles
    outputPrecision = 7;
    printNodeIndex = true;      // print node indices of tree nodes as comments
    collapseSampledAncestors = true;

#if defined( RB_BEAGLE )
    useBeagle                     = false;      // don't use BEAGLE by default
    beagleDevice                  = "auto";     // auto select BEAGLE device by default
    beagleResource                = 0;          // the default BEAGLE resource
    beagleUseDoublePrecision      = true;       // BEAGLE will use double precision by default
    beagleMaxCPUThreads           = -1;          // no max set, auto threading up to number of cores
    //beagleScalingMode            = "dynamic";   // dynamic rescale as needed plus fixed frequency
    beagleScalingMode             = "manual";   // manually rescale as needed
    beagleDynamicScalingFrequency = 100;        // dynamic rescale every 100 evaluations by default
#endif /* RB_BEAGLE */
    
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
    std::cout << "collapseSampledAncestors = " << (collapseSampledAncestors ? "true" : "false") << std::endl;

#if defined( RB_BEAGLE )
    std::cout << "useBeagle = " << (useBeagle ? "true" : "false") << std::endl;
    std::cout << "beagleDevice = " << beagleDevice << std::endl;
    std::cout << "beagleResource = " << beagleResource << std::endl;
    std::cout << "beagleUseDoublePrecision = " << (beagleUseDoublePrecision ? "true" : "false") << std::endl;
    std::cout << "beagleMaxCPUThreads = " << beagleMaxCPUThreads << std::endl;
    std::cout << "beagleScalingMode = " << beagleScalingMode << std::endl;
    std::cout << "beagleDynamicScalingFrequency = " << beagleDynamicScalingFrequency << std::endl;
#endif /* RB_BEAGLE */
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


void RbSettings::setCollapseSampledAncestors(bool w)
{
    // replace the internal value with this new value
    collapseSampledAncestors = w;

    // save the current settings for the future.
    writeUserSettings();
}

#if defined( RB_BEAGLE )
void RbSettings::setUseBeagle(bool w)
{
    // replace the internal value with this new value
    useBeagle = w;

    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleDevice(const std::string &d)
{
    // replace the internal value with this new value
    beagleDevice = d;

    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleResource(size_t w)
{
    // replace the internal value with this new value
    beagleResource = w;

    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleUseDoublePrecision(bool s)
{
    // replace the internal value with this new value
    beagleUseDoublePrecision = s;
    
    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleMaxCPUThreads(size_t w)
{
    // replace the internal value with this new value
    beagleMaxCPUThreads = w;
    
    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleScalingMode(const std::string &bsm)
{
    // replace the internal value with this new value
    beagleScalingMode = bsm;
    
    // save the current settings for the future.
    //writeUserSettings();
}

void RbSettings::setBeagleDynamicScalingFrequency(size_t w)
{
    // replace the internal value with this new value
    beagleDynamicScalingFrequency = w;
    
    // save the current settings for the future.
    //writeUserSettings();
}
#endif /* RB_BEAGLE */


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
    else if ( key == "collapseSampledAncestors" )
    {
        collapseSampledAncestors = value == "true";
    }
#if defined( RB_BEAGLE )
    else if ( key == "useBeagle" )
    {
        useBeagle = value == "true";
    }
    else if ( key == "beagleDevice" )
    {
	if ( value == "cpu" || value == "cpu_sse" || value == "cpu_avx" || value == "gpu_opencl" || value == "gpu_cuda" || value == "auto" )
        {
            beagleDevice = value;
        }
        else
        {
            throw(RbException("beagleDevice must be set to cpu, cpu_sse, cpu_avx, gpu_opencl, gpu_cuda, or auto"));
        }
    }
    else if ( key == "beagleResource" )
    {
        size_t w = atoi(value.c_str());
        if (w < 0)
            throw(RbException("beagleResource must be a positive integer"));
        
        beagleResource = atoi(value.c_str());
    }
    else if ( key == "beagleUseDoublePrecision" )
    {
        beagleUseDoublePrecision = value == "true";
    }
    else if ( key == "beagleMaxCPUThreads" )
    {
        size_t w = atoi(value.c_str());
        if (w < 0)
            throw(RbException("beagleMaxCPUThreads must be a positive integer"));
        
        beagleMaxCPUThreads = atoi(value.c_str());
    }
    else if ( key == "beagleScalingMode" )
    {
	if ( value == "dynamic" || value == "auto" || value == "always" || value == "manual")
        {
            beagleScalingMode = value;
        }
        else
        {
            throw(RbException("beagleScalingMode must be set to dynamic, auto, always, or manual"));
        }
    }
    else if ( key == "beagleDynamicScalingFrequency" )
    {
        size_t w = atoi(value.c_str());
        if (w < 0)
            throw(RbException("beagleDynamicScalingFrequency must be a positive integer"));
        
        beagleDynamicScalingFrequency = atoi(value.c_str());
    }
#endif /* RB_BEAGLE */
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
    writeStream << "collapseSampledAncestors=" << (collapseSampledAncestors ? "true" : "false") << std::endl;
    writeStream.close();

}

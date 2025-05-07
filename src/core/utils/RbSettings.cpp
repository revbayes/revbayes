#include "RbSettings.h"

#include <unistd.h>
#include <cstdlib> //includes std::atof
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <boost/lexical_cast.hpp>

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
    return scaling_density;
}

const std::string& RbSettings::getScalingMethod( void ) const
{
    // return the internal value
    return scaling_method;
}

bool RbSettings::getScalingPerMixture( void ) const
{
    // return the internal value
    return scaling_per_mixture;
}

bool RbSettings::getUseScaling( void ) const
{
    // return the internal value
    return use_scaling;
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
    else if ( key == "partialLikelihoodStoring" )
    {
        return partial_likelihood_storing;
    }
    else if ( key == "scalingDensity" )
    {
        return StringUtilities::to_string(scaling_density);
    }
    else if ( key == "scalingMethod" )
    {
        return scaling_method;
    }
    else if ( key == "scalingPerMixture" )
    {
        return scaling_per_mixture ? "true" : "false";
    }
    else if ( key == "useScaling" )
    {
        return use_scaling ? "true" : "false";
    }
    else if ( key == "debugMCMC" )
    {
        return std::to_string(debugMCMC);
    }
    else if ( key == "logMCMC" )
    {
        return std::to_string(logMCMC);
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


const std::string& RbSettings::getPartialLikelihoodStoring( void ) const
{
    // return the internal value
    return partial_likelihood_storing;
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
    std::cout << "partialLikelihoodStoring = " << partial_likelihood_storing << std::endl;
    std::cout << "useScaling = " << (use_scaling ? "true" : "false") << std::endl;
    std::cout << "scalingMethod = " << scaling_method << std::endl;
    std::cout << "scalingDensity = " << scaling_density << std::endl;
    std::cout << "scalingPerMixture = " << (scaling_per_mixture ? "true" : "false") << std::endl;
    std::cout << "debugMCMC = " << debugMCMC << std::endl;
    std::cout << "logMCMC = " << logMCMC << std::endl;

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

void RbSettings::setPartialLikelihoodStoring(const std::string s)
{
    // replace the internal value with this new value
    partial_likelihood_storing = s;
    
    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setUseScaling(bool w)
{
    // replace the internal value with this new value
    use_scaling = w;

    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setScalingMethod(const std::string s)
{
    // replace the internal value with this new value
    scaling_method = s;
    
    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setScalingDensity(size_t w)
{
    // replace the internal value with this new value
    scaling_density = w;
    
    // save the current settings for the future.
    writeUserSettings();
}

void RbSettings::setScalingPerMixture(bool tf)
{
    // replace the internal value with this new value
    scaling_per_mixture = tf;
    
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
    else if ( key == "partialLikelihoodStoring" )
    {
        partial_likelihood_storing = value;
    }
    else if ( key == "useScaling" )
    {
        use_scaling = (value == "true");
    }
    else if ( key == "scalingMethod" )
    {
        scaling_method = value;
    }
    else if ( key == "scalingDensity" )
    {
        size_t w = boost::lexical_cast<int>(value);
        if (w < 1)
            throw(RbException("scalingDensity must be an integer greater than 0"));
        
        scaling_density = boost::lexical_cast<int>(value);
    }
    else if ( key == "scalingPerMixture" )
    {
        scaling_per_mixture = (value == "true");
    }
    else if ( key == "debugMCMC" )
    {
        debugMCMC = boost::lexical_cast<int>(value);
    }
    else if ( key == "logMCMC" )
    {
        logMCMC = boost::lexical_cast<int>(value);
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
//    assert( moduleDir == "modules" or is_directory(moduleDir) );
    writeStream << "moduledir=" << moduleDir << std::endl;
    writeStream << "outputPrecision=" << outputPrecision << std::endl;
    writeStream << "printNodeIndex=" << (printNodeIndex ? "true" : "false") << std::endl;
    writeStream << "tolerance=" << tolerance << std::endl;
    writeStream << "linewidth=" << lineWidth << std::endl;
    writeStream << "partialLikelihoodStoring=" << partial_likelihood_storing << std::endl;
    writeStream << "useScaling=" << (use_scaling ? "true" : "false") << std::endl;
    writeStream << "scalingDensity=" << scaling_density << std::endl;
    writeStream << "scalingMethod=" << scaling_method << std::endl;
    writeStream << "scalingPerMixture=" << (scaling_per_mixture ? "true" : "false") << std::endl;
    writeStream.close();

}

void showDebug(const std::string& s, int level)
{
    if (RbSettings::userSettings().getLogMCMC() >= level)
	std::cerr<<s<<"\n";
}

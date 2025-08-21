#ifndef RbSettings_H
#define RbSettings_H

#include <cstddef>
#include <iosfwd>
#include <string> // IWYU pragma: keep
#include <sstream>

#include "RbFileManager.h"

class RbSettings {

public:
    static RbSettings&          userSettings(void)                                  //!< Get a reference to the singleton RbSettings object
    {
        static RbSettings settings;
        return settings;
    }

    void                        listOptions(void) const;                            //!< Retrieve a list of all user options and their current values
    void                        readUserSettings(void);                             //!< Read settings from the config file.
    void                        writeUserSettings(void);                            //!< Write the current settings into a file.
    
    // Access functions
    size_t                      getLineWidth(void) const;                           //!< Retrieve the line width that will be used for the screen width when printing
    const RevBayesCore::path&   getModuleDir(void) const;                           //!< Retrieve the module directory name
    std::string                 getOption(const std::string &k) const;              //!< Retrieve a user option
    size_t                      getOutputPrecision(void) const;                     //!< Retrieve the default output precision width
    bool                        getPrintNodeIndex(void) const;                      //!< Retrieve the flag whether we should print node indices
    size_t                      getScalingDensity(void) const;                      //!< Retrieve the scaling density that determines how often to scale the likelihood in CTMC models
    double                      getTolerance(void) const;                           //!< Retrieve the tolerance for comparing doubles
    bool                        getUseScaling(void) const;                          //!< Retrieve the flag whether we should scale the likelihood in CTMC models
    int                         getDebugMCMC(void) const;                           //!< How much work should we perform to check MCMC?
    int                         getLogMCMC(void) const;                             //!< How much logging should we perform to check MCMC?
    bool                        getEcho() const;                                    //!< Should we print user commands to the screen (in a script or -e expr).
    bool                        getContinueOnError() const;                         //!< Should we continue after an error (in a script or -e expr).

    // setters
    void                        setLineWidth(size_t w);                             //!< Set the line width that will be used for the screen width when printing
    void                        setModuleDir(const RevBayesCore::path &md);         //!< Set the module directory name
    void                        setOutputPrecision(size_t p);                       //!< Set the default output precision width
    void                        setOption(const std::string &k, const std::string &v, bool write);  //!< Set the key value pair.
    void                        setPrintNodeIndex(bool tf);                         //!< Set the flag whether we should print node indices
    void                        setScalingDensity(size_t w);                        //!< Set the scaling density n, where CTMC likelihoods are scaled every n-th node (min 1)
    void                        setTolerance(double t);                             //!< Set the tolerance for comparing double
    void                        setUseScaling(bool s);                              //!< Set the flag whether we should scale the likelihood in CTMC models
    void                        setDebugMCMC(int d);                                //!< How much work should we perform to check MCMC?
    void                        setLogMCMC(int d);                                  //!< How much logging should we perform to check MCMC?
    void                        setEcho(bool b);                                    //!< Should we print user commands to the screen (in a script or -e expr)
    void                        setContinueOnError(bool b);                         //!< Should we continue after an error (int a script or -e expr)
    
private:
    RbSettings(void);                                   //!< Default constructor
    RbSettings(const RbSettings&) = delete;             //!< Prevent copy
    ~RbSettings(void) {};                               //!< Delete function table
    RbSettings&                 operator=(const RbSettings& s) = delete;                     //!< Prevent assignment

    // Variables that have user settings
    size_t                      lineWidth = 160;
    RevBayesCore::path          moduleDir = "modules";
    size_t                      outputPrecision = 7;
    bool                        printNodeIndex = true;                                     //!< Should the node index of a tree be printed as a comment?
    size_t                      scalingDensity = 1;
    double                      tolerance=10e-10;                                          //!< Tolerance for comparison of doubles
    bool                        useScaling=true;
    int                         debugMCMC = 0;
    int                         logMCMC = 0;

    bool                        echo=false;
    bool                        continue_on_error=false;
};

void showDebug(const std::string& s, int level=1);

class withReason
{
    double value;
    int level = 2;
    std::ostringstream reason;

public:
    template <typename T>
    withReason& operator<<(const T& t) {reason<<t; return *this;}

    operator double() const {showDebug(reason.str(), level); return value;}

    withReason(double d):value(d) {}
    withReason(double d, int l):value(d), level(l) {}
};

#endif


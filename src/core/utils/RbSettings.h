#ifndef RbSettings_H
#define RbSettings_H

#include <cstddef>
#include <iosfwd>
#include <string> // IWYU pragma: keep

#include "RbFileManager.h"

class RbSettings {

    public:
        static RbSettings&          userSettings(void)                                  //!< Get a reference to the singleton RbSettings object
                                    {
                                       static RbSettings settings = RbSettings();
									   return settings;
                                    }
   
        void                        initializeUserSettings(void);                       //!< Initialize the user settings to default values
    
    
        // Access functions
        bool                        getCollapseSampledAncestors(void) const;            //!< Retrieve the whether to should display sampled ancestors as 2-degree nodes when printing
        size_t                      getLineWidth(void) const;                           //!< Retrieve the line width that will be used for the screen width when printing
        const RevBayesCore::path&   getModuleDir(void) const;                           //!< Retrieve the module directory name
        std::string                 getOption(const std::string &k) const;              //!< Retrieve a user option
        size_t                      getOutputPrecision(void) const;                     //!< Retrieve the default output precision width
        const std::string&          getPartialLikelihoodStoring(void) const;            //!< Retrieve the method for partial likelihood storing
        bool                        getPrintNodeIndex(void) const;                      //!< Retrieve the flag whether we should print node indices
        size_t                      getScalingDensity(void) const;                      //!< Retrieve the scaling density that determines how often to scale the likelihood in CTMC models
        const std::string&          getScalingMethod(void) const;                       //!< Retrieve the scaling method for the likelihood in CTMC models
        bool                        getScalingPerMixture(void) const;                   //!< Retrieve the flag whether we should scale the likelihood in CTMC models per mixture category or all categories jointly
        double                      getTolerance(void) const;                           //!< Retrieve the tolerance for comparing doubles
        bool                        getUseScaling(void) const;                          //!< Retrieve the flag whether we should scale the likelihood in CTMC models
        void                        listOptions(void) const;                            //!< Retrieve a list of all user options and their current values

    bool                            getUseBeagle(void) const;                           //!< Retrieve the flag whether we should use the BEAGLE library in CTMC models
#if defined( RB_BEAGLE )
    const std::string&              getBeagleDevice(void) const;                        //!< Retrieve the BEAGLE device that is being used
    size_t                          getBeagleResource(void) const;                      //!< Retrieve the BEAGLE resource to be used
    bool                            getBeagleUseDoublePrecision(void) const;            //!< Retrieve the flag whether BEAGLE will use double precision floating point format
    size_t                          getBeagleMaxCPUThreads(void) const;                 //!< Retrieve the maximum number of CPU threads BEAGLE is set to use        
    const std::string&              getBeagleScalingMode(void) const;                   //!< Retrieve the BEAGLE numerical scaling mode
    size_t                          getBeagleDynamicScalingFrequency(void) const;       //!< Retrieve the BEAGLE evaluation frequency for calculation of updated numerical scaling factors
#endif /* RB_BEAGLE */


    // setters
    void                            setCollapseSampledAncestors(bool);                  //!< Set whether to should display sampled ancestors as 2-degree nodes when printing
    void                            setLineWidth(size_t w);                             //!< Set the line width that will be used for the screen width when printing
    void                            setModuleDir(const RevBayesCore::path &md);         //!< Set the module directory name
    void                            setOutputPrecision(size_t p);                       //!< Set the default output precision width
    void                            setOption(const std::string &k, const std::string &v, bool write);  //!< Set the key value pair.
    void                            setPartialLikelihoodStoring(const std::string s);   //!< Set the method for partial likelihood storing (node|branch|both)
    void                            setPrintNodeIndex(bool tf);                         //!< Set the flag whether we should print node indices
    void                            setScalingDensity(size_t w);                        //!< Set the scaling density n, where CTMC likelihoods are scaled every n-th node (min 1)
    void                            setScalingMethod(const std::string s);              //!< Set the scaling method, either node or threshold
    void                            setScalingPerMixture(bool tf);                      //!< Set the flag whether we should scale the likelihood in CTMC models per mixture category or all categories together
    void                            setTolerance(double t);                             //!< Set the tolerance for comparing double
    void                            setUseScaling(bool tf);                             //!< Set the flag whether we should scale the likelihood in CTMC models
    void                            setWorkingDirectory(const std::string &wd);         //!< Set the current working directory
       
#if defined( RB_BEAGLE )
    void                            setUseBeagle(bool s);                               //!< Set the flag whether we should use the BEAGLE library in CTMC models
    void                            setBeagleDevice(const std::string &bsm);            //!< Set the BEAGLE device to use
    void                            setBeagleResource(size_t w);                        //!< Set the BEAGLE resource to be used
    void                            setBeagleUseDoublePrecision(bool s);                //!< Set the flag whether BEAGLE will use double precision floating point format
    void                            setBeagleMaxCPUThreads(size_t w);                   //!< Set the maximum number of CPU threads BEAGLE is set to use        
    void                            setBeagleScalingMode(const std::string &bsm);       //!< Set the BEAGLE numerical scaling mode
    void                            setBeagleDynamicScalingFrequency(size_t w);         //!< Set the BEAGLE evaluation frequency for calculation of updated numerical scaling factors
#endif /* RB_BEAGLE */


    private:
                                    RbSettings(void);                                   //!< Default constructor
                                    RbSettings(const RbSettings&) {}                    //!< Prevent copy
                                   ~RbSettings(void) {}                                 //!< Delete function table
        RbSettings&                 operator=(const RbSettings& s);                     //!< Prevent assignment


        void                        writeUserSettings(void);                            //!< Write the current settings into a file.
    
    // Variables that have user settings
        bool                        collapseSampledAncestors;
        size_t                      lineWidth;
        RevBayesCore::path          moduleDir;
        size_t                      outputPrecision;
        bool                        printNodeIndex;                                     //!< Should the node index of a tree be printed as a comment?
        double                      tolerance;                                          //!< Tolerance for comparison of doubles
        std::string                 partial_likelihood_storing;                         //!< How should we store partial likelihoods (branch|node|both)?
    // related to rescaling to avoid underflow
        bool                        use_scaling;
        size_t                      scaling_density;
        bool                        scaling_per_mixture;
        std::string                 scaling_method;

#if defined( RB_BEAGLE )
    bool                        useBeagle;
    std::string                 beagleDevice;
    size_t                      beagleResource;        
    bool                        beagleUseDoublePrecision;
    size_t                      beagleMaxCPUThreads;
    std::string                 beagleScalingMode;
    size_t                      beagleDynamicScalingFrequency;
#endif /* RB_BEAGLE */
    
};

#endif


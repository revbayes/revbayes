#ifndef RbSettings_H
#define RbSettings_H

#include <stddef.h>
#include <iosfwd>
#include <string> // IWYU pragma: keep

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
        const std::string&          getModuleDir(void) const;                           //!< Retrieve the module directory name
        std::string                 getOption(const std::string &k) const;              //!< Retrieve a user option
        size_t                      getOutputPrecision(void) const;                     //!< Retrieve the default output precision width
        bool                        getPrintNodeIndex(void) const;                      //!< Retrieve the flag whether we should print node indices
        size_t                      getScalingDensity(void) const;                      //!< Retrieve the scaling density that determines how often to scale the likelihood in CTMC models
        double                      getTolerance(void) const;                           //!< Retrieve the tolerance for comparing doubles
        bool                        getUseScaling(void) const;                          //!< Retrieve the flag whether we should scale the likelihood in CTMC models
        const std::string&          getWorkingDirectory(void) const;                    //!< Retrieve the current working directory
        void                        listOptions(void) const;                            //!< Retrieve a list of all user options and their current values
    
#if defined( RB_BEAGLE )
    bool                        getUseBeagle(void) const;                           //!< Retrieve the flag whether we should use the BEAGLE library in CTMC models
    bool                        getBeagleAuto(void) const;                          //!< Retrieve the flag whether we should automatically select the BEAGLE resource
    size_t                      getBeagleResource(void) const;                      //!< Retrieve the BEAGLE resource to be used
    bool                        getBeagleUseDoublePrecision(void) const;            //!< Retrieve the flag whether BEAGLE will use double precision floating point format
    size_t                      getBeagleMaxCPUThreads(void) const;                 //!< Retrieve the maximum number of CPU threads BEAGLE is set to use        
    const std::string&          getBeagleScalingMode(void) const;                   //!< Retrieve the BEAGLE numerical scaling mode
    size_t                      getBeagleDynamicScalingFrequency(void) const;       //!< Retrieve the BEAGLE evaluation frequency for calculation of updated numerical scaling factors
#endif /* RB_BEAGLE */


        // setters
        void                        setCollapseSampledAncestors(bool);                  //!< Set whether to should display sampled ancestors as 2-degree nodes when printing
        void                        setLineWidth(size_t w);                             //!< Set the line width that will be used for the screen width when printing
        void                        setModuleDir(const std::string &md);                //!< Set the module directory name
        void                        setOutputPrecision(size_t p);                       //!< Set the default output precision width
        void                        setOption(const std::string &k, const std::string &v, bool write);  //!< Set the key value pair.
        void                        setPrintNodeIndex(bool tf);                         //!< Set the flag whether we should print node indices
        void                        setScalingDensity(size_t w);                        //!< Set the scaling density n, where CTMC likelihoods are scaled every n-th node (min 1)
        void                        setTolerance(double t);                             //!< Set the tolerance for comparing double
        void                        setUseScaling(bool s);                              //!< Set the flag whether we should scale the likelihood in CTMC models
        void                        setWorkingDirectory(const std::string &wd);         //!< Set the current working directory
       
#if defined( RB_BEAGLE )
    void                        setUseBeagle(bool s);                               //!< Set the flag whether we should use the BEAGLE library in CTMC models
    void                        setBeagleAuto(bool s);                              //!< Set the flag whether we should automatically select the BEAGLE resource
    void                        setBeagleResource(size_t w);                        //!< Set the BEAGLE resource to be used
    void                        setBeagleUseDoublePrecision(bool s);                //!< Set the flag whether BEAGLE will use double precision floating point format
    void                        setBeagleMaxCPUThreads(size_t w);                   //!< Set the maximum number of CPU threads BEAGLE is set to use        
    void                        setBeagleScalingMode(const std::string &bsm);       //!< Set the BEAGLE numerical scaling mode
    void                        setBeagleDynamicScalingFrequency(size_t w);         //!< Set the BEAGLE evaluation frequency for calculation of updated numerical scaling factors
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
        std::string                 moduleDir;
        size_t                      outputPrecision;
        bool                        printNodeIndex;                                     //!< Should the node index of a tree be printed as a comment?
        size_t                      scalingDensity;
        double                      tolerance;                                          //!< Tolerance for comparison of doubles
        bool                        useScaling;
        std::string                 workingDirectory;

#if defined( RB_BEAGLE )
    bool                        useBeagle;
    bool                        beagleAuto;
    size_t                      beagleResource;        
    bool                        beagleUseDoublePrecision;
    size_t                      beagleMaxCPUThreads;
    std::string                 beagleScalingMode;
    size_t                      beagleDynamicScalingFrequency;
#endif /* RB_BEAGLE */
    
};

#endif


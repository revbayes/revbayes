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
        size_t                      getLineWidth(void) const;                           //!< Retrieve the line width that will be used for the screen width when printing
        const RevBayesCore::path&   getModuleDir(void) const;                           //!< Retrieve the module directory name
        std::string                 getOption(const std::string &k) const;              //!< Retrieve a user option
        size_t                      getOutputPrecision(void) const;                     //!< Retrieve the default output precision width
        bool                        getPrintNodeIndex(void) const;                      //!< Retrieve the flag whether we should print node indices
        size_t                      getScalingDensity(void) const;                      //!< Retrieve the scaling density that determines how often to scale the likelihood in CTMC models
        double                      getTolerance(void) const;                           //!< Retrieve the tolerance for comparing doubles
        bool                        getUseScaling(void) const;                          //!< Retrieve the flag whether we should scale the likelihood in CTMC models
        void                        listOptions(void) const;                            //!< Retrieve a list of all user options and their current values

        // setters
        void                        setLineWidth(size_t w);                             //!< Set the line width that will be used for the screen width when printing
        void                        setModuleDir(const RevBayesCore::path &md);         //!< Set the module directory name
        void                        setOutputPrecision(size_t p);                       //!< Set the default output precision width
        void                        setOption(const std::string &k, const std::string &v, bool write);  //!< Set the key value pair.
        void                        setPrintNodeIndex(bool tf);                         //!< Set the flag whether we should print node indices
        void                        setScalingDensity(size_t w);                        //!< Set the scaling density n, where CTMC likelihoods are scaled every n-th node (min 1)
        void                        setTolerance(double t);                             //!< Set the tolerance for comparing double
        void                        setUseScaling(bool s);                              //!< Set the flag whether we should scale the likelihood in CTMC models
    
    private:
                                    RbSettings(void);                                   //!< Default constructor
                                    RbSettings(const RbSettings&) {}                    //!< Prevent copy
                                   ~RbSettings(void) {}                                 //!< Delete function table
        RbSettings&                 operator=(const RbSettings& s);                     //!< Prevent assignment


        void                        writeUserSettings(void);                            //!< Write the current settings into a file.
    
		// Variables that have user settings
        size_t                      lineWidth;
        RevBayesCore::path          moduleDir;
        size_t                      outputPrecision;
        bool                        printNodeIndex;                                     //!< Should the node index of a tree be printed as a comment?
        size_t                      scalingDensity;
        double                      tolerance;                                          //!< Tolerance for comparison of doubles
        bool                        useScaling;
};

#endif


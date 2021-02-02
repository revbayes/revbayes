#ifndef RlMonitor_H
#define RlMonitor_H

#include "Monitor.h"
#include "RevVariable.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <cstring>

namespace RevLanguage {

    /** @brief Base class for all monitors in Rev Language
     *
     * Monitors are tasked with saving information about one or several variable DAG node(s).
     * */
    
    class Monitor : public WorkspaceToCoreWrapperObject<RevBayesCore::Monitor> {
        
    public:
        
        Monitor(void);  //!< Default constructor
        Monitor(RevBayesCore::Monitor *m);  //!< Copy constructor with internal value
        
        // Basic utility functions
        virtual Monitor*                            clone(void) const = 0;
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        virtual void                                printValue(std::ostream& o, bool user) const;
        virtual void                                printValue(std::ostream& o) const;


        // Member method functions
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);
        
        std::vector<std::string>                    getConstructorFunctionAliases(void) const;  //!< Get the alternative constructor functions in Rev.
        std::string                                 getConstructorFunctionName(void) const;  //!< Get the full constructor function in Rev
        virtual std::vector<std::string>            getMonitorAliases(void) const { return std::vector<std::string>(); }  //!< Get the alternative names used for the constructor function in Rev.
        virtual std::string                         getMonitorName(void) const = 0; //!< Get name of the monitor for constructor function in Rev
        
    protected:

        //! compare variable names according to alphabetical order
        static bool                                 compareVarNames(RevPtr<const RevVariable> i, RevPtr<const RevVariable> j) { return strcmp(i->getName().c_str(), j->getName().c_str()) < 0; }

    };
    
}

#endif

#ifndef RlFileMonitor_H
#define RlFileMonitor_H

#include "RlMonitor.h"

namespace RevLanguage {

    /** @brief Base class for all file monitors in Rev Language
     *
     * File onitors are tasked with saving information about one or several variable DAG node(s) to a file.
     * */
    
    class FileMonitor : public Monitor {
        
    public:
        
        FileMonitor(void);  //!< Default constructor
        
        // Basic utility functions
        virtual FileMonitor*                        clone(void) const = 0;
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);

        virtual const MemberRules&                  getParameterRules(void) const;                                        //!< Get member rules (const)

    protected:

        virtual void                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                   append;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   printgen;
        RevPtr<const RevVariable>                   separator;
        RevPtr<const RevVariable>                   version;

    };
    
}

#endif

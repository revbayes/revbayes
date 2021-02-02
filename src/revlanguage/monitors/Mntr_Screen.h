#ifndef Mntr_Screen_H
#define Mntr_Screen_H

#include "RlMonitor.h"

namespace RevLanguage {

    /**
     * @copybrief RevBayesCore::ScreenMonitor
     * @see RevBayesCore::ScreenMonitor for the internal object
     **/
    
    class Mntr_Screen : public Monitor {
        
    public:
        
        Mntr_Screen(void);
        
        // Basic utility functions
        virtual Mntr_Screen*                        clone(void) const;
        void                                        constructInternalObject(void);  //!< Build a new internal ScreenMonitor.
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMonitorName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);
        
        std::vector<RevPtr<const RevVariable> >     vars;  //!< variables to monitor
        RevPtr<const RevVariable>                   printgen;  //!< print frequency
        RevPtr<const RevVariable>                   prior; //!< whether to print the prior
        RevPtr<const RevVariable>                   posterior;  //!< whether to print the posterior
        RevPtr<const RevVariable>                   likelihood;  //!< whether to print the likelihood
        
    };
    
}

#endif

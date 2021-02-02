#ifndef SRC_REVLANGUAGE_MONITORS_MNTRNEXUSFILE_H_
#define SRC_REVLANGUAGE_MONITORS_MNTRNEXUSFILE_H_

#include <iosfwd>
#include <vector>

#include "RlMonitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;

/** @copybrief RevBayesCore::NexusMonitor
 * @see Mntr_ExtendedNewickFile for output in table format
**/
class Mntr_NexusFile: public Monitor {

public:

    Mntr_NexusFile();
    virtual Mntr_NexusFile*                     clone(void) const;
    void                                        constructInternalObject(void);  //!< Build a new internal NexusMonitor.

    static const std::string&                   getClassType(void);
    static const TypeSpec&                      getClassTypeSpec(void);
    virtual const TypeSpec&                     getTypeSpec(void) const;

    std::string                                 getMonitorName(void) const;
    const MemberRules&                          getParameterRules(void) const;

    virtual void                                printValue(std::ostream& o) const;

protected:

    void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

    std::vector<RevPtr<const RevVariable> >     vars; //!< data associated with the tree
    RevPtr<const RevVariable>                   filename; //!< output file
    RevPtr<const RevVariable>                   tree; //!< monitored tree
    RevPtr<const RevVariable>                   isNodeParameter; //!< whether data is on the nodes or branches
    RevPtr<const RevVariable>                   writeTaxa; //!< whether to write a taxa block
    RevPtr<const RevVariable>                   printgen; //!< print frequency
};

} /* namespace RevLanguage */

#endif /* SRC_REVLANGUAGE_MONITORS_MNTRNEXUSFILE_H_ */

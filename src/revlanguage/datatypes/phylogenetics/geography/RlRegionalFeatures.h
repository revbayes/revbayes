#ifndef RlRegionalFeatures_H
#define RlRegionalFeatures_H

#include "ModelObject.h"
#include "RegionalFeatures.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class RlRegionalFeatures : public ModelObject<RevBayesCore::RegionalFeatures> {
        
    public:
        RlRegionalFeatures(void);                                                                          //!< Default constructor
        RlRegionalFeatures(RevBayesCore::RegionalFeatures *m);                                                    //!< Default constructor
                                            RlRegionalFeatures(RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatures> *d);                        //!< Default constructor
//        RlRegionalFeatures(const WorkspaceVector<WorkspaceVector<RlString> >& m);
        
        
        // Basic utility functions
        RlRegionalFeatures*                            clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions

        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
        
    private:
        RevBayesCore::RegionalFeatures*            rf;
    };
    
}
#endif /* defined(RlRegionalFeatures_H) */

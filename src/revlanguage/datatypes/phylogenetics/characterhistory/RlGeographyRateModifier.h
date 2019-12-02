//#ifndef RlGeographyRateModifier_H
//#define RlGeographyRateModifier_H
//
//#include "ModelObject.h"
//#include "GeographyRateModifier.h"
//#include "TypedDagNode.h"
//
//#include <ostream>
//#include <string>
//
//namespace RevLanguage {
//    
//    class GeographyRateModifier : public ModelObject<RevBayesCore::GeographyRateModifier> {
//        
//    public:
//        
//        GeographyRateModifier(void);                                                                                  //!< Default constructor
//        GeographyRateModifier(RevBayesCore::GeographyRateModifier *m);                                       //!< Default constructor
//        GeographyRateModifier(RevBayesCore::TypedDagNode<RevBayesCore::GeographyRateModifier> *d);           //!< Default constructor
//        
//        // Basic utility functions
//        GeographyRateModifier*              clone(void) const;                                                          //!< Clone object
//        static const std::string&           getClassType(void);                                                         //!< Get Rev type
//        static const TypeSpec&              getClassTypeSpec(void);                                                     //!< Get class type spec
//        const TypeSpec&                     getTypeSpec(void) const;                                                    //!< Get language type of the object
//        
//        std::string                         getGuiName(void) { return ""; }
//        std::string                         getGuiUnicodeSymbol(void) { return ""; }
//        std::string                         getGuiInfo(void) { return ""; }
//    };
//}
//#endif

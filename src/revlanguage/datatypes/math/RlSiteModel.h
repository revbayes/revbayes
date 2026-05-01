/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the SiteModel class.
 *
 */

#ifndef RlSiteModel_H
#define RlSiteModel_H

#include <ostream>
#include <vector>

#include "RlRateGenerator.h"
#include "SiteModel.h"
#include "RevPtr.h"

namespace RevBayesCore { class SiteModel; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;

    class SiteModel : public ModelObject<RevBayesCore::SiteModel> {

    public:

        SiteModel(void);                                                                                             //!< Default constructor
        SiteModel(const RevBayesCore::SiteModel& m);                                                                 //!< Default constructor
        SiteModel(RevBayesCore::SiteModel *m);                                                                       //!< Default constructor
        SiteModel(RevBayesCore::TypedDagNode<RevBayesCore::SiteModel> *m);                                           //!< Default constructor

        // Basic utility functions
        SiteModel*                          clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object

        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions

    protected:
        void                                initMethods(void);
    };

}

#endif

/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the SiteMixtureModel class.
 *
 * @brief Declaration of RlSiteMixtureModel
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-06 20:14:22 +0200 (Mon, 06 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: Real.h 1746 2012-08-06 18:14:22Z hoehna $
 */

#ifndef RlSiteMixtureModel_H
#define RlSiteMixtureModel_H

#include <ostream>
#include <vector>

#include "RlRateGenerator.h"
#include "SiteMixtureModel.h"
#include "RevPtr.h"

namespace RevBayesCore { class SiteMixtureModel; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;

    class SiteMixtureModel : public ModelObject<RevBayesCore::SiteMixtureModel> {

    public:

        SiteMixtureModel(void);                                                                                             //!< Default constructor
        SiteMixtureModel(const RevBayesCore::SiteMixtureModel& m);                                                  //!< Default constructor
        SiteMixtureModel(RevBayesCore::SiteMixtureModel *m);                                                        //!< Default constructor
        SiteMixtureModel(RevBayesCore::TypedDagNode<RevBayesCore::SiteMixtureModel> *m);                            //!< Default constructor

        // Basic utility functions
        SiteMixtureModel*           clone(void) const;                                                                      //!< Clone object
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

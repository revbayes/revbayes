/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the MixtureModel class.
 *
 * @brief Declaration of RlMixtureModel
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

#ifndef RlMixtureModel_H
#define RlMixtureModel_H

#include <ostream>
#include <vector>

#include "RlRateGenerator.h"
#include "MixtureModel.h"
#include "RevPtr.h"

namespace RevBayesCore { class MixtureModel; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
    class MixtureModel : public ModelObject<RevBayesCore::MixtureModel> {
        
    public:
        
        MixtureModel(void);                                                                                                         //!< Default constructor
        MixtureModel(const RevBayesCore::MixtureModel& m);                                                                          //!< Default constructor
        MixtureModel(RevBayesCore::MixtureModel *m);                                                                                //!< Default constructor
        MixtureModel(RevBayesCore::TypedDagNode<RevBayesCore::MixtureModel> *m);                                                    //!< Default constructor
        
        // Basic utility functions
        MixtureModel*                       clone(void) const;                                                                      //!< Clone object
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

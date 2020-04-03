#ifndef RlDiscretizedContinuousCharacterData_H
#define RlDiscretizedContinuousCharacterData_H

#include <iostream>
#include <string>
#include <vector>

#include "DiscretizedContinuousCharacterData.h"
#include "ModelObject.h"
#include "RlHomologousCharacterData.h"
#include "TypedDagNode.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlConstantNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
    class DiscretizedContinuousCharacterData : public ModelObject<RevBayesCore::DiscretizedContinuousCharacterData>, public HomologousCharacterData {
        
    public:
        DiscretizedContinuousCharacterData(void);                                                                                                                      //!< Default constructor
        DiscretizedContinuousCharacterData(const RevBayesCore::DiscretizedContinuousCharacterData &d);                                                                            //!< Constructor from new core data type pointer
        DiscretizedContinuousCharacterData(RevBayesCore::DiscretizedContinuousCharacterData *d);                                                                                  //!< Constructor from new core data type pointer
        DiscretizedContinuousCharacterData(RevBayesCore::TypedDagNode<RevBayesCore::DiscretizedContinuousCharacterData>*d);                                                       //!< Constructor with DAG node
        
        virtual ~DiscretizedContinuousCharacterData();
        
        
        // Basic utility functions
        DiscretizedContinuousCharacterData*                     clone(void) const;                                                                          //!< Clone object
        void                                                    concatenate(const RevObject& d, std::string type = "") const;                               //!< Concatenate two sequences
        void                                                    concatenate(const DiscretizedContinuousCharacterData& d, std::string type = "") const;                 //!< Concatenate two sequences
        static const std::string&                               getClassType(void);                                                                         //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                     //!< Get class type spec
        const TypeSpec&                                         getTypeSpec(void) const;                                                                    //!< Get language type of the object
        
        RevPtr<RevVariable>                                     executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found);     //!< Execute member method
        
        
    private:
        
        void                                                    initMethods(void);
        
    };
    
}

#endif 

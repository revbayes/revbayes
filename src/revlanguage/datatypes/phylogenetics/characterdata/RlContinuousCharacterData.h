#ifndef RlContinuousCharacterData_H
#define RlContinuousCharacterData_H

#include <iostream>
#include <string>
#include <vector>

#include "ContinuousCharacterData.h"
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
#include "RlConstantNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
    class ContinuousCharacterData : public ModelObject<RevBayesCore::ContinuousCharacterData>, public HomologousCharacterData {
        
    public:
        ContinuousCharacterData(void);                                                                                                                      //!< Default constructor
        ContinuousCharacterData(const RevBayesCore::ContinuousCharacterData &d);                                                                            //!< Constructor from new core data type pointer
        ContinuousCharacterData(RevBayesCore::ContinuousCharacterData *d);                                                                                  //!< Constructor from new core data type pointer
        ContinuousCharacterData(RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData>*d);                                                       //!< Constructor with DAG node
        
        virtual ~ContinuousCharacterData();
        
        
        // Basic utility functions
        ContinuousCharacterData*                                clone(void) const;                                                                          //!< Clone object
        void                                                    concatenate(const RevObject& d, std::string type = "") const;                               //!< Concatenate two sequences
        void                                                    concatenate(const ContinuousCharacterData& d, std::string type = "") const;                 //!< Concatenate two sequences
        static const std::string&                               getClassType(void);                                                                         //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                     //!< Get class type spec
        const TypeSpec&                                         getTypeSpec(void) const;                                                                    //!< Get language type of the object
        
        RevPtr<RevVariable>                                     executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found);     //!< Execute member method
        
        
    private:
        
        void                                                    initMethods(void);
        
    };
    
}

#endif 

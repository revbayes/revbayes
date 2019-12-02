#ifndef AddWorkspaceVectorType_H
#define AddWorkspaceVectorType_H

#include <string>
#include "Func_modelVector.h"
#include "Func__vectorIndexOperator.h"
#include "ModelVector.h"
#include "Workspace.h"

namespace RevLanguage {
    
    template <typename objType, int d>
    class AddWorkspaceVectorType {
        
        public:
            //!< Add the given object and a vector version of it.
            static bool                 addTypeToWorkspace( Workspace &w, objType *o )
                                            {
                                            w.addType( o );
                                            w.addFunction( new Func_modelVector<objType>() );
                                            w.addFunction( new Func__vectorIndexOperator<objType>() );
                                            AddWorkspaceVectorType<ModelVector<objType>, d-1>::addTypeToWorkspace( w, new ModelVector<objType>() );
                                            return true;
                                            }
    };
    
    template <typename objType>
    class AddWorkspaceVectorType<objType,0> {
        
        public:
            //!< Add the given object.
            static bool                 addTypeToWorkspace( Workspace &w, objType *o )
                                            {
                                            w.addType( o );
                                            return true;
                                            }
    };
}

#endif


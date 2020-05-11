#ifndef RLMatrixRealSymmetric_H
#define RLMatrixRealSymmetric_H

#include <iostream>
#include <vector>

#include "ModelObject.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"


namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the symmetrix matrix real.
     *
     * The RevLanguage wrapper of the symmetric matrix real simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the SymmetricMatrixReal.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot and Sebastian Hoehna)
     * @since 2014-03-27, version 1.0
     *
     */
    class MatrixRealSymmetric : public MatrixReal {
        
    public:
        
        MatrixRealSymmetric(void);
        MatrixRealSymmetric(const RevBayesCore::MatrixReal &m);
        MatrixRealSymmetric(RevBayesCore::MatrixReal* m);
        MatrixRealSymmetric(RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> *mat);                                                            //!< Construct from DAG node
        
        
        // the value type definition
        virtual MatrixRealSymmetric*    clone(void) const;                                                      //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        virtual const TypeSpec&         getTypeSpec(void) const;                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>     executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions

    protected:

    private:
        void                            initializeMethods(void);

    };

}

#endif /* defined(__revbayes__MatrixRealSymmetric__) */

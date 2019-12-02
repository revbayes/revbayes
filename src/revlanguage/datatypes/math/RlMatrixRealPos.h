#ifndef RlMatrixRealPos_H
#define RlMatrixRealPos_H

#include <iostream>



#include <iostream>
#include <vector>

#include "MatrixReal.h"
#include "ModelObject.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"


namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the matrix real.
     *
     * The RevLanguage wrapper of the matrix real simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the MatrixRealPos.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot and Sebastian Hoehna)
     * @since 2014-03-27, version 1.0
     *
     */
    class MatrixRealPos : public MatrixReal  {
        
    public:
        MatrixRealPos(void);                                                                       //!< Default constructor
        MatrixRealPos(const RevBayesCore::MatrixReal& m);
        MatrixRealPos(RevBayesCore::MatrixReal* m);
        MatrixRealPos(RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> *mat);                  //!< Construct from DAG node
        
        // the value type definition
        virtual MatrixRealPos*              clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        virtual const TypeSpec&             getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
        
    private:
        void                                initializeMethods(void);
    };
    
}

#endif /* defined(__revbayes__MatrixRealPos__) */

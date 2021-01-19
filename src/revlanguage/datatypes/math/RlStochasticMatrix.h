#ifndef RlStochasticMatrix_H
#define RlStochasticMatrix_H

#include <iostream>
#include <vector>

#include "MatrixReal.h"
#include "ModelObject.h"
#include "RlMatrixRealPos.h"
#include "TypedDagNode.h"


namespace RevLanguage {
    
    
    /**
     * The RevLanguage datatype for the (row) stochastic matrix. A (row)
     * stochastic matrix has rows that sum to 1.
     *
     */
    class StochasticMatrix : public MatrixRealPos  {
        
    public:
        StochasticMatrix(void);                                                                       //!< Default constructor
        StochasticMatrix(const RevBayesCore::MatrixReal& m);
        StochasticMatrix(RevBayesCore::MatrixReal* m);
        StochasticMatrix(RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> *mat);                  //!< Construct from DAG node
        
        // the value type definition
        virtual StochasticMatrix*           clone(void) const;                                                  //!< Clone object
        void                                constructInternalObject(void);                                      //!< We construct the a new internal MCMC object.
        static const std::string&           getClassType(void);                                                 //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                             //!< Get class type spec
        virtual const TypeSpec&             getTypeSpec(void) const;                                            //!< Get language type of the object
        std::string                         getConstructorFunctionName(void) const;                             //!< Get the name used for the constructor function in Rev.
        std::vector<std::string>            getConstructorFunctionAliases(void) const;                          //!< Get the aliases used for the constructor function in Rev.

        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
        const MemberRules&                  getParameterRules(void) const;                                      //!< Get member rules (const)

        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
        
    protected:
        void                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable

        RevPtr<const RevVariable>           arg;

    private:
        void                                initializeMethods(void);
    };
    
}

#endif /* defined(__revbayes__StochasticMatrix__) */

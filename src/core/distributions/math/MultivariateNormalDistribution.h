/**
 * @file
 * This file contains the declaration of the MultivariateNormally distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a normal distribution in the model graph.
 *
 * @brief Declaration of the stochastic DAG node base class.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date:$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id:$
 */



#ifndef MultivariateNormalDistribution_H
#define MultivariateNormalDistribution_H

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "MatrixReal.h"

#include <vector>

namespace RevBayesCore {
    
    class MultivariateNormalDistribution : public TypedDistribution< RbVector<double> > {
        
    public:
        MultivariateNormalDistribution(const TypedDagNode< RbVector<double> > *m, const TypedDagNode<MatrixReal>* cov, const TypedDagNode<MatrixReal>* prec, const TypedDagNode<double>* sc);
        virtual                                                         ~MultivariateNormalDistribution(void);                                                //!< Virtual destructor
        
        
        // public member functions
        MultivariateNormalDistribution*                                 clone(void) const;
        
        MatrixReal                                                      computeContrasts(void);
        RevLanguage::RevPtr<RevLanguage::RevVariable>                   executeProcedure(const std::string &n, const std::vector<DagNode*> args, bool &f);  //!< execute the procedure
        void                                                            clampAt(size_t i, double v);
        double                                                          computeLnProbability(void);
        void                                                            redrawValue(void);

    protected:
        // Parameter management functions
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*                         mean;
        const TypedDagNode< MatrixReal >*                               covariance;
        const TypedDagNode< MatrixReal >*                               precision;
        const TypedDagNode< double >*                                   scale;

        std::vector<bool>                                               observed;
        std::vector<double>                                             observations;
        
    };
    
}

#endif

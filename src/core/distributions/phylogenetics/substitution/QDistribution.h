#ifndef QDistribution_H
#define QDistribution_H

#include "RateMatrix_Rational.h"

#include "RbBoolean.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Distribution on a rate matrix.
     *
     * The ....
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class QDistribution : public TypedDistribution<RateGenerator>, public MemberObject< Boolean > {
//    class QDistribution : public TypedDistribution<RateGenerator> {

    public:
        QDistribution (const TypedDagNode<RbVector<double> >* alpha,
                       const TypedDagNode<double>* rho);                                                                        //!< Constructor
        
        // public member functions
        QDistribution*               clone(void) const;                                                                         //!< Create an independent clone

    protected:

        // Parameter management functions
        double                                          computeLnProbability(void);                                             //!< Compute the log-transformed probability of the current value.

        // Parameter management functions
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Boolean &rv) const;     //!< Map the member methods to internal function calls
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter

        void                                            keepSpecialization(DagNode *toucher);
        void                                            restoreSpecialization(DagNode *toucher);
        void                                            touchSpecialization(DagNode *toucher, bool touchAll);

    private:
        
        // helper functions
        void                                            redrawValue(void);

        bool                                            is_reversible;                  //!< Indicates whether the current value is reversible
        const TypedDagNode<RbVector<double> >*          alpha;                          //!< parameter for the prior on the stationary frequencies
        const TypedDagNode<double>*                     rho;                            //!< parameter for the probability of being time reversible

    };
}

#endif

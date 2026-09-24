#ifndef QDistribution_H
#define QDistribution_H

#include "RateMatrix_MPQ.h"

#include "MatrixReal.h"
#include "MemberObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Distribution on a rate matrix.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna and John Huelsenbeck)
     * @since 2014-03-18, version 1.0
     *
     */
    class QDistribution : public TypedDistribution<RateGenerator>, public MemberObject< Boolean >, public MemberObject< Simplex >, public MemberObject< RbVector<double> >, public MemberObject< double > {

    public:
        enum FIXED_MODEL { FREE, REVERSIBLE_ONLY, NON_REVERSIBLE_ONLY };

                                                QDistribution (const TypedDagNode<RbVector<double> >* alpha, double log_rho_rev, double log_rho_nr, FIXED_MODEL fm = FREE);   //!< Constructor
        
                                                // public member functions
        QDistribution*                          clone(void) const;                                                                         //!< Create an independent clone
        double                                  getLnRhoReversible(void) const { return log_rho_reversible; }
        double                                  getLnRhoNonReversible(void) const { return log_rho_non_reversible; }
        double                                  getLnPriorOdds(void) const { return log_rho_reversible - log_rho_non_reversible; }
        void                                    setLnPriorOdds(double delta);                                                              //!< Set log(rho_R/rho_N), keeping rho_R + rho_N = 1
        FIXED_MODEL                             getFixedModel(void) const { return fixed_model; }
        bool                                    isModelFixed(void) const { return fixed_model != FREE; }

    protected:

                                                // parameter management functions
        double                                  computeLnProbability(void);                                             //!< Compute the log-transformed probability of the current value.
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Boolean &rv) const;     //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Simplex &rv) const;     //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Map the member methods to internal function calls
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter

    private:
        
                                                // helper functions
        void                                    redrawValue(void);

        const TypedDagNode<RbVector<double> >*  alpha;                          //!< parameter for the prior on the stationary frequencies
        double                                  log_rho_reversible;               //!< log of the prior probability of the time-reversible model
        double                                  log_rho_non_reversible;            //!< log of the prior probability of the non-reversible model
        FIXED_MODEL                             fixed_model;                     //!< confine the chain to one model, for measuring them separately
    };
}

#endif

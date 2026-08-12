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
     * The ....
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class QDistribution : public TypedDistribution<RateGenerator>, public MemberObject< Boolean >, public MemberObject< Simplex >, public MemberObject< RbVector<double> >, public MemberObject< double > {

    public:
                                                QDistribution (const TypedDagNode<RbVector<double> >* alpha, double log_rho_rev, double log_rho_nr);                                                                        //!< Constructor
        
                                                // public member functions
        QDistribution*                          clone(void) const;                                                                         //!< Create an independent clone

        /* The prior on the model indicator.

           These are ordinary data members rather than DAG parameters, so that a
           move can retune them: with a large Bayes factor the chain would
           otherwise never visit the disfavoured model and its posterior
           probability could not be estimated at all. Tilting the prior until the
           chain splits its time evenly between the two models fixes that, and
           costs nothing, because

               BF_NR  =  (p_N / p_R) * (rho_R / rho_N)

           recovers the Bayes factor from any prior odds whatsoever. Only the
           variance of that estimate depends on the tilt, and it is smallest at
           50:50. Whatever tilt is used has to be reported, which is what the
           getters below are for. */
        double                                  getLnRhoReversible(void) const { return log_rho_reversible; }
        double                                  getLnRhoNonReversible(void) const { return log_rho_non_reversible; }
        double                                  getLnPriorOdds(void) const { return log_rho_reversible - log_rho_non_reversible; }
        void                                    setLnPriorOdds(double delta);                                                              //!< Set log(rho_R/rho_N), keeping rho_R + rho_N = 1

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
        double                                  log_rho_reversible;             //!< log of the prior probability of the time-reversible model
        double                                  log_rho_non_reversible;         //!< log of the prior probability of the non-reversible model
    };
}

#endif

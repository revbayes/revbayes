#ifndef AVMVNProposal_H
#define AVMVNProposal_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "Proposal.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class ContinuousStochasticNode;
class DagNode;
class RandomNumberGenerator;
class Simplex;
template <class valueType> class RbVector;
template <class valueType> class StochasticNode;

    /**
     * @brief Smart Multivariate Proposal for several parameters jointly.
     *
     *
     * This proposal learns the covariance matrix of a set of variables on
     * a transformed space and uses that to propose values intelligently.
     *
     * @author The RevBayes Development Core Team
     * @copyright GPL version 3
     * @since 2015-05-25, version 1.0
     *
     */
    class AVMVNProposal : public Proposal {

    public:
        AVMVNProposal( double s, double e, size_t n0, size_t c0, size_t m );                                 //!< Constructor

        void                                        addUntransformedScalar(StochasticNode<double> *v);                                    //!< Add an up-scaling variable
        void                                        addLogScalar(StochasticNode<double> *v);                                    //!< Add an up-scaling variable
        void                                        addLogitScalar(ContinuousStochasticNode *v);                                    //!< Add an up-scaling variable
        void                                        addUntransformedVector(StochasticNode<RbVector<double> > *v);                         //!< Add an up-scaling variable
        void                                        addLogConstrainedSumVector(RevBayesCore::StochasticNode<RevBayesCore::Simplex> *v);                         //!< Add an up-scaling variable
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        AVMVNProposal*                              clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        removeUntransformedScalar(StochasticNode<double> *v);                                    //!< Add an up-scaling variable
        void                                        removeLogScalar(StochasticNode<double> *v);                                    //!< Add an up-scaling variable
        void                                        removeLogitScalar(ContinuousStochasticNode *v);                                    //!< Add an up-scaling variable
        void                                        removeUntransformedVector(StochasticNode<RbVector<double> > *v);                         //!< Add an up-scaling variable
        void                                        removeLogConstrainedSumVector(RevBayesCore::StochasticNode<RevBayesCore::Simplex> *v);                         //!< Add an up-scaling variable
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal

    protected:

        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on


    private:
        // parameters

        std::vector<StochasticNode<double> *>               noTransformScalarVariables;
        std::vector<StochasticNode<double> *>               logTransformScalarVariables;
        std::vector<ContinuousStochasticNode *>             logitTransformScalarVariables;
        std::vector<StochasticNode<RbVector<double> > *>    noTransformVectorVariables;
        std::vector<StochasticNode<Simplex> *>              logConstrainedSumTransformVectorVariables;

        MatrixReal                                          C_emp;                                                   //!< The empirical covariance matrix of the samples
        MatrixReal                                          AVMVN_cholesky_L;                                       //!< Lower diagonal of Cholesky decomposition of

        size_t                                              waitBeforeLearning;                                      //!< How long to wait before tracking empirical covariances
        size_t                                              waitBeforeUsing;                                         //!< How long to wait before using the empirical covariances
        size_t                                              nTried;                                                  //!< How many times has this move been used?
        size_t                                              updates;                                                 //!< How many updates have been tried?
        size_t                                              maxUpdates;                                              //!< How many updates until we stop monitoring the covariances?

        double                                              sigma;                                                  //!< Variance of pre-learned (independent) normal proposal, also scales the MVN (proportional to variance, not SD)
        double                                              epsilon;                                                //!< Controls the weighting of the learned VCV and an Identity matrix ( eps * I + (1-eps) * empirical)
        size_t                                              dim;                                                    //!< Dimension of proposal
        double                                              lnHastingsratio;                                        //!< The Hastings ratio, so that helper functions can adjust as needed

        std::vector<double>                                 storedValues;                                           //!< The values before proposing the move, for Hastings ratio calculations
        std::vector<double>                                 storedValuesUntransformed;                              //!< The values before proposing the move, for resetting
        std::vector<double>                                 proposedValues;                                         //!< The values proposed by the move
        std::vector<double>                                 x_bar;                                                  //!< The averages in transformed space

        // functions
        void                                                calculateHastingsRatio(std::vector<double> x_prime, std::vector<double> x);
        void                                                getAVMVNMemberVariableValues(std::vector<double> *x);
        void                                                getAVMVNMemberVariableValuesUntransformed(std::vector<double> *x);
        void                                                resetAVMVNMemberVariableValues(std::vector<double> x);
        std::vector<double>                                 rMVNCholesky(std::vector<double> mu, MatrixReal L, RandomNumberGenerator& rng, double scale);


    };

}

#endif

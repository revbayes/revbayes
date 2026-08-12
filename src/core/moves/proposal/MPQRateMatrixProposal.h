#ifndef MPQRateMatrixProposal_H
#define MPQRateMatrixProposal_H

#include <cstddef>
#include <iosfwd>
#include <gmpxx.h>

#include "Polyhedron.h"
#include "Proposal.h"
#include "RateGenerator.h"
#include "RateMatrix_MPQ.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The time-reversible and non-reversible rate matrix proposal.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John & Sebastian)
     * @since 2009-09-08, version 1.0
     *
     */
    class MPQRateMatrixProposal : public Proposal {
        
    public:

        /* The sub-moves this proposal is a mixture of.

           They are enumerated because each one needs its own tuning parameter and
           its own acceptance rate. RevBayes hands tune() a single acceptance rate
           for the move as a whole, which is useless here: a rate of 0.4 could be
           six sub-moves at 0.4, or three at 0.8 and three at 0. So this class
           keeps its own per-sub-move counters instead. undoProposal is called on
           rejection and only on rejection, so rejections can be attributed to
           whichever sub-move was last proposed, and acceptance follows by
           subtraction from the number tried. */
        enum SUB_MOVE { PI_ALL,             //!< all four stationary frequencies, weights held fixed
                        PI_SINGLE,          //!< one stationary frequency
                        REV_BB_ALL,         //!< all six backbone weights, time-reversible model
                        REV_BB_SINGLE,      //!< one backbone weight, time-reversible model
                        NR_BB_ALL,          //!< all six backbone weights, non-reversible model
                        NR_BB_SINGLE,       //!< one backbone weight, non-reversible model
                        NR_U,               //!< the point (u1,u2,u3) inside the polyhedron
                        JUMP_TO_NR,         //!< reversible jump, to the non-reversible model
                        JUMP_TO_REV,        //!< reversible jump, to the time-reversible model
                        NUM_SUB_MOVES };


                                                            MPQRateMatrixProposal(StochasticNode<RateGenerator> *n);                           //!<  constructor
                                                            MPQRateMatrixProposal( const MPQRateMatrixProposal& p);                             //!<  copy constructor

                                                            // utility functions
        void                                                cleanProposal(void);                                                                //!< Clean up proposal
        MPQRateMatrixProposal*                              clone(void) const;                                                                  //!< Clone object
        double                                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                                setProposalTuningParameter(double tp);
        void                                                setTuneModelPrior(bool tf) { tune_model_prior = tf; }                                //!< Adapt the prior on the model indicator toward 50:50
        void                                                setTuneSubMoves(bool tf) { tune_sub_moves = tf; }                                    //!< Adapt each sub-move's tuning parameter against its own acceptance rate
        double                                              getSubMoveTuningParameter(SUB_MOVE m) const { return tuning[m]; }
        double                                              getSubMoveAcceptanceRate(SUB_MOVE m) const;                                          //!< Acceptance rate since the last tuning call, or -1 if untried
        static const char*                                  getSubMoveName(SUB_MOVE m);
        double                                              getLnPriorOdds(void) const { return ln_prior_odds; }                                 //!< log(rho_R / rho_N) currently in force
        double                                              getProportionReversible(void) const;                                                 //!< Fraction of visits to the reversible model since the last tuning call
        void                                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        double                                              applySubMove(SUB_MOVE m);                                                            //!< Carry out the chosen sub-move
        double                                              updateToNonReversible(void);
        double                                              updateToReversible(void);
        Polyhedron                                          poly;
                
                                                            // parameters
        StochasticNode<RateGenerator>*                      variable;
        
        /* One tuning parameter per sub-move.

           For the Dirichlet sub-moves this is the concentration, and note that it
           acts in the opposite direction to the usual scale parameter: a LARGER
           concentration means a proposal more tightly centred on the current
           value, hence smaller steps and a HIGHER acceptance rate. Tuning it
           therefore has to run the other way round from the textbook rule, which
           is written for a step size. For NR_U it is a window width and the usual
           direction applies. The two jump sub-moves have no tuning parameter; the
           polyhedron determines their proposal entirely. */
        double                                              tuning[NUM_SUB_MOVES];
        double                                              tune_target[NUM_SUB_MOVES];
        long                                                n_tried[NUM_SUB_MOVES];
        long                                                n_rejected[NUM_SUB_MOVES];
        SUB_MOVE                                            last_sub_move;
        bool                                                tune_sub_moves;

        /* Adaptation of the prior on the model indicator.

           When one model is strongly supported the chain stops visiting the other
           one and its posterior probability cannot be estimated. Tilting the prior
           until the chain divides its time evenly restores the estimate without
           biasing it: the Bayes factor recovered as (p_N/p_R)(rho_R/rho_N) does not
           depend on the tilt, only its variance does.

           Adaptation changes the target distribution, so it must not run while
           samples are being kept. Tune during burnin, then freeze, and report the
           value that was frozen. */
        bool                                                tune_model_prior;
        double                                              ln_prior_odds;                                                                      //!< log(rho_R / rho_N)
        long                                                n_reversible;                                                                       //!< visits to the reversible model since the last tuning call
        long                                                n_visits;                                                                           //!< visits to either model since the last tuning call
        long                                                n_tuning_calls;

        //!< The two indices of the last modified element.
        RateMatrix_MPQ                                      stored_Q;
        std::vector<mpq_class>                              W;
    };
    
}

#endif


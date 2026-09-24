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
                                                            MPQRateMatrixProposal(const MPQRateMatrixProposal& p);                             //!<  copy constructor

        // Basic utility functions
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
                                                            /* Concentration of the point drawn from the polyhedron.

                                                               The point (u1,u2,u3) is drawn from a triangulation of the polyhedron
                                                               about its center, which is the time reversible matrix at
                                                               (1/2,1/2,1/2). At alphaT = 1 the draw is uniform over the polyhedron,
                                                               which is an independence proposal from something close to the prior:
                                                               fine when the data say little about u, hopeless when they say a lot,
                                                               because a draw from the whole polyhedron then lands nowhere near where
                                                               the non-reversible model wants to be and the jump is refused.

                                                               Raising alphaT concentrates the draw near the center, which makes the
                                                               jump a small move: the proposed non-reversible matrix is close to the
                                                               current reversible one and their likelihoods are close, so the jump
                                                               turns on the prior and the Jacobian rather than on a large likelihood
                                                               difference. Lowering it below one pushes the draw out toward the
                                                               facets instead.

                                                               Which direction helps is a property of the data and cannot be
                                                               predicted, so measure it: monitor Q.getU() to see where the
                                                               non-reversible model actually sits, and read the jump acceptance rates
                                                               out of the operator summary. Nothing here affects correctness. The
                                                               proposal density is computed from the same triangulation that produced
                                                               the point, so any alphaT gives a valid move; only the acceptance rate
                                                               changes. */
        void                                                setPolyhedronAlpha(double x);
        double                                              getPolyhedronAlpha(void) const { return poly.getAlphaT(); }
                                                            /* Where to center the draw from the polyhedron.

                                                               Establish it with a short run confined to the non-reversible model,
                                                               monitoring Q.getU(), and pass the posterior mean of the three
                                                               coordinates. Anything the polyhedron cannot accommodate on a given
                                                               proposal is pulled back toward (1/2,1/2,1/2) for that proposal only;
                                                               see Polyhedron::chooseCenter. */
        void                                                setPolyhedronCenter(double u1, double u2, double u3);
        void                                                useReversiblePolyhedronCenter(void) { poly.useReversibleCenter(); }
                                                            /* Adapt alphaT during burnin against the jump acceptance rate.

                                                               There is no target acceptance rate to aim at here, as there is for an
                                                               ordinary random walk: this is an independence proposal and more
                                                               acceptance is simply better. So the tuning hill-climbs, keeping the
                                                               direction that improved matters and reversing when it did not, with a
                                                               step that shrinks as it goes. Like every other adaptation here it must
                                                               finish before sampling begins. */
        void                                                setTunePolyhedronAlpha(bool tf) { tune_poly_alpha = tf; }
        double                                              getSubMoveTuningParameter(SUB_MOVE m) const { return tuning[m]; }
        double                                              getSubMoveAcceptanceRate(SUB_MOVE m) const;                                          //!< Acceptance rate since the last tuning call, or -1 if untried
        double                                              getSubMoveAcceptanceRateTotal(SUB_MOVE m) const;                                     //!< Acceptance rate over the whole run, or -1 if untried
                                                            /* Whether the move just proposed was a jump between models, and which way.
                                                               MetropolisHastingsMove asks this, under MPQ_DEBUG_JUMPS, to print the
                                                               likelihood and prior ratios alongside the Hastings ratio for jumps only:
                                                               the proposal cannot see those ratios itself, and they are the pieces of
                                                               the acceptance probability that decide whether a jump ever succeeds. */
        bool                                                lastMoveWasJump(void) const { return last_sub_move == JUMP_TO_NR || last_sub_move == JUMP_TO_REV; }
        std::string                                         lastSubMoveName(void) const { return getSubMoveName( last_sub_move ); }
        double                                              lastSubMoveTuning(void) const { return tuning[last_sub_move]; }
                                                            /* Largest relative change any entry of the rate matrix underwent in the
                                                               last proposal, in double precision. This is what the likelihood sees.
                                                               A move that changes the matrix by a part in ten thousand and costs
                                                               tens of log-likelihood units is the signature of a likelihood that is
                                                               not a smooth function of the matrix; see the MOVEACC diagnostic. */
        double                                              lastMoveMaxRelativeChange(void) const;
        int                                                 lastJumpDirection(void) const { return last_sub_move == JUMP_TO_NR ? 1 : (last_sub_move == JUMP_TO_REV ? -1 : 0); }
        long                                                getSubMoveTriedTotal(SUB_MOVE m) const { return n_tried_total[m]; }
        static const char*                                  getSubMoveName(SUB_MOVE m);
        double                                              getLnPriorOdds(void) const { return ln_prior_odds; }                                 //!< log(rho_R / rho_N) currently in force
        double                                              getProportionReversible(void) const;
        long                                                getNumFailedToReversible(void) const { return num_failed_to_reversible; }
        long                                                getNumFailedToNonReversible(void) const { return num_failed_to_nonreversible; }                                                 //!< Fraction of visits to the reversible model since the last tuning call
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
                                                            /* The counters above are cleared by every tuning call, because that is the
                                                               window the tuning has to respond to. Reporting them is misleading: an
                                                               operator summary printed after the last tuning call describes a handful
                                                               of moves rather than the run. These two accumulate over the whole run
                                                               and are what should be looked at when judging whether a move works. */
        long                                                n_tried_total[NUM_SUB_MOVES];
        long                                                n_rejected_total[NUM_SUB_MOVES];
        long                                                num_dead_intervals[NUM_SUB_MOVES];                                                  //!< tuning intervals in which a sub-move was refused every time; see tune()
        long                                                tune_calls_total;
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

                                                            /* State for the search phase of the prior tuning; see tune().

                                                               A Robbins-Monro step is proportional to the observed log odds of the two
                                                               models, and that quantity saturates: once an interval contains no visits
                                                               at all to one of them, the observed log odds is bounded by the length of
                                                               the interval and reports the same thing whether the prior is off by ten
                                                               units or a thousand. The step size is then limited by ignorance rather
                                                               than by the schedule, and the tilt crawls. So while one model is
                                                               unvisited the tuning ignores the log odds entirely and simply searches:
                                                               it steps in the obvious direction, doubling the step each time it keeps
                                                               going the same way and halving it whenever it overshoots and has to turn
                                                               round. That brackets a tilt of any magnitude in a number of intervals
                                                               proportional to its logarithm rather than its square. Robbins-Monro takes
                                                               over, from the start of its schedule, as soon as both models are being
                                                               visited and the log odds means something again. */
        double                                              tune_step;
        int                                                 tune_last_dir;

                                                            /* The settings for the polyhedron are kept here as well as on the
                                                               Polyhedron itself.

                                                               Polyhedron cannot be copied: it owns a pool of vertices and its copy
                                                               constructor is deleted, so the copy constructor of this class leaves
                                                               poly default constructed rather than copying it. RevBayes clones a
                                                               proposal when it builds the MCMC, which means anything set on the
                                                               original never reaches the copy that actually runs. Holding the values
                                                               here and re-applying them in the copy constructor is what makes them
                                                               survive the clone.

                                                               This is easy to miss because it fails silently: the move still works,
                                                               it simply runs at the default concentration of one about the
                                                               time-reversible center, which is exactly the configuration one is
                                                               trying to change. */
        double                                              poly_alpha;
        double                                              poly_center[3];
        bool                                                poly_center_set;

                                                            // hill-climbing state for the alphaT tuning; see tune()
        bool                                                tune_poly_alpha;
        double                                              poly_alpha_factor;
        double                                              poly_alpha_last_rate;
        long                                                poly_alpha_calls;

                                                            /* Jumps refused because the move could not be defined, as opposed to jumps
                                                               rejected on the evidence. See reportProposalFailure. */
        long                                                num_failed_to_reversible;
        long                                                num_failed_to_nonreversible;
        static void                                         reportProposalFailure(const char* what, long count);
        double                                              currentLnPriorOdds(void) const;                                                     //!< the tilt actually in force, read from the distribution

#       ifdef MPQ_DEBUG_JUMPS
                                                            /* Temporary instrumentation. Compile with -DMPQ_DEBUG_JUMPS to have every
                                                               attempted reversible-jump move print the pieces of its Hastings ratio and
                                                               whether it was accepted. Nothing here is compiled in otherwise. */
        double                                              dbg_ln_jacobian;
        double                                              dbg_ln_poly;
        double                                              dbg_ln_hastings;
        double                                              dbg_u[3];                   //!< the point in the polyhedron: proposed on R->N, current on N->R
        int                                                 dbg_direction;              //!< +1 for R to N, -1 for N to R, 0 for no jump pending
        long                                                dbg_printed;
        void                                                dbgReport(const char* outcome);
#       endif

                                                            //!< The two indices of the last modified element.
        RateMatrix_MPQ                                      stored_Q;
        std::vector<mpq_class>                              W;
    };
}

#endif


#ifndef MPQRateMatrixProposal_H
#define MPQRateMatrixProposal_H

#include <cstddef>
#include <iosfwd>
#include "RbRational.h"

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

        /* The sub-moves for updating a MPQ rate matrix. They are enumerated because 
           each one needs its own tuning parameter and its own acceptance rate.  */
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

                                                            // basic utility functions
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
                                                               facets instead. */
        void                                                setPolyhedronAlpha(double x);
        double                                              getPolyhedronAlpha(void) const { return poly.getAlphaT(); }
                                                            /* Where to center the draw from the polyhedron.

                                                               Establish it with a short run confined to the non-reversible model,
                                                               monitoring Q.getU(), and pass the posterior mean of the three
                                                               coordinates. Anything the polyhedron cannot accommodate on a given
                                                               proposal is pulled back toward (1/2,1/2,1/2) for that proposal only. */
        void                                                setPolyhedronCenter(double u1, double u2, double u3);
        void                                                useReversiblePolyhedronCenter(void) { poly.useReversibleCenter(); }
        void                                                setTunePolyhedronAlpha(bool tf) { tune_poly_alpha = tf; }
        double                                              getSubMoveTuningParameter(SUB_MOVE m) const { return tuning[m]; }
        double                                              getSubMoveAcceptanceRate(SUB_MOVE m) const;                                          //!< Acceptance rate since the last tuning call, or -1 if untried
        double                                              getSubMoveAcceptanceRateTotal(SUB_MOVE m) const;                                     //!< Acceptance rate over the whole run, or -1 if untried
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
        
        double                                              tuning[NUM_SUB_MOVES];
        double                                              tune_target[NUM_SUB_MOVES];
        long                                                n_tried[NUM_SUB_MOVES];
        long                                                n_rejected[NUM_SUB_MOVES];
        long                                                n_tried_total[NUM_SUB_MOVES];
        long                                                n_rejected_total[NUM_SUB_MOVES];
        long                                                num_dead_intervals[NUM_SUB_MOVES];                                                  //!< tuning intervals in which a sub-move was refused every time; see tune()
        long                                                tune_calls_total;
        SUB_MOVE                                            last_sub_move;
        bool                                                tune_sub_moves;
        bool                                                tune_model_prior;
        double                                              ln_prior_odds;                                                                      //!< log(rho_R / rho_N)
        long                                                n_reversible;                                                                       //!< visits to the reversible model since the last tuning call
        long                                                n_visits;                                                                           //!< visits to either model since the last tuning call
        long                                                n_tuning_calls;
        double                                              tune_step;
        double                                              tune_bracket_lo;                                                                    //!< largest lnPriorOdds known to be too small (only the non-reversible model visited)
        double                                              tune_bracket_hi;                                                                    //!< smallest lnPriorOdds known to be too large (only the reversible model visited)
        double                                              poly_alpha;
        double                                              poly_center[3];
        bool                                                poly_center_set;

                                                            // hill-climbing state for the alphaT tuning; see tune()
        bool                                                tune_poly_alpha;
        double                                              poly_alpha_factor;
        double                                              poly_alpha_last_rate;
        long                                                poly_alpha_calls;

        long                                                num_failed_to_reversible;
        long                                                num_failed_to_nonreversible;
        static void                                         reportProposalFailure(const char* what, long count);
        double                                              currentLnPriorOdds(void) const;                                                     //!< the tilt actually in force, read from the distribution

                                                            //!< The two indices of the last modified element.
        RateMatrix_MPQ                                      stored_Q;
        std::vector<mpq_class>                              W;
    };
}

#endif


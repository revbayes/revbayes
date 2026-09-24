#include "MPQRateMatrixProposal.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbConstants.h"
#include "RateMatrix_MPQ.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "QDistribution.h"


#define MIN_FREQ    10e-4
#define MAX_LN_PRIOR_ODDS  1000.0
#define MIN_TUNE_STEP      0.001
#define MPQ_DEBUG_JUMP_LIMIT  200
#define MAX_TUNE_STEP      256.0
#define MIN_CONCENTRATION  1.0
#define MAX_CONCENTRATION  1.0e6
#define MIN_WINDOW         1.0e-3
#define MIN_TRIES_TO_CALL_DEAD  200
#define MAX_WINDOW         2.0
#define A           0
#define C           1
#define G           2
#define T           3

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MPQRateMatrixProposal::MPQRateMatrixProposal(StochasticNode<RateGenerator> *n) : Proposal(),
variable(n) {
    
    // starting values.
    tuning[PI_ALL]        = 100.0;
    tuning[PI_SINGLE]     =  50.0;
    tuning[REV_BB_ALL]    = 100.0;
    tuning[REV_BB_SINGLE] =  50.0;
    tuning[NR_BB_ALL]     = 500.0;
    tuning[NR_BB_SINGLE]  =  50.0;
    tuning[NR_U]          =   0.2;
    tuning[JUMP_TO_NR]    =   0.0;   // not tuned
    tuning[JUMP_TO_REV]   =   0.0;   // not tuned

    /* 0.44 is the usual target for a move on one coordinate; 0.234 is the usual
       target for a move on many at once. */
    tune_target[PI_ALL]        = 0.234;
    tune_target[PI_SINGLE]     = 0.44;
    tune_target[REV_BB_ALL]    = 0.234;
    tune_target[REV_BB_SINGLE] = 0.44;
    tune_target[NR_BB_ALL]     = 0.234;
    tune_target[NR_BB_SINGLE]  = 0.44;
    tune_target[NR_U]          = 0.234;
    tune_target[JUMP_TO_NR]    = 0.0;
    tune_target[JUMP_TO_REV]   = 0.0;

    for (int i=0; i<NUM_SUB_MOVES; i++)
        {
        n_tried[i]          = 0;
        n_rejected[i]       = 0;
        n_tried_total[i]    = 0;
        n_rejected_total[i] = 0;
        }
    last_sub_move  = PI_ALL;
    tune_sub_moves = true;

    tune_model_prior = false;
    ln_prior_odds    = 0.0;
    n_reversible     = 0;
    n_visits         = 0;
    n_tuning_calls   = 0;
    tune_step        = 1.0;
    tune_last_dir    = 0;
    tune_calls_total = 0;
    for (int i=0; i<NUM_SUB_MOVES; i++)
        num_dead_intervals[i] = 0;

    poly_alpha       = 1.0;
    poly_center[0]   = 0.5;
    poly_center[1]   = 0.5;
    poly_center[2]   = 0.5;
    poly_center_set  = false;

    num_failed_to_reversible    = 0;
    num_failed_to_nonreversible = 0;
#   ifdef MPQ_DEBUG_JUMPS
    dbg_ln_jacobian = 0.0;
    dbg_ln_poly     = 0.0;
    dbg_ln_hastings = 0.0;
    dbg_u[0]        = 0.0;
    dbg_u[1]        = 0.0;
    dbg_u[2]        = 0.0;
    dbg_direction   = 0;
    dbg_printed     = 0;
#   endif
    
    W.resize( 6 );

    // tell the base class to add the node
    addNode( variable );
}

MPQRateMatrixProposal::MPQRateMatrixProposal( const MPQRateMatrixProposal& p ) : Proposal( p ),
    variable( p.variable ),
    last_sub_move( p.last_sub_move ),
    tune_sub_moves( p.tune_sub_moves ),
    tune_model_prior( p.tune_model_prior ),
    ln_prior_odds( p.ln_prior_odds ),
    n_reversible( p.n_reversible ),
    n_visits( p.n_visits ),
    n_tuning_calls( p.n_tuning_calls ),
    tune_step( p.tune_step ),
    tune_last_dir( p.tune_last_dir ),
    tune_calls_total( p.tune_calls_total ),
    poly_alpha( p.poly_alpha ),
    poly_center_set( p.poly_center_set ),
    num_failed_to_reversible( p.num_failed_to_reversible ),
    num_failed_to_nonreversible( p.num_failed_to_nonreversible ),
    stored_Q( p.stored_Q ),
    W( p.W ) {

    /* poly is default constructed by this copy, because a Polyhedron cannot be
       copied. Re-apply the settings so that the clone RevBayes actually runs is
       configured the way the original was. */
    for (int i=0; i<3; i++)
        poly_center[i] = p.poly_center[i];
    poly.setAlphaT( poly_alpha );
    if ( poly_center_set == true )
        poly.setCenter( poly_center[0], poly_center[1], poly_center[2] );

    for (int i=0; i<NUM_SUB_MOVES; i++)
        {
        tuning[i]           = p.tuning[i];
        tune_target[i]      = p.tune_target[i];
        n_tried[i]          = p.n_tried[i];
        n_rejected[i]       = p.n_rejected[i];
        n_tried_total[i]    = p.n_tried_total[i];
        n_rejected_total[i] = p.n_rejected_total[i];
        num_dead_intervals[i] = p.num_dead_intervals[i];
        }
        
    // tell the base class to add the node
    addNode(variable);
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MPQRateMatrixProposal::cleanProposal(void) {

#   ifdef MPQ_DEBUG_JUMPS
    dbgReport("ACCEPTED");
#   endif
}

void MPQRateMatrixProposal::reportProposalFailure(const char* what, long count) {

    long p = 1;
    while (p < count)
        p *= 10;
    if (p != count)
        return;

    std::cerr << "MPQRateMatrixProposal: " << what;
    if (count == 1)
        std::cerr << " (first occurrence)";
    else
        std::cerr << " (occurrence " << count << ")";
    std::cerr << std::endl;
}

#ifdef MPQ_DEBUG_JUMPS
void MPQRateMatrixProposal::dbgReport(const char* outcome) {

    if ( dbg_direction == 0 || dbg_printed >= MPQ_DEBUG_JUMP_LIMIT )
        {
        dbg_direction = 0;
        return;
        }
    dbg_printed++;

    /* The point in the polyhedron is printed with the rest: on a jump to the
       non-reversible model it is the point that was proposed, and on a jump back it
       is the point the chain was sitting at. All three coordinates are one half
       exactly at time reversibility, so how far they sit from a half says how
       non-reversible the proposal was, which is what polyhedronCenter and
       polyhedronAlpha exist to control. */
    std::cerr << "JUMP " << (dbg_direction > 0 ? "R->N" : "N->R")
              << "  u = ("           << dbg_u[0] << ", " << dbg_u[1] << ", " << dbg_u[2] << ")"
              << "  lnJacobian = "   << dbg_ln_jacobian
              << "  lnPolyDensity = "<< dbg_ln_poly
              << "  lnHastings = "   << dbg_ln_hastings
              << "  lnPriorOdds = "  << currentLnPriorOdds()
              << "  -> "             << outcome << std::endl;

    if ( dbg_printed == MPQ_DEBUG_JUMP_LIMIT )
        std::cerr << "JUMP (further jump reports suppressed)" << std::endl;

    dbg_direction = 0;
}
#endif

double MPQRateMatrixProposal::currentLnPriorOdds(void) const {

    const QDistribution* q_dist = dynamic_cast<const QDistribution*>( &variable->getDistribution() );
    if (q_dist != NULL)
        return q_dist->getLnPriorOdds();
    return ln_prior_odds;
}

void MPQRateMatrixProposal::setPolyhedronAlpha(double x) {

    poly_alpha = x;
    poly.setAlphaT(x);
}

void MPQRateMatrixProposal::setPolyhedronCenter(double u1, double u2, double u3) {

    poly_center[0]  = u1;
    poly_center[1]  = u2;
    poly_center[2]  = u3;
    poly_center_set = true;
    poly.setCenter(u1, u2, u3);
}

double MPQRateMatrixProposal::lastMoveMaxRelativeChange(void) const {

    const RateMatrix_MPQ& Q = static_cast<const RateMatrix_MPQ&>(variable->getValue());
    double worst = 0.0;
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i == j)
                continue;
            double a = Q(i,j).get_d();
            double b = stored_Q(i,j).get_d();
            double d = fabs(a - b) / ( fabs(b) > 1e-300 ? fabs(b) : 1e-300 );
            if (d > worst)
                worst = d;
            }
        }
    return worst;
}

MPQRateMatrixProposal* MPQRateMatrixProposal::clone(void) const {
    
    return new MPQRateMatrixProposal(*this);
}

/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MPQRateMatrixProposal::getProposalName(void) const {

    static std::string name = "MPQRateMatrixProposal";
    return name;
}


double MPQRateMatrixProposal::getProposalTuningParameter(void) const {

    /* There is no single tuning parameter for this move; there are seven. This
       returns the stationary-frequency concentration so that the RevBayes
       machinery has something to report, but printParameterSummary is where the
       tuned values actually are. */
    return tuning[PI_ALL];
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MPQRateMatrixProposal::doProposal(void)  {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());

    n_visits++;
    if (Q.getIsReversible() == true)
        n_reversible++;

    bool jumps_allowed = true;
    const QDistribution* q_dist = dynamic_cast<const QDistribution*>( &variable->getDistribution() );
    if ( q_dist != NULL && q_dist->isModelFixed() == true )
        jumps_allowed = false;

    if (Q.getIsReversible() == true)
        {
        // the time-reversible model: stationary frequencies at fixed weights, or
        // the six backbone weights at fixed frequencies
        double u = rng->uniform01();
        double jump_threshold = ( jumps_allowed == true ) ? 0.90 : 1.0;
        if (u < 0.5 * jump_threshold)
            last_sub_move = ( rng->uniform01() < 0.25 ) ? PI_ALL : PI_SINGLE;
        else if (u < jump_threshold)
            last_sub_move = ( rng->uniform01() < 0.166 ) ? REV_BB_ALL : REV_BB_SINGLE;
        else
            last_sub_move = JUMP_TO_NR;
        }
    else
        {
        /* The non-reversible model. Its eight free weights are carried as a
           backbone plus (u1,u2,u3), so it takes two sub-moves to reach all of
           them; the stationary frequencies are a third, separate direction. */
        double u = rng->uniform01();
        double jump_threshold = ( jumps_allowed == true ) ? 0.90 : 1.0;
        if (u < (1.0/3.0) * jump_threshold)
            last_sub_move = ( rng->uniform01() < 0.25 ) ? PI_ALL : PI_SINGLE;
        else if (u < (2.0/3.0) * jump_threshold)
            last_sub_move = ( rng->uniform01() < 0.166 ) ? NR_BB_ALL : NR_BB_SINGLE;
        else if (u < jump_threshold)
            last_sub_move = NR_U;
        else
            last_sub_move = JUMP_TO_REV;
        }

    n_tried[last_sub_move]++;
    n_tried_total[last_sub_move]++;
    return applySubMove( last_sub_move );
}

double MPQRateMatrixProposal::applySubMove(SUB_MOVE m) {

    RandomNumberGenerator* rng = GLOBAL_RNG;
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());

    if (m != JUMP_TO_NR && m != JUMP_TO_REV)
        {
        // the two jump sub-moves store the matrix themselves, further down
        stored_Q = Q;
        }

    switch (m)
        {
        case PI_ALL:
            /* Holds all twelve weights fixed and moves only pi. Legal in both
               models: reversibility and stationarity are properties of the
               weights alone, so there is deliberately no reversibility guard. */
            return Q.updateStationaryFrequencies(rng, tuning[PI_ALL], 0.5);

        case PI_SINGLE:
            return Q.updateStationaryFrequenciesSingle(rng, tuning[PI_SINGLE], 0.5);

        case REV_BB_ALL:
            /* The backbone is a state coordinate of the time-reversible model and
               the prior is flat on it, so the Hastings ratio is the Dirichlet
               proposal ratio alone. */
            return Q.updateBackboneWeights(rng, tuning[REV_BB_ALL], 0.5);

        case REV_BB_SINGLE:
            return Q.updateBackboneWeightsSingle(rng, tuning[REV_BB_SINGLE], 0.5);

        case NR_BB_ALL:
            return Q.updateNonReversibleBackbone(rng, tuning[NR_BB_ALL], 0.5);

        case NR_BB_SINGLE:
            return Q.updateNonReversibleBackbone(rng, tuning[NR_BB_SINGLE], 0.5);

        case NR_U:
            return Q.updateNonReversibleU(rng, tuning[NR_U]);

        case JUMP_TO_NR:
            return updateToNonReversible();

        case JUMP_TO_REV:
            return updateToReversible();

        default:
            throw RbException("Unknown sub-move in MPQRateMatrixProposal");
        }
}

void MPQRateMatrixProposal::prepareProposal(void)  {
    
}

/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void MPQRateMatrixProposal::printParameterSummary(std::ostream &o, bool name_only) const  {
    
    o << "lnPriorOdds = ";
    if (name_only == false)
        {
        o << currentLnPriorOdds();
        /* The eigensystem route's own diagnostics. These are only populated when the
           build uses that route (MPQ_TIPROBS_EIGEN); with the default scaling-and-
           squaring route there is nothing to clamp and they stay at zero. */
        if ( RateMatrix_MPQ::getNumClampedTiProbs() > 0 || RateMatrix_MPQ::getWorstRowSumError() > 1e-12 )
            {
            o << ", EIGEN ROUTE: clampedNegatives = " << RateMatrix_MPQ::getNumClampedTiProbs()
              << " (worst " << RateMatrix_MPQ::getWorstNegativeTiProb() << ")"
              << ", worstRowSumError = " << RateMatrix_MPQ::getWorstRowSumError();
            }
        for (int i=0; i<NUM_SUB_MOVES; i++)
            {
            if ( num_dead_intervals[i] > 0 )
                o << ", DEAD INTERVALS " << getSubMoveName( (SUB_MOVE)i ) << " = " << num_dead_intervals[i];
            }
        if ( num_failed_to_reversible > 0 || num_failed_to_nonreversible > 0 || poly.getNumFailures() > 0 )
            {
            o << ", FAILED JUMPS: toReversible = "    << num_failed_to_reversible
              << ", toNonReversible = "               << num_failed_to_nonreversible
              << ", polyhedron = "                    << poly.getNumFailures();
            }
        for (int i=0; i<NUM_SUB_MOVES; i++)
            {
            if ( i == JUMP_TO_NR || i == JUMP_TO_REV )
                continue;
            o << ", " << getSubMoveName( (SUB_MOVE)i ) << " = " << tuning[i]
              << " [acc " << getSubMoveAcceptanceRateTotal( (SUB_MOVE)i )
              << ", n "   << n_tried_total[i] << "]";
            }

        o << ", polyhedronAlpha = " << poly.getAlphaT();
        o << ", polyhedronCenter = (" << poly.getCenter().getX().get_d() << ", "
                                      << poly.getCenter().getY().get_d() << ", "
                                      << poly.getCenter().getZ().get_d() << ")";
        if (poly.getNumCenterClipped() > 0)
            o << ", centerClipped = " << poly.getNumCenterClipped();
        o << ", jumpAcceptance R->N = "  << getSubMoveAcceptanceRateTotal( JUMP_TO_NR )
          << " (n "                      << n_tried_total[JUMP_TO_NR] << ")"
          << ", N->R = "                 << getSubMoveAcceptanceRateTotal( JUMP_TO_REV )
          << " (n "                      << n_tried_total[JUMP_TO_REV] << ")";
        }
}

/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MPQRateMatrixProposal::undoProposal(void)  {
    
    n_rejected[last_sub_move]++;
    n_rejected_total[last_sub_move]++;

#   ifdef MPQ_DEBUG_JUMPS
    dbgReport("rejected");
#   endif

    RateMatrix_MPQ& v = static_cast<RateMatrix_MPQ&>( variable->getValue() );

    v = stored_Q;
}

/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MPQRateMatrixProposal::swapNodeInternal(DagNode *oldN, DagNode *newN) {
    
    variable = static_cast< StochasticNode<RateGenerator>* >(newN) ;
}

void MPQRateMatrixProposal::setProposalTuningParameter(double tp) {

    tuning[PI_ALL] = tp;
}

static double tuneConcentration(double alpha0, double rate, double target) {

    if (rate > target)
        alpha0 /= ( 1.0 + ((rate - target) / (1.0 - target)) );
    else if (rate < 0.1 * target)
        {
        alpha0 *= 4.0;
        }
    else
        alpha0 *= ( 2.0 - rate / target );

    if ( alpha0 < MIN_CONCENTRATION )
        alpha0 = MIN_CONCENTRATION;
    if ( alpha0 > MAX_CONCENTRATION )
        alpha0 = MAX_CONCENTRATION;
    return alpha0;
}

/* Adjust a random-walk window toward a target acceptance rate. This one is a step
   size, so the usual direction applies. */
static double tuneWindow(double delta, double rate, double target) {

    if (rate > target)
        delta *= ( 1.0 + ((rate - target) / (1.0 - target)) );
    else if (rate < 0.1 * target)
        delta /= 4.0;                       // see tuneConcentration
    else
        delta /= (2.0 - rate / target);

    if (delta < MIN_WINDOW)
        delta = MIN_WINDOW;
    if (delta > MAX_WINDOW)
        delta = MAX_WINDOW;
    return delta;
}

void MPQRateMatrixProposal::tune(double rate) {

    tune_calls_total++;

    if ( tune_sub_moves == true )
        {
        for (int i=0; i<NUM_SUB_MOVES; i++)
            {
            // the two jump sub-moves have nothing to tune
            if ( i == JUMP_TO_NR || i == JUMP_TO_REV )
                continue;

            // a sub-move that was never proposed tells us nothing; leave it alone
            if ( n_tried[i] == 0 )
                continue;

            double r = 1.0 - (double)n_rejected[i] / (double)n_tried[i];
            double before = tuning[i];

            bool at_limit;
            if ( i == NR_U )
                at_limit = ( tuning[i] <= 10.0 * MIN_WINDOW );
            else
                at_limit = ( tuning[i] >= 0.1 * MAX_CONCENTRATION );

            if ( n_rejected[i] == n_tried[i] && n_tried[i] >= MIN_TRIES_TO_CALL_DEAD && at_limit == true )
                {
                num_dead_intervals[i]++;
                long p = 1;
                while (p < num_dead_intervals[i])
                    p *= 10;
                if ( p == num_dead_intervals[i] )
                    std::cerr << "MPQRateMatrixProposal: sub-move " << getSubMoveName( (SUB_MOVE)i )
                              << " was refused on all " << n_tried[i] << " attempts in this tuning interval"
                              << " (tuning parameter " << tuning[i] << "; interval " << tune_calls_total
                              << "; this has now happened " << num_dead_intervals[i] << " times)."
                              << " Its tuning parameter is being left alone." << std::endl;
                continue;
                }

            if ( i == NR_U )
                tuning[i] = tuneWindow( tuning[i], r, tune_target[i] );
            else
                tuning[i] = tuneConcentration( tuning[i], r, tune_target[i] );

#           ifdef MPQ_DEBUG_TUNING
            std::cerr << "TUNE interval " << tune_calls_total
                      << "  " << getSubMoveName( (SUB_MOVE)i )
                      << "  tried " << n_tried[i] << "  accepted " << (n_tried[i] - n_rejected[i])
                      << "  rate " << r << "  target " << tune_target[i]
                      << "  tuning " << before << " -> " << tuning[i] << std::endl;
#           endif
            }
        }

    for (int i=0; i<NUM_SUB_MOVES; i++)
        {
        n_tried[i]    = 0;
        n_rejected[i] = 0;
        }

    if ( tune_model_prior == false || n_visits == 0 )
        {
        n_reversible = 0;
        n_visits     = 0;
        return;
        }

    QDistribution* q_dist = dynamic_cast<QDistribution*>( &variable->getDistribution() );
    if ( q_dist == NULL || q_dist->isModelFixed() == true )
        {
        /* Not a dnQ node, or one confined to a single model. Either way there is no
           model indicator to balance, and adapting its prior would be meaningless. */
        n_reversible = 0;
        n_visits     = 0;
        return;
        }

    long n_non_reversible = n_visits - n_reversible;

    if ( n_reversible == 0 || n_non_reversible == 0 )
        {
        /* Search phase. One of the models was never visited, so the observed log
           odds carries no information about HOW far off the prior is, only about
           which way. Step in that direction, doubling while the direction holds and
           halving on a reversal, which brackets the answer geometrically. */
        int dir = ( n_reversible == 0 ) ? 1 : -1;

        if ( tune_last_dir != 0 && dir != tune_last_dir )
            tune_step *= 0.5;
        else
            tune_step *= 2.0;

        if ( tune_step < MIN_TUNE_STEP )
            tune_step = MIN_TUNE_STEP;
        if ( tune_step > MAX_TUNE_STEP )
            tune_step = MAX_TUNE_STEP;

        ln_prior_odds  = q_dist->getLnPriorOdds() + dir * tune_step;
        tune_last_dir  = dir;
#       ifdef MPQ_DEBUG_TUNING
        std::cerr << "TUNE interval " << tune_calls_total << "  lnPriorOdds SEARCH  reversible visits "
                  << n_reversible << " of " << n_visits << "  step " << tune_step
                  << "  tilt -> " << ln_prior_odds << std::endl;
#       endif

        n_tuning_calls = 0;
        }
    else
        {
        // refinement phase
        n_tuning_calls++;

        double p_rev = (double)n_reversible / (double)n_visits;
        double observed_ln_odds = log(p_rev) - log(1.0 - p_rev);

        double gamma = 1.0 / sqrt( (double)n_tuning_calls );
        ln_prior_odds = q_dist->getLnPriorOdds() - gamma * observed_ln_odds;
#       ifdef MPQ_DEBUG_TUNING
        std::cerr << "TUNE interval " << tune_calls_total << "  lnPriorOdds REFINE  reversible visits "
                  << n_reversible << " of " << n_visits << "  p_R " << p_rev
                  << "  tilt -> " << ln_prior_odds << std::endl;
#       endif
        }

    // a hard stop, in case the data are so decisive that no reachable tilt evens
    // the chain out; better a pinned prior that is reported than a runaway one
    if ( ln_prior_odds >  MAX_LN_PRIOR_ODDS )
        ln_prior_odds =  MAX_LN_PRIOR_ODDS;
    if ( ln_prior_odds < -MAX_LN_PRIOR_ODDS )
        ln_prior_odds = -MAX_LN_PRIOR_ODDS;

    q_dist->setLnPriorOdds( ln_prior_odds );

    variable->touch();
    variable->getLnProbability();
    variable->keep();

    n_reversible = 0;
    n_visits     = 0;
}

const char* MPQRateMatrixProposal::getSubMoveName(SUB_MOVE m) {

    switch (m)
        {
        case PI_ALL:        return "pi(all)";
        case PI_SINGLE:     return "pi(one)";
        case REV_BB_ALL:    return "revBackbone(all)";
        case REV_BB_SINGLE: return "revBackbone(one)";
        case NR_BB_ALL:     return "nonrevBackbone(all)";
        case NR_BB_SINGLE:  return "nonrevBackbone(one)";
        case NR_U:          return "polyhedronPoint";
        case JUMP_TO_NR:    return "jumpToNonreversible";
        case JUMP_TO_REV:   return "jumpToReversible";
        default:            return "unknown";
        }
}

double MPQRateMatrixProposal::getSubMoveAcceptanceRate(SUB_MOVE m) const {

    if (n_tried[m] == 0)
        return -1.0;
    return 1.0 - (double)n_rejected[m] / (double)n_tried[m];
}

double MPQRateMatrixProposal::getSubMoveAcceptanceRateTotal( SUB_MOVE m ) const {

    if (n_tried_total[m] == 0)
        return -1.0;
    return 1.0 - (double)n_rejected_total[m] / (double)n_tried_total[m];
}

double MPQRateMatrixProposal::getProportionReversible(void) const {

    if ( n_visits == 0 )
        return 0.0;
    return (double)n_reversible / (double)n_visits;
}

double MPQRateMatrixProposal::updateToNonReversible(void) {
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == false)
        throw RbException("Can only move to a non-reversible model from a reversible one");

    // store old values
    stored_Q = Q;

    // calculate the Jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double wCG = piC * stored_Q(C,G).get_d();
    double wCT = piC * stored_Q(C,T).get_d();
    double wGT = piG * stored_Q(G,T).get_d();

    /* The Jacobian takes the logarithm of three of the backbone weights, so a
       weight of zero would make it negative infinity and the jump would be refused
       every time, silently. Refuse it explicitly instead, and say so. */
    if (wCG <= 0.0 || wCT <= 0.0 || wGT <= 0.0)
        {
        num_failed_to_nonreversible++;
        reportProposalFailure("cannot move to the non-reversible model: one of the C-G, C-T or G-T backbone weights is zero, so the Jacobian is undefined", num_failed_to_nonreversible);
        return RbConstants::Double::neginf;
        }

    double lnJacobian = (log(64.0) + log(wCG) + log(wCT) + log(wGT));
    //lnJacobian -= log(piC) + 2.0 * log(piG);

    // polyhedron density for forward move
    Q.calculateWeights(W);
    Vector randomPoint;
    double lnFwd = poly.lnProbabilityForward(W, randomPoint);
    if (std::isfinite(lnFwd) == false)
        {
        num_failed_to_nonreversible++;
        reportProposalFailure("cannot move to the non-reversible model: the polyhedron could not supply a proposal density", num_failed_to_nonreversible);
        return RbConstants::Double::neginf;
        }
    double lnRv = -lnFwd;
    
    // update rates (off diagonal components)
    mpq_class u1 = randomPoint.getX();
    mpq_class u2 = randomPoint.getY();
    mpq_class u3 = randomPoint.getZ();
#   if 1
    Q.nonreversibilize(u1, u2, u3);
#   else
    Q.setIsReversible(false);
    Q(A,C) = stored_Q(A,C) + (pi[C]/pi[A]) * stored_Q(C,G) * (2 * u1 - 1) + (pi[C]/pi[A]) * stored_Q(C,T) * (2 * u2 - 1);
    Q(A,G) = stored_Q(A,G) - (pi[C]/pi[A]) * stored_Q(C,G) * (2 * u1 - 1) + (pi[G]/pi[A]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(A,T) = stored_Q(A,T) - (pi[C]/pi[A]) * stored_Q(C,T) * (2 * u2 - 1) - (pi[G]/pi[A]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(C,G) = 2 * stored_Q(C,G) * u1;
    Q(C,T) = 2 * stored_Q(C,T) * u2;
    Q(G,T) = 2 * stored_Q(G,T) * u3;
    Q(C,A) = (pi[A]/pi[C]) * stored_Q(A,C) - stored_Q(C,G) * (2 * u1 - 1) - stored_Q(C,T) * (2 * u2 - 1);
    Q(G,A) = (pi[A]/pi[G]) * stored_Q(A,G) + (pi[C]/pi[G]) * stored_Q(C,G) * (2 * u1 - 1) - stored_Q(G,T) * (2 * u3 - 1);
    Q(T,A) = (pi[A]/pi[T]) * stored_Q(A,T) + (pi[C]/pi[T]) * stored_Q(C,T) * (2 * u2 - 1) + (pi[G]/pi[T]) * stored_Q(G,T) * (2 * u3 - 1);
    Q(G,C) = 2 * (pi[C]/pi[G]) * stored_Q(C,G) * (1 - u1);
    Q(T,C) = 2 * (pi[C]/pi[T]) * stored_Q(C,T) * (1 - u2);
    Q(T,G) = 2 * (pi[G]/pi[T]) * stored_Q(G,T) * (1 - u3);
    
    // update the diagonal components of the rate matrix
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
#   endif
    
    double lnProb = lnRv + lnJacobian;
#   ifdef MPQ_DEBUG_JUMPS
    dbg_direction   = 1;
    dbg_ln_jacobian = lnJacobian;
    dbg_ln_poly     = lnRv;
    dbg_ln_hastings = lnProb;
    dbg_u[0]        = u1.get_d();
    dbg_u[1]        = u2.get_d();
    dbg_u[2]        = u3.get_d();
#   endif
    if (std::isfinite(lnProb) == false)
        {
        num_failed_to_nonreversible++;
        reportProposalFailure("cannot move to the non-reversible model: the Hastings ratio is not finite", num_failed_to_nonreversible);
        return RbConstants::Double::neginf;
        }
    return lnProb;
}

double MPQRateMatrixProposal::updateToReversible(void) {
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());
    
    if (Q.getIsReversible() == true)
        throw RbException("Can only move to a time reversible model from a non-reversible one");

    // store old values
    stored_Q = Q;

    // average rates with same stationary frequencies
#   if 1
    Q.reversibilize();
#   else
    Q.setIsReversible(true);
    std::vector<mpq_class>& pi = Q.getPi();
    mpq_class averageRate;
    mpq_class sum;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                Q(i,j) = (pi[i] * stored_Q(i,j) + pi[j] * stored_Q(j,i)) / (2 * pi[i]);
                sum += Q(i,j);
                averageRate += pi[i] * Q(i,j);
                }
            }
        Q(i,i) = -sum;
        }
        
    Q.setExchangeabilityRates();
    
    // make certain average rate is one
    if (averageRate != 1)
        {
        mpq_class factor = 1 / averageRate;
        std::cout << "Warning: the average rate should be one in updateToReversible" << std::endl;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                Q(i,j) *= factor;
        }
#   endif

    // jacobian
    std::vector<mpq_class>& pi = Q.getPi();
    double piC = pi[C].get_d();
    double piG = pi[G].get_d();
    double piT = pi[T].get_d();
    double wCG = piC * stored_Q(C,G).get_d();
    double wGC = piG * stored_Q(G,C).get_d();
    double wCT = piC * stored_Q(C,T).get_d();
    double wTC = piT * stored_Q(T,C).get_d();
    double wGT = piG * stored_Q(G,T).get_d();
    double wTG = piT * stored_Q(T,G).get_d();
    //double lnJacobian = (piC * piG * piG) / (8.0 * (wCG + wGC) * (wCT + wTC) * (wGT + wTG));
    //double lnJacobian = log(piC) + 2.0 * log(piG);
    //lnJacobian -= (log(8.0) + log(wCG + wGC) + log(wCT + wTC) + log(wGT + wTG));
    if (wCG + wGC <= 0.0 || wCT + wTC <= 0.0 || wGT + wTG <= 0.0)
        {
        num_failed_to_reversible++;
        reportProposalFailure("cannot move to the time-reversible model: one of the C-G, C-T or G-T weight pairs sums to zero, so the Jacobian is undefined", num_failed_to_reversible);
        return RbConstants::Double::neginf;
        }

    double lnJacobian = -(log(8.0) + log(wCG + wGC) + log(wCT + wTC) + log(wGT + wTG));
    
    // polyhedron parameters
    Q.calculateWeights(W);
    
    // figure out where in the (u1,u2,u3) space the non-reversible model lives
    if (Q(C,G) <= 0 || Q(C,T) <= 0 || Q(G,T) <= 0)
        {
        num_failed_to_reversible++;
        reportProposalFailure("cannot move to the time-reversible model: a reversibilized rate is zero, so the point in the polyhedron cannot be recovered", num_failed_to_reversible);
        return RbConstants::Double::neginf;
        }
    mpq_class u1 = stored_Q(C,G) / (2 * Q(C,G));
    mpq_class u2 = stored_Q(C,T) / (2 * Q(C,T));
    mpq_class u3 = stored_Q(G,T) / (2 * Q(G,T));
    Vector pt(u1, u2, u3);
    double lnRv = poly.lnProbabilityReverse(W, pt);

    double lnProb = lnRv + lnJacobian;
#   ifdef MPQ_DEBUG_JUMPS
    dbg_direction   = -1;
    dbg_ln_jacobian = lnJacobian;
    dbg_ln_poly     = lnRv;
    dbg_ln_hastings = lnProb;
    dbg_u[0]        = u1.get_d();
    dbg_u[1]        = u2.get_d();
    dbg_u[2]        = u3.get_d();
#   endif
    if (std::isfinite(lnProb) == false)
        {
        num_failed_to_reversible++;
        reportProposalFailure("cannot move to the time-reversible model: the Hastings ratio is not finite", num_failed_to_reversible);
        return RbConstants::Double::neginf;
        }
    return lnProb;
}



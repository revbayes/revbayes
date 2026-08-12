#include "MPQRateMatrixProposal.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RateMatrix_MPQ.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "QDistribution.h"


#define MIN_FREQ    10e-4
#define MAX_LN_PRIOR_ODDS  100.0
#define MIN_CONCENTRATION  1.0
#define MAX_CONCENTRATION  1.0e7
#define MIN_WINDOW         1.0e-4
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
MPQRateMatrixProposal::MPQRateMatrixProposal( StochasticNode<RateGenerator> *n ) : Proposal(),
variable( n )
{
    
    /* Starting values. These matter much less than they used to, because each one
       is now tuned against its own acceptance rate, but a sensible start shortens
       burnin. The concentrations are large because a Dirichlet centred on the
       current value with a large concentration is a small step. */
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
       target for a move on many at once, and is what RevBayes uses for its own
       multi-dimensional proposals. */
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
        n_tried[i]    = 0;
        n_rejected[i] = 0;
        }
    last_sub_move  = PI_ALL;
    tune_sub_moves = true;

    tune_model_prior = false;
    ln_prior_odds    = 0.0;
    n_reversible     = 0;
    n_visits         = 0;
    n_tuning_calls   = 0;
    
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
stored_Q( p.stored_Q ),
W( p.W )
{
    for (int i=0; i<NUM_SUB_MOVES; i++)
        {
        tuning[i]      = p.tuning[i];
        tune_target[i] = p.tune_target[i];
        n_tried[i]     = p.n_tried[i];
        n_rejected[i]  = p.n_rejected[i];
        }
        
    // tell the base class to add the node
    addNode( variable );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MPQRateMatrixProposal::cleanProposal( void ) {

    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MPQRateMatrixProposal* MPQRateMatrixProposal::clone( void ) const {
    
    return new MPQRateMatrixProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MPQRateMatrixProposal::getProposalName( void ) const {

    static std::string name = "MPQRateMatrixProposal";
    return name;
}


double MPQRateMatrixProposal::getProposalTuningParameter( void ) const
{
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
double MPQRateMatrixProposal::doProposal( void ) 
{
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());

    /* Record which model the chain is in, for the prior adaptation in tune().

       Counting here rather than once per generation is deliberate and unbiased:
       the move schedule picks moves in proportion to their weights, independent
       of the state, so the states seen on entry to this function are a thinned
       sample of the states the chain occupies. */
    n_visits++;
    if (Q.getIsReversible() == true)
        n_reversible++;

    /* Choose the sub-move. The probability of proposing a jump is 0.10 from
       either model, so that term cancels in the acceptance ratio. Keep it that
       way, or put the ratio back in explicitly. */
    if (Q.getIsReversible() == true)
        {
        // the time-reversible model: stationary frequencies at fixed weights, or
        // the six backbone weights at fixed frequencies
        double u = rng->uniform01();
        if (u < 0.45)
            last_sub_move = ( rng->uniform01() < 0.25 ) ? PI_ALL : PI_SINGLE;
        else if (u < 0.90)
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
        if (u < 0.30)
            last_sub_move = ( rng->uniform01() < 0.25 ) ? PI_ALL : PI_SINGLE;
        else if (u < 0.60)
            last_sub_move = ( rng->uniform01() < 0.166 ) ? NR_BB_ALL : NR_BB_SINGLE;
        else if (u < 0.90)
            last_sub_move = NR_U;
        else
            last_sub_move = JUMP_TO_REV;
        }

    n_tried[last_sub_move]++;
    return applySubMove( last_sub_move );
}


/* Carry out one sub-move. Every one of them stores the current rate matrix first,
   because undoProposal restores from that copy. */
double MPQRateMatrixProposal::applySubMove( SUB_MOVE m )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(variable->getValue());

    if ( m != JUMP_TO_NR && m != JUMP_TO_REV )
        {
        // the two jump sub-moves store the matrix themselves, further down
        stored_Q = Q;
        }

    switch ( m )
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
            /* Not free: the state is the eight free weights and the move is
               expressed in (backbone, u) coordinates, so it carries the
               reparameterization Jacobian 64 w_CG w_CT w_GT, applied inside
               RateMatrix_MPQ. */
            return Q.updateNonReversibleBackbone(rng, tuning[NR_BB_ALL], 0.5);

        case NR_BB_SINGLE:
            return Q.updateNonReversibleBackbone(rng, tuning[NR_BB_SINGLE], 0.5);

        case NR_U:
            /* The target does not depend on (u1,u2,u3) at fixed backbone, so a
               symmetric window proposal needs no Hastings term at all. */
            return Q.updateNonReversibleU(rng, tuning[NR_U]);

        case JUMP_TO_NR:
            return updateToNonReversible();

        case JUMP_TO_REV:
            return updateToReversible();

        default:
            throw RbException("Unknown sub-move in MPQRateMatrixProposal");
        }
}


/**
 *
 */
void MPQRateMatrixProposal::prepareProposal( void ) 
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void MPQRateMatrixProposal::printParameterSummary(std::ostream &o, bool name_only) const 
{
    
    /* The tilt applied to the model prior has to be reported, because the Bayes
       factor is only recoverable from the sampled model frequencies together with
       the prior odds that produced them:

           BF_NR = (p_N / p_R) * (rho_R / rho_N) = (p_N / p_R) * exp(lnPriorOdds)

       and with equal priors P(non-reversible | X) = BF_NR / (1 + BF_NR). */
    o << "lnPriorOdds = ";
    if (name_only == false)
    {
        o << ln_prior_odds;
        for (int i=0; i<NUM_SUB_MOVES; i++)
            {
            if ( i == JUMP_TO_NR || i == JUMP_TO_REV )
                continue;
            o << ", " << getSubMoveName( (SUB_MOVE)i ) << " = " << tuning[i];
            }
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MPQRateMatrixProposal::undoProposal( void ) 
{
    
    /* This is called on rejection and only on rejection, which is what makes the
       per-sub-move acceptance rates in tune() possible: the rejection belongs to
       whichever sub-move doProposal last chose. Counting acceptances directly is
       not an option, since cleanProposal is not guaranteed to be called only on
       acceptance, so acceptance is obtained by subtracting from the number tried. */
    n_rejected[last_sub_move]++;

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


void MPQRateMatrixProposal::setProposalTuningParameter(double tp)
{
    tuning[PI_ALL] = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
/* Adjust a Dirichlet concentration toward a target acceptance rate.

   The direction is the reverse of the textbook rule, and this is the easy thing
   to get backwards. The textbook rule is written for a step size, where a larger
   value means bigger steps and hence a lower acceptance rate. A Dirichlet
   concentration is the other way round: proposing from Dirichlet(alpha0 * x) puts
   more mass near the current point x as alpha0 grows, so a larger concentration
   means SMALLER steps and a HIGHER acceptance rate. Accepting too often therefore
   calls for a smaller concentration, not a larger one.

   Getting this backwards does not produce an obvious failure. The chain still
   targets the right distribution, since tuning changes only the proposal; it just
   drifts to a concentration where nothing is ever accepted, or where every
   proposal is a step of no consequence, and mixing quietly collapses. */
static double tuneConcentration(double alpha0, double rate, double target)
{
    if ( rate > target )
        alpha0 /= ( 1.0 + ((rate - target) / (1.0 - target)) );
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
static double tuneWindow(double delta, double rate, double target)
{
    if ( rate > target )
        delta *= ( 1.0 + ((rate - target) / (1.0 - target)) );
    else
        delta /= ( 2.0 - rate / target );

    if ( delta < MIN_WINDOW )
        delta = MIN_WINDOW;
    if ( delta > MAX_WINDOW )
        delta = MAX_WINDOW;
    return delta;
}


void MPQRateMatrixProposal::tune( double rate ) {

    /* The rate passed in is ignored on purpose. It is the acceptance rate of this
       move as a whole, and this move is a mixture of nine sub-moves with quite
       different behaviour: an aggregate of 0.4 is equally consistent with every
       sub-move sitting at 0.4 and with half of them at 0.8 and half at 0. Each
       sub-move is tuned against its own rate instead, accumulated in doProposal
       and undoProposal. */

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

            if ( i == NR_U )
                tuning[i] = tuneWindow( tuning[i], r, tune_target[i] );
            else
                tuning[i] = tuneConcentration( tuning[i], r, tune_target[i] );
            }
        }

    for (int i=0; i<NUM_SUB_MOVES; i++)
        {
        n_tried[i]    = 0;
        n_rejected[i] = 0;
        }

    /* ------------------------------------------------------------------
       Adaptation of the prior on the model indicator.
       ------------------------------------------------------------------ */

    if ( tune_model_prior == false || n_visits == 0 )
        {
        n_reversible = 0;
        n_visits     = 0;
        return;
        }

    QDistribution* q_dist = dynamic_cast<QDistribution*>( &variable->getDistribution() );
    if ( q_dist == NULL )
        {
        // not a dnQ node, so there is no model prior here to adapt
        n_reversible = 0;
        n_visits     = 0;
        return;
        }

    n_tuning_calls++;

    /* Observed odds of the two models. The half-counts keep the logit finite when
       one of the models went unvisited over the interval, which is exactly the
       situation the tuning exists to escape from and so must not divide by zero. */
    double p_rev = (n_reversible + 0.5) / (n_visits + 1.0);
    double observed_ln_odds = log(p_rev) - log(1.0 - p_rev);

    /* Since the posterior odds are the Bayes factor times the prior odds, and the
       target is even odds, the exactly correcting step is

           delta_new = delta_old - log(p_R / p_N)

       which would land on the answer in one move if p_R were known rather than
       estimated. It is estimated, so damp by 1/sqrt(k): the first call takes the
       full Newton step, which is what gets a badly tilted chain back into range
       quickly, and later calls average out the noise in p_R rather than chasing
       it. This is the usual Robbins-Monro schedule and it is why adaptation has
       to stop before sampling begins. */
    double gamma = 1.0 / sqrt( (double)n_tuning_calls );
    ln_prior_odds = q_dist->getLnPriorOdds() - gamma * observed_ln_odds;

    // a hard stop, in case the data are so decisive that no reachable tilt evens
    // the chain out; better a pinned prior that is reported than a runaway one
    if ( ln_prior_odds >  MAX_LN_PRIOR_ODDS )
        ln_prior_odds =  MAX_LN_PRIOR_ODDS;
    if ( ln_prior_odds < -MAX_LN_PRIOR_ODDS )
        ln_prior_odds = -MAX_LN_PRIOR_ODDS;

    q_dist->setLnPriorOdds( ln_prior_odds );

    n_reversible = 0;
    n_visits     = 0;
}


const char* MPQRateMatrixProposal::getSubMoveName( SUB_MOVE m ) {

    switch ( m )
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


double MPQRateMatrixProposal::getSubMoveAcceptanceRate( SUB_MOVE m ) const {

    if ( n_tried[m] == 0 )
        return -1.0;
    return 1.0 - (double)n_rejected[m] / (double)n_tried[m];
}


double MPQRateMatrixProposal::getProportionReversible( void ) const {

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
    double lnJacobian = (log(64.0) + log(wCG) + log(wCT) + log(wGT));
    //lnJacobian -= log(piC) + 2.0 * log(piG);

    // polyhedron density for forward move
    Q.calculateWeights(W);
    Vector randomPoint;
    double lnRv = -poly.lnProbabilityForward(W, randomPoint);
    
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
    
    return lnRv + lnJacobian;
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
    double lnJacobian = -(log(8.0) + log(wCG + wGC) + log(wCT + wTC) + log(wGT + wTG));
    
    // polyhedron parameters
    Q.calculateWeights(W);
    
    // figure out where in the (u1,u2,u3) space the non-reversible model lives
    mpq_class u1 = stored_Q(C,G) / (2 * Q(C,G));
    mpq_class u2 = stored_Q(C,T) / (2 * Q(C,T));
    mpq_class u3 = stored_Q(G,T) / (2 * Q(G,T));
    Vector pt(u1, u2, u3);
    double lnRv = poly.lnProbabilityReverse(W, pt);

    return lnRv + lnJacobian;
}



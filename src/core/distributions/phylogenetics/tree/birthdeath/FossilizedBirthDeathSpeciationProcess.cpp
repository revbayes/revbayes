#include <optional>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevBayesCore;

/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    c              Fossil observation counts.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    a              Anagenetic speciation rates.
 * \param[in]    b              Symmetric speciation probability.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
FossilizedBirthDeathSpeciationProcess::FossilizedBirthDeathSpeciationProcess(const TypedDagNode<double> *ra,
                                                           const DagNode *inspeciation,
                                                           const DagNode *inextinction,
                                                           const DagNode *inpsi,
                                                           const TypedDagNode<double> *inrho,
                                                           const DagNode *inlambda_a,
                                                           const DagNode *inbeta,
                                                           const TypedDagNode< RbVector<double> > *intimes,
                                                           const std::string &incondition,
                                                           const std::vector<Taxon> &intaxa,
                                                           bool comp,
                                                           bool ext) :
    AbstractBirthDeathProcess(ra, incondition, intaxa, true, NULL),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, comp),
    extended( ext )
{

    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    homogeneous_lambda_a             = NULL;
    homogeneous_beta                 = NULL;
    heterogeneous_lambda_a           = NULL;
    heterogeneous_beta               = NULL;

    heterogeneous_lambda_a = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inlambda_a);
    homogeneous_lambda_a   = dynamic_cast<const TypedDagNode<double >*>(inlambda_a);
    heterogeneous_beta     = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inbeta);
    homogeneous_beta       = dynamic_cast<const TypedDagNode<double >*>(inbeta);

    addParameter( homogeneous_lambda_a );
    addParameter( heterogeneous_lambda_a );
    addParameter( homogeneous_beta );
    addParameter( heterogeneous_beta );

    RbException no_timeline_err = RbException("No time intervals provided for heterogeneous fossilized birth death process");

    RbException inconsistent_rates = RbException("Inconsistent number of rates in fossilized birth death process.");

    if ( heterogeneous_lambda_a != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_lambda_a->getValue().size() != num_intervals ) throw inconsistent_rates;
    }
    if ( heterogeneous_beta != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_beta->getValue().size() != num_intervals ) throw inconsistent_rates;
    }

    species        = std::vector<SpeciesEntry>(taxa.size());
    next_species   = std::vector<SpeciesEntry>(taxa.size());
    stored_species = std::vector<SpeciesEntry>(taxa.size());

    anagenetic    = std::vector<double>(num_intervals, 0.0);
    symmetric     = std::vector<double>(num_intervals, 0.0);

    
    redrawValue();

    // nothing has touched this node yet, so the first density reads what the pull leaves here
    updateRanges();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
FossilizedBirthDeathSpeciationProcess* FossilizedBirthDeathSpeciationProcess::clone( void ) const
{
    return new FossilizedBirthDeathSpeciationProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    // the cast is for the cached per-taxon terms, which this const override has to refresh
    double lnProb = const_cast<FossilizedBirthDeathSpeciationProcess*>(this)->computeLnProbabilityRanges();

    // no tree this process can produce, so reject rather than score it
    if ( invalid_continuation == true ) return RbConstants::Double::neginf;

    lnProb += computeLnProbabilityTimes();

    // every budding event carries 1-beta; the symmetric ones pay where they close a range
    lnProb += budding_lnProb;

    return lnProb;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityTimes( void ) const
{
    double lnProb = 0.0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        // the parent species is a sampled ancestor, which means different things in the two trees
        if ( species[i].parent_is_sa == true )
        {
            double y_a  = ranges[i].birth;
            size_t y_ai = findIndex(y_a);

            if ( extended == true )
            {
                // the ancestor's range ended by anagenetic speciation, not by extinction
                lnProb -= log( birth[y_ai] );
                lnProb -= log( death[y_ai] );
                lnProb += log( anagenetic[y_ai] );
            }
            else
            {
                // the speciation separating the two ranges is unobserved, so integrate it over
                // (y_a, first_min): the bracket is P(at least one species change along it)
                double o  = ranges[i].first;
                size_t oi = findIndex(o);

                // first_min is younger than y_a, so the intermediate terms subtract here
                double ln_q  = q(oi, o) - q(y_ai, y_a);
                double ln_qt = q(oi, o, true) - q(y_ai, y_a, true);
                for (size_t j = oi; j < y_ai; j++)
                {
                    ln_q  -= q_i[j];
                    ln_qt -= q_tilde_i[j];
                }

                lnProb -= log( birth[y_ai] );
                lnProb += log( 1.0 - exp( ln_q - ln_qt ) );
            }
        }
    }

    return lnProb;
}


/** Close range i at d: an extinction, a symmetric speciation, or the marginalization limit. */
double FossilizedBirthDeathSpeciationProcess::rangeEndTerm( size_t i, size_t di, double d ) const
{
    // the event contributes lambda*beta, and its daughters already carry a lambda each
    if ( species[i].ends_symmetric == true ) return log( symmetric[di] ) - log( birth[di] );

    if ( extended == true ) return log( death[di] );

    if ( species[i].ends_at_sa == true ) return 0.0;

    return log( p( di, d ) );
}


double FossilizedBirthDeathSpeciationProcess::getMaxTaxonAge( const TopologyNode& node ) const
{
    if( node.isTip() )
    {
        return ranges[node.getIndex()].first;
    }
    else
    {
        double max = 0;

        for( size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            max = std::max( getMaxTaxonAge( node.getChild(i) ), max );
        }

        return max;
    }
}


double FossilizedBirthDeathSpeciationProcess::lnProbTreeShape(void) const
{
    // the divergence-time density is for an unlabeled oriented tree, so convert to a labeled one
    // by multiplying by 1 / n!

    return - RbMath::lnFactorial( value->getNumberOfExtantTips() );
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Probability of survival.
 */
double FossilizedBirthDeathSpeciationProcess::pSurvival(double start, double end) const
{
    return AbstractFossilizedBirthDeathRangeProcess::p(findIndex(start), start, true);
}


/**
 * q_i(t)
 */
double FossilizedBirthDeathSpeciationProcess::q( size_t i, double t, bool tilde ) const
{

    if ( t == times[i] ) return 0.0;

    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];

    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    double q = log(4.0) + ln_e - 2.0*log(tmp);

    if (tilde)
    {
        q = 0.5 * (q - (b+d+f)*dt);

        double a = anagenetic[i];
        double s = symmetric[i];

        q = - a * dt - s * (b + d + f) * dt + (1.0 - s) * q;
    }

    return q;
}


/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::redrawValue(SimulationCondition c)
{
    redrawValue();
}


/** Set the tree, adopting its annotations where it carries them and drawing the rest. */
void FossilizedBirthDeathSpeciationProcess::setValue(Tree *v, bool force)
{
    AbstractBirthDeathProcess::setValue(v, force);

    // a tree built here indexes its tips by taxon; one set from outside carries whatever order it
    // was written in, and the per-taxon terms are read by index
    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for (size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
    }
    this->getValue().orderNodesByIndex();

    // before the ages either way: the flags name each origination time, which bounds them
    if ( adoptSpeciesLabels() == false )
    {
        normalizeContinuationFlags();
    }

    prepareProbComputation();
    updateRanges();

    // drawn against the tree the constructor built, so draw them again inside this one
    if ( adoptAppearances() == false )
    {
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            initializeFirstLast(i);
        }
    }
}


/** The tree annotated with the range labels and appearances a plain Newick drops. */
std::string FossilizedBirthDeathSpeciationProcess::getHiddenStateString( void ) const
{
    Tree t = *this->value;

    labelSpecies( t.getRoot(), continuingSpecies( t.getRoot() ) );

    for (size_t i = 0; i < taxa.size(); ++i)
    {
        TopologyNode &tip = t.getNode(i);

        std::stringstream first, last;
        first.precision( 17 );
        last.precision( 17 );
        first << ranges[i].first;
        last  << ranges[i].last;

        tip.setNodeParameter( "FAD", first.str() );
        tip.setNodeParameter( "LAD", last.str() );
    }

    return t.getNewickRepresentation( false );
}


void FossilizedBirthDeathSpeciationProcess::setHiddenStateFromString( const std::string &s )
{
    // through the value's own type, so the Newick is read back as a time tree
    Tree *t = this->value->clone();
    t->initFromString( s );

    setValue( t );
}


/** The range each node belongs to, stamped top down. */
void FossilizedBirthDeathSpeciationProcess::labelSpecies( TopologyNode &node, int s ) const
{
    if ( s >= 0 ) node.setNodeParameter( "range", StringUtilities::to_string(s) );

    std::vector<TopologyNode*> children = node.getChildren();
    for (size_t c = 0; c < children.size(); c++)
    {
        // a child that does not continue begins its own species at this node
        int cs = children[c]->continuesParentSpecies() ? s : continuingSpecies( *children[c] );

        labelSpecies( *children[c], cs );
    }
}


/** A child continues its parent's range exactly when the two labels match. */
bool FossilizedBirthDeathSpeciationProcess::adoptSpeciesLabels( void )
{
    Tree &t = this->getValue();

    size_t n = t.getNumberOfNodes();

    std::vector<std::string> label( n );
    bool annotated = false;

    for (size_t i = 0; i < n; ++i)
    {
        const TopologyNode &node = t.getNode(i);

        std::optional<std::string> s = node.getNodeParameter( "range" );
        if ( s.has_value() == false ) continue;

        label[ node.getIndex() ] = *s;
        annotated = true;
    }

    if ( annotated == false ) return false;

    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode &node = t.getNode(i);

        node.eraseNodeParameter( "range" );

        if ( node.isRoot() == true ) continue;

        const std::string &mine   = label[ node.getIndex() ];
        const std::string &parent = label[ node.getParent().getIndex() ];

        node.setContinuesParentSpecies( mine.empty() == false && mine == parent );
    }

    return true;
}


bool FossilizedBirthDeathSpeciationProcess::adoptAppearances( void )
{
    Tree &t = this->getValue();

    bool annotated = false;

    for (size_t i = 0; i < taxa.size(); ++i)
    {
        TopologyNode &tip = t.getNode(i);

        std::optional<std::string> f = tip.eraseNodeParameter( "FAD" );
        std::optional<std::string> l = tip.eraseNodeParameter( "LAD" );

        if ( f.has_value() == true )
        {
            ranges[i].first = atof( f->c_str() );
            annotated = true;
        }
        if ( l.has_value() == true )
        {
            ranges[i].last = atof( l->c_str() );
            annotated = true;
        }

        if ( annotated == true ) dirty_taxa[i] = true;
    }

    return annotated;
}


void FossilizedBirthDeathSpeciationProcess::setOccurrences( const std::vector<Taxon> &t )
{
    AbstractFossilizedBirthDeathRangeProcess::setOccurrences( t );

    // the tree distribution keeps its own taxon copy to build tips from, so a later redraw uses
    // the clamped record rather than the one the process was constructed with
    AbstractRootedTreeDistribution::taxa = taxa;
}


void FossilizedBirthDeathSpeciationProcess::setMcmcMode(bool tf)
{
    AbstractBirthDeathProcess::setMcmcMode(tf);
    if ( tf == true ) warnIfNoResampleMove();
    if ( tf == true ) warnIfNoReportingNode();
}


void FossilizedBirthDeathSpeciationProcess::redrawValue(void)
{
    // draw the ranges as the matrix process does, then hang a budding topology on them:
    // conditional on the ranges the tree shape is uniform, so a uniform attachment is a valid draw
    drawRanges();

    RandomNumberGenerator* rng = GLOBAL_RNG;
    size_t n = taxa.size();

    // the tip is the extinction and must sit below the oldest appearance, which independent
    // draws can violate under age uncertainty
    double present = times.front();
    for (size_t i = 0; i < n; ++i)
    {
        if ( ranges[i].death > ranges[i].first ) ranges[i].death = rng->uniform01() * (ranges[i].first - present) + present;
    }

    // one tip per taxon at its drawn end age; top[i] tracks the current top of lineage i's subtree
    std::vector<TopologyNode*> top(n);
    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode* tip = new TopologyNode( taxa[i], i );
        tip->setTipAgeUnconstrained( true );
        tip->setAge( ranges[i].death );
        top[i] = tip;
    }

    // independent births need not form a tree, so draw them jointly: oldest first, each inside
    // the span of a still-living lineage already placed
    double origin = getOriginAge();

    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    std::sort( order.begin(), order.end(), [this](size_t a, size_t b){ return ranges[a].first > ranges[b].first; } );

    size_t root_lineage = order[0];
    std::vector<size_t> parent(n, n);   // parent[k] = lineage k buds off; n marks the origin
    ranges[root_lineage].birth = origin;
    for (size_t idx = 1; idx < n; ++idx)
    {
        size_t k = order[idx];

        // already-placed lineages that are old enough to bud k
        std::vector<size_t> cand;
        for (size_t p = 0; p < idx; ++p)
        {
            if ( ranges[order[p]].birth > ranges[k].first ) cand.push_back( order[p] );
        }

        // an origin younger than k's oldest age leaves none, so bud from the origin lineage and
        // let the density reject the birth, which lets the caller draw another origin
        size_t j = root_lineage;
        if ( cand.empty() == false )
        {
            size_t pick = size_t( rng->uniform01() * cand.size() );
            if ( pick >= cand.size() ) pick = cand.size() - 1;

            j = cand[pick];
        }

        // birth in (max(first_k, d_j), b_j): after j is born and while it is still alive. The
        // birth must stay below b_j whatever the window, or the tree is assembled out of order.
        double lo = std::max( ranges[k].first, ranges[j].death );
        if ( lo >= ranges[j].birth ) lo = 0.0;

        ranges[k].birth = rng->uniform01()*(ranges[j].birth - lo) + lo;
        parent[k] = j;
    }

    // build the tree youngest birth first, so each lineage's subtree is complete before it attaches
    std::vector<size_t> byBirth(n);
    for (size_t i = 0; i < n; ++i) byBirth[i] = i;
    std::sort( byBirth.begin(), byBirth.end(), [this](size_t a, size_t b){ return ranges[a].birth < ranges[b].birth; } );

    for (size_t idx = 0; idx < n; ++idx)
    {
        size_t k = byBirth[idx];
        if ( parent[k] == n ) continue;   // the origin lineage stays the root
        size_t j = parent[k];

        // budding node at b_k: lineage j continues, k is the new species. Recorded on the children
        // rather than in their order, so no later move can reassign it by permuting them
        TopologyNode* node = new TopologyNode();
        node->setAge( ranges[k].birth );
        node->addChild( top[j] );
        node->addChild( top[k] );
        top[j]->setParent( node );
        top[k]->setParent( node );
        top[j]->setContinuesParentSpecies( true );
        top[k]->setContinuesParentSpecies( false );
        top[j] = node;
    }

    Tree *psi = new Tree();
    psi->setRooted( true );
    psi->setRoot( top[root_lineage], true );

    delete this->value;
    this->value = psi;

    // index tips by taxon order, as the likelihood expects
    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for( size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
    }
    this->getValue().orderNodesByIndex();

    // whatever path built this tree, leave it with exactly one continuation per node
    normalizeContinuationFlags();
}



/** beta in the interval containing age, read from the parameter rather than the prepared cache. */
double FossilizedBirthDeathSpeciationProcess::symmetricAt( double age ) const
{
    if ( homogeneous_beta != NULL )
    {
        return homogeneous_beta->getValue();
    }

    const RbVector<double>& probs = heterogeneous_beta->getValue();
    size_t i = findIndex( age );

    return i < probs.size() ? probs[i] : 0.0;
}


bool FossilizedBirthDeathSpeciationProcess::hasSymmetricSpeciation( void ) const
{
    if ( homogeneous_beta != NULL )
    {
        return homogeneous_beta->getValue() > 0.0;
    }
    if ( heterogeneous_beta != NULL )
    {
        const RbVector<double>& probs = heterogeneous_beta->getValue();
        for (size_t i = 0; i < probs.size(); ++i)
        {
            if ( probs[i] > 0.0 ) return true;
        }
    }

    return false;
}


bool FossilizedBirthDeathSpeciationProcess::hasAnagenesis( void ) const
{
    if ( homogeneous_lambda_a != NULL )
    {
        return homogeneous_lambda_a->getValue() > 0.0;
    }
    if ( heterogeneous_lambda_a != NULL )
    {
        const RbVector<double>& rates = heterogeneous_lambda_a->getValue();
        for (size_t i = 0; i < rates.size(); ++i)
        {
            if ( rates[i] > 0.0 ) return true;
        }
    }

    return false;
}


/**
 * Redraw the budding topology with the ranges fixed. Every compatible tree has the same density,
 * so this is a Gibbs step. Returns false where the premise fails and leaves the tree alone.
 */
bool FossilizedBirthDeathSpeciationProcess::redrawTopology( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // an anagenetic attachment sits at ranges[j].death == ranges[k].birth, which the candidate test below cannot reach
    if ( hasAnagenesis() == true ) return false;

    // the draw is pure budding, so under beta > 0 it does not target the conditional: the current
    // tree may hold symmetric nodes the replacement cannot
    if ( hasSymmetricSpeciation() == true ) return false;

    updateRanges();

    size_t n = taxa.size();
    if ( n < 2 ) return false;

    size_t root_lineage = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if ( ranges[i].birth > ranges[root_lineage].birth ) root_lineage = i;
    }

    // choose every attachment before touching the tree, so an impossible draw costs nothing
    std::vector<size_t> parent(n, n);
    for (size_t k = 0; k < n; ++k)
    {
        if ( k == root_lineage ) continue;

        std::vector<size_t> cand;
        for (size_t j = 0; j < n; ++j)
        {
            if ( j != k && ranges[j].birth > ranges[k].birth && ranges[j].death < ranges[k].birth ) cand.push_back( j );
        }

        if ( cand.empty() == true ) return false;

        size_t pick = size_t( rng->uniform01() * cand.size() );
        if ( pick >= cand.size() ) pick = cand.size() - 1;

        parent[k] = cand[pick];
    }

    std::vector<TopologyNode*> top(n);
    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode* tip = new TopologyNode( taxa[i], i );
        tip->setTipAgeUnconstrained( true );
        tip->setAge( ranges[i].death );
        top[i] = tip;
    }

    // youngest birth first, so a lineage's subtree is complete before it attaches
    std::vector<size_t> byBirth(n);
    for (size_t i = 0; i < n; ++i) byBirth[i] = i;
    std::sort( byBirth.begin(), byBirth.end(), [this](size_t a, size_t b){ return ranges[a].birth < ranges[b].birth; } );

    for (size_t idx = 0; idx < n; ++idx)
    {
        size_t k = byBirth[idx];
        if ( parent[k] == n ) continue;

        size_t j = parent[k];

        // budding node at b_k: lineage j continues, k is the new species. Recorded on the children
        // rather than in their order, so no later move can reassign it by permuting them
        TopologyNode* node = new TopologyNode();
        node->setAge( ranges[k].birth );
        node->addChild( top[j] );
        node->addChild( top[k] );
        top[j]->setParent( node );
        top[k]->setParent( node );
        top[j]->setContinuesParentSpecies( true );
        top[k]->setContinuesParentSpecies( false );
        top[j] = node;
    }

    Tree *psi = new Tree();
    psi->setRooted( true );
    psi->setRoot( top[root_lineage], true );

    delete this->value;
    this->value = psi;

    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for( size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
    }
    this->getValue().orderNodesByIndex();

    // whatever path built this tree, leave it with exactly one continuation per node
    normalizeContinuationFlags();

    return true;
}

/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn)
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the minimum birth age
    std::vector<double> first_occurrences;

    double current_age = RbConstants::Double::inf;
    double minimum_age = 0.0;
    double max_age = getOriginAge();

    for (size_t i = 0; i < n.size(); ++i)
    {
        // make sure the tip age is equal to the last occurrence
        if( n[i]->isTip() )
        {
            bool extinct = n[i]->getTaxon().isExtinct();

            double present = times.front();

            n[i]->setAge( extinct * rng->uniform01() * (n[i]->getTaxon().getMinAge() - present) + present );

            size_t j = find(taxa.begin(), taxa.end(), n[i]->getTaxon()) - taxa.begin();

            double minmax = std::max(ranges[i].first_min, n[i]->getAge());
            ranges[j].first = GLOBAL_RNG->uniform01()*(std::min(max_age, taxa[i].getMaxAge()) - minmax) + minmax;
        }

        double first_occurrence = getMaxTaxonAge( *n[i] );

        if( first_occurrence > minimum_age )
        {
            minimum_age = first_occurrence;
        }

        first_occurrences.push_back( first_occurrence );

        if ( current_age > n[i]->getAge() )
        {
            current_age = n[i]->getAge();
        }
    }

    if( minimum_age > max_age )
    {
        throw RbException() << "Tree age is " << max_age << " but minimum fossil origin is " << minimum_age;
    }


    if ( age <= minimum_age )
    {
        age = rng->uniform01() * ( max_age - minimum_age ) + minimum_age;
    }


    std::vector<double> ages;
    while ( n.size() > 2 && current_age < age )
    {

        // get all the nodes with first occurrences younger than the current age
        std::vector<TopologyNode*> active_nodes;
        std::vector<TopologyNode*> active_right_nodes;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age >= n[i]->getAge() )
            {
                active_nodes.push_back( n[i] );

                if( current_age >= first_occurrences[i] )
                {
                    active_right_nodes.push_back( n[i] );
                }
            }

        }

        // we need to get the next node age older than the current age
        double next_node_age = age;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age < n[i]->getAge() && n[i]->getAge() < next_node_age )
            {
                next_node_age = n[i]->getAge();
            }
            if ( current_age < first_occurrences[i] && first_occurrences[i] < next_node_age )
            {
                next_node_age = first_occurrences[i];
            }

        }

        // only simulate if there are at least two valid/active nodes and one active right node
        if ( active_nodes.size() <= 2 || active_right_nodes.empty() )
        {
            current_age = next_node_age;
        }
        else
        {
            // now we simulate new ages
            double next_sim_age = simulateNextAge(active_nodes.size()-2, age, present, current_age, alwaysReturn);

            if ( next_sim_age < next_node_age )
            {
                // randomly pick two nodes
                size_t index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                TopologyNode* left_child = active_nodes[index_left];

                size_t index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                TopologyNode* right_child = active_right_nodes[index_right];

                while( left_child == right_child )
                {
                    index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                    left_child = active_nodes[index_left];

                    index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                    right_child = active_right_nodes[index_right];
                }

                // erase the nodes also from the origin nodes vector
                std::vector<TopologyNode *>::iterator child_it_left = std::find( n.begin(), n.end(), left_child );
                std::vector<double>::iterator fa_it_left = first_occurrences.begin() + std::distance( n.begin(), child_it_left );
                double fa_left = *fa_it_left;

                first_occurrences.erase(fa_it_left);
                n.erase(child_it_left);

                std::vector<TopologyNode *>::iterator child_it_right = std::find( n.begin(), n.end(), right_child );
                std::vector<double>::iterator fa_it_right = first_occurrences.begin() + std::distance( n.begin(), child_it_right );
                double fa_right = *fa_it_right;

                first_occurrences.erase(fa_it_right);
                n.erase(child_it_right);


                // create a parent for the two
                TopologyNode *parent = new TopologyNode();
                parent->addChild( left_child );
                parent->addChild( right_child );
                left_child->setParent( parent );
                right_child->setParent( parent );
                left_child->setContinuesParentSpecies( true );
                right_child->setContinuesParentSpecies( false );
                parent->setAge( next_sim_age );

                // insert the parent to our list
                n.push_back( parent );
                first_occurrences.push_back( std::max(fa_left, fa_right) );

                current_age = next_sim_age;
                ages.push_back( next_sim_age );
            }
            else
            {
                current_age = next_node_age;
            }

        }

        if ( n.size() > 2 && current_age >= age  )
            throw RbException() << "Unexpected number of taxa (remaining #taxa was " << n.size() << " and age was " << current_age << " with maximum age of " << age << ") in tree simulation";

    }


    if ( n.size() == 2 )
    {

        // pick two nodes
        TopologyNode* left_child = n[0];
        TopologyNode* right_child = n[1];

        // make sure the speciation event is older than the new species first occurrence
        if( first_occurrences[1] > age )
        {
            if( first_occurrences[0] > age )
            {
                throw RbException() << "Cannot simulate clade of age " << age << ", minimum age is " << minimum_age;
            }
            else
            {
                std::swap( left_child, right_child );
            }
        }
        else if( age > first_occurrences[0] )
        {
            if( rng->uniform01() < 0.5 )
            {
                std::swap( left_child, right_child );
            }
        }

        // erase the nodes also from the origin nodes vector
        n.clear();

        // create a parent for the two
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left_child );
        parent->addChild( right_child );
        left_child->setParent( parent );
        right_child->setParent( parent );
        // this simulator has no notion of which lineage keeps the ancestral species, so name one:
        // the density needs exactly one continuation per node, and mvRotateNode explores the choice
        left_child->setContinuesParentSpecies( true );
        right_child->setContinuesParentSpecies( false );
        parent->setAge( age );

        // insert the parent to our list
        n.push_back( parent );
    }
    else
    {
        throw RbException() << "Unexpected number of taxa (" << n.size() << ") in tree simulation";
    }


}

/**
 * @param alwaysReturn whether the simulation can return times which are not valid draws from the distribution (for initial values)
*/
std::vector<double> FossilizedBirthDeathSpeciationProcess::simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const
{
    if(!alwaysReturn) throw RbException("Impossible to simulate under a true FossilizedBirthDeathSpeciation process");

    std::vector<double> t(n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        t[i] = simulateDivergenceTime(origin, min);
    }

    // finally sort the times
    std::sort(t.begin(), t.end());

    return t;
}

/**
 * Simulate new speciation times.
 */
double FossilizedBirthDeathSpeciationProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant SSBDP
    // direct forward simulation under FBDRP is not feasible

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = findIndex(present);

    // get the parameters
    double age = origin - present;
    double b = birth[i];
    double d = death[i];
    double p_e = i == 0 ? homogeneous_rho->getValue() : 0.0;


    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(p_e*b+(b*(1-p_e)-d)*exp((d-b)*age) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (exp(age*(d-b)) - 1) / exp(age*(d-b)) ) / (b-d);
        }
    }
    else
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(p_e*b*exp((b-d)*age)+(b*(1-p_e)-d) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (1 - exp(age*(b-d)))  ) / (b-d);
        }
    }

    return present + t;
}


void FossilizedBirthDeathSpeciationProcess::normalizeContinuationFlags( void )
{
    // the legal repair depends on beta at each node, so the interval boundaries have to be current
    prepareProbComputation();

    normalizeContinuationFlags( getValue().getRoot() );
}


/** Make every node name exactly one continuing child. Construction only, never the density. */
/** The taxon whose range this node belongs to, found by following the continuations down. */
int FossilizedBirthDeathSpeciationProcess::continuingSpecies( const TopologyNode &node ) const
{
    if ( node.isTip() == true ) return int( node.getIndex() );

    const std::vector<TopologyNode*> &children = node.getChildren();
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->continuesParentSpecies() == true ) return continuingSpecies( *children[c] );
    }

    return -1;
}


void FossilizedBirthDeathSpeciationProcess::normalizeContinuationFlags( const TopologyNode &node )
{
    if ( node.isTip() == true ) return;

    std::vector<TopologyNode*> children = node.getChildren();

    for (size_t c = 0; c < children.size(); c++)
    {
        normalizeContinuationFlags( *children[c] );
    }

    // a sampled ancestor tip is a sample of this node's species and always continues
    bool sa_tip = false;
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == true )
        {
            sa_tip = true;
            children[c]->setContinuesParentSpecies( true );
        }
    }

    size_t n_cont = 0;
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == false && children[c]->continuesParentSpecies() == true ) ++n_cont;
    }

    // zero continuations is symmetric speciation, legal where this node's own beta is positive.
    // The interval matters here: repairing to one it forbids starts the chain at -inf
    bool legal = ( n_cont == 1 ) ||
                 ( n_cont == 0 && ( sa_tip == true || symmetricAt( node.getAge() ) > 0.0 ) );

    if ( legal == false )
    {
        // every other child's range starts here, and none may be born younger than its own
        // oldest occurrence; lowest index among the allowed, so presentations agree
        size_t keep = children.size();

        for (size_t c = 0; c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;

            bool feasible = true;
            for (size_t o = 0; o < children.size() && feasible == true; o++)
            {
                if ( o == c || children[o]->isSampledAncestorTip() == true ) continue;

                int s = continuingSpecies( *children[o] );
                if ( s < 0 || ranges[s].first_min >= node.getAge() ) feasible = false;
            }

            if ( feasible == false ) continue;
            if ( keep == children.size() || children[c]->getIndex() < children[keep]->getIndex() ) keep = c;
        }

        // no choice fits the record, so take the lowest index and let the density reject it
        for (size_t c = 0; keep == children.size() && c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;
            keep = c;
        }

        for (size_t c = 0; c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;
            children[c]->setContinuesParentSpecies( c == keep );
        }
    }
}


FossilizedBirthDeathSpeciationProcess::RangeFlow FossilizedBirthDeathSpeciationProcess::updateRanges( const TopologyNode& node )
{
    if( node.isTip() )
    {
        return RangeFlow{ int(node.getIndex()), 0.0 };
    }

    RangeFlow flow = { -1, 0.0 };

    std::vector<TopologyNode* > children = node.getChildren();

    bool sa = node.isSampledAncestorTipOrParent();

    // a sampled ancestor tip always continues. Among the rest, one continuing is a budding event
    // and none is symmetric speciation. Read only: a density must not write to the value it scores
    size_t n_cont = 0;
    bool   sa_tip = false;

    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == true )
        {
            sa_tip = true;
            if ( children[c]->continuesParentSpecies() == false ) invalid_continuation = true;
        }
        else if ( children[c]->continuesParentSpecies() == true )
        {
            ++n_cont;
        }
    }

    if ( n_cont > 1 ) invalid_continuation = true;

    // screens against a timeline with no symmetric speciation; the beta in rangeEndTerm settles
    // which interval permits it
    if ( sa_tip == true )
    {
        // the sampled ancestor already carries the species. A sibling that continues it as well is
        // a sample taken within a range that runs on past it, which only ends by symmetric speciation
        if ( n_cont == 1 && hasSymmetricSpeciation() == false ) invalid_continuation = true;
    }
    else
    {
        // no child carries the species, so it ends here by symmetric speciation
        if ( n_cont == 0 && hasSymmetricSpeciation() == false ) invalid_continuation = true;
    }

    // stop before the assignment loop, which would index the range table through species == -1
    if ( invalid_continuation == true ) return flow;

    // a bifurcation that is not symmetric is a budding event; the symmetric ones pay their beta
    // through rangeEndTerm on the range they close
    if ( sa_tip == false && n_cont == 1 )
    {
        budding_lnProb += log( 1.0 - symmetricAt( node.getAge() ) );
    }

    // a species that ends below reaches this node unnamed; the sampled ancestor here names it
    int pending_species = -1;
    double pending_end  = 0.0;

    for(int c = 0; c < children.size(); c++)
    {
        const TopologyNode& child = *children[c];

        RangeFlow sub = updateRanges(child);

        // a subtree that failed propagates up rather than writing through a negative index
        if ( invalid_continuation == true ) return flow;

        // the subtree's species ended below and is still unnamed; carry it past this node
        if ( sub.species < 0 )
        {
            pending_end = sub.end_age;
            continue;
        }

        int i = sub.species;

        // is child a new species?
        // set start time at this node
        if( child.continuesParentSpecies() == false )
        {
            double age = node.getAge(); // y_{a(i)}

            if ( age != ranges[i].birth )
            {
                ranges[i].birth = age;
                dirty_taxa[i] = true;
            }

            next_species[i].parent_is_sa = sa;
        }
        // child is the ancestral species
        else
        {
            // propagate species index
            flow.species = i;

            // a sampled ancestor tip is a sample of this species, so it can name one that ended below
            if ( child.isSampledAncestorTip() == true ) pending_species = i;

            // its range ends here and the lineage carries on below
            if ( sa == true ) next_species[i].ends_at_sa = true;

            // if this is the root
            // set the start time to the origin
            if( node.isRoot() )
            {
                double age = getOriginAge();

                if ( age != ranges[i].birth )
                {
                    ranges[i].birth = age;
                    origin = age;
                    dirty_taxa[i] = true;
                }
            }
        }
    }

    // a species that ended below is named by the sampled ancestor here, which is a sample of it
    if ( pending_end > 0.0 && pending_species >= 0 )
    {
        if ( ranges[pending_species].death != pending_end )
        {
            ranges[pending_species].death = pending_end;
            dirty_taxa[pending_species] = true;
        }
        next_species[pending_species].ends_symmetric = true;
        pending_end = 0.0;
    }

    // no child carries this node's species, so it ends here and the node above has to name it
    if ( flow.species < 0 && sa_tip == false )
    {
        flow.end_age = node.getAge();
    }
    else
    {
        flow.end_age = pending_end;
    }

    return flow;
}

/**
 *
 *
 */
void FossilizedBirthDeathSpeciationProcess::prepareProbComputation( void ) const 
{
    AbstractFossilizedBirthDeathRangeProcess::prepareProbComputation();

    if ( homogeneous_lambda_a != NULL )
    {
        anagenetic = std::vector<double>(num_intervals, homogeneous_lambda_a->getValue() );
    }
    else
    {
        anagenetic = heterogeneous_lambda_a->getValue();
    }
    if ( homogeneous_beta != NULL )
    {
        symmetric = std::vector<double>(num_intervals, homogeneous_beta->getValue() );
    }
    else
    {
        symmetric = heterogeneous_beta->getValue();
    }

    for (size_t i = 0; i < num_intervals; i++)
    {
        if ( i < num_intervals-1 )
        {
            double dt = times[i+1] - times[i];

            q_tilde_i[i] = - anagenetic[i] * dt - symmetric[i] * (birth[i] + death[i] + fossil[i]) * dt + (1.0 - symmetric[i]) * q_tilde_i[i];
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
/** Re-take from each tip the ages it is authoritative for, before the recursion can bail. */
void FossilizedBirthDeathSpeciationProcess::repairRanges( void )
{
    // an extended tip is an extinction and may sit below its bins; a non-extended one is the
    // youngest appearance and must stay inside them
    for (size_t i = 0; i < getValue().getNumberOfNodes(); i++)
    {
        TopologyNode &node = getValue().getNode(i);
        // an extant tip is pinned at the present, which lies outside the fossil age range whenever
        // the taxon's youngest occurrence is old, so the range must not constrain it
        if ( node.isTip() )
        {
            size_t ti = node.getIndex();
            if ( ti >= taxa.size() ) continue;

            bool constrained = ( extended == false && taxa[ti].isExtinct() == true );
            node.setTipAgeUnconstrained( constrained == false );

            double age = node.getAge();

            // the recursion sets this too, but bails on an invalid continuation and leaves it stale
            if ( ranges[ti].death != age )
            {
                ranges[ti].death = age;
                dirty_taxa[ti] = true;
            }

            if ( constrained == true && ( ranges[ti].last != age || ( ranges[ti].singleton && ranges[ti].first != age ) ) )
            {
                ranges[ti].last = age;
                if ( ranges[ti].singleton == true ) ranges[ti].first = age;
                dirty_taxa[ti] = true;
            }
        }
    }
}


void FossilizedBirthDeathSpeciationProcess::updateRanges( void )
{
    // the recursion only sets this and bails once set, so a stale one would freeze the chain
    invalid_continuation = false;

    repairRanges();

    // the root lineage never reaches the assignment, so build a fresh table rather than edit one
    next_species.assign( taxa.size(), SpeciesEntry() );
    budding_lnProb = 0.0;

    const TopologyNode &root = getValue().getRoot();

    updateRanges(root);

    // partial_likelihood reads this table through rangeEndTerm, so an entry that moved is dirty
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( next_species[i] == species[i] ) continue;

        species[i] = next_species[i];
        dirty_taxa[i] = true;
    }

    // a lone lineage is both root and tip, so the recursion sets neither of its times
    if ( root.isTip() )
    {
        size_t i = root.getIndex();

        if ( getOriginAge() != ranges[i].birth )
        {
            ranges[i].birth = getOriginAge();
            dirty_taxa[i] = true;
        }
    }

    max_birth = 0;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( ranges[i].birth > ranges[max_birth].birth ) max_birth = i;
    }

    origin = getOriginAge();
}


void FossilizedBirthDeathSpeciationProcess::keepSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::restoreSpecialization(const DagNode *toucher)
{
    // or the next pull compares against the rejected tree and reads a moved entry as unchanged
    species = stored_species;

    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    bool first_touch = ( touched == false );

    if ( toucher == dag_node )
    {
        if ( first_touch )
        {
            stored_likelihood = partial_likelihood;
            stored_ranges     = ranges;
        }

        touched = true;
    }
    else
    {
        // every per-taxon term reads the rates and the timeline, so a parameter dirties all of them
        AbstractBirthDeathProcess::touchSpecialization(toucher, touchAll);
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }

    if ( first_touch == true ) stored_species = species;

    // the proposal has already written the value, so pull now; the pull reads the interval times
    prepareProbComputation();
    updateRanges();
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathSpeciationProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda_a)
    {
        heterogeneous_lambda_a = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_beta)
    {
        heterogeneous_beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda_a)
    {
        homogeneous_lambda_a = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_beta)
    {
        homogeneous_beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
        AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
    }
}

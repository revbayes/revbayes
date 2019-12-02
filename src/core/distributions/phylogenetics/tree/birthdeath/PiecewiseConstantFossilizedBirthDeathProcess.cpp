#include "DistributionExponential.h"
#include "PiecewiseConstantFossilizedBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

#include <algorithm>
#include <cmath>

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
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 */
PiecewiseConstantFossilizedBirthDeathProcess::PiecewiseConstantFossilizedBirthDeathProcess(const TypedDagNode<double> *ra,
                                                                                           const DagNode *inspeciation,
                                                                                           const DagNode *inextinction,
                                                                                           const DagNode *inpsi,
                                                                                           const DagNode *incounts,
                                                                                           const TypedDagNode<double> *inrho,
                                                                                           const TypedDagNode< RbVector<double> > *intimes,
                                                                                           const std::string &incondition,
                                                                                           const std::vector<Taxon> &intaxa,
                                                                                           bool uo,
                                                                                           bool pa,
                                                                                           bool ex) :
    AbstractBirthDeathProcess(ra, incondition, intaxa, uo),
    AbstractPiecewiseConstantFossilizedRangeProcess(inspeciation, inextinction, inpsi, incounts, inrho, intimes, intaxa, pa),
    extended(ex)
{
    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    redrawValue();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
PiecewiseConstantFossilizedBirthDeathProcess* PiecewiseConstantFossilizedBirthDeathProcess::clone( void ) const
{
    return new PiecewiseConstantFossilizedBirthDeathProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double PiecewiseConstantFossilizedBirthDeathProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    // variable declarations and initialization
    double lnProbTimes = computeLnProbabilityTimes();

    return lnProbTimes;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double PiecewiseConstantFossilizedBirthDeathProcess::computeLnProbabilityTimes( void ) const
{
    double lnProb = computeLnProbabilityRanges();

    // integrate speciation times for descendants of sampled ancestors
    for(size_t j = 0; j < I.size(); j++)
    {
        size_t i = I[j];

        double y_a = b_i[i];
        double o   = taxa[i].getAgeRange().getMax();

        size_t y_ai = l(y_a);

        size_t oi = presence_absence ? oldest_intervals[i] : l(o);

        // offset speciation density
        lnProb -= log( birth[y_ai] );

        // evaluate antiderivative at oi

        // replace q with q~ at the birth time
        double x = q(y_ai, y_a, true) - q(y_ai, y_a);

        // replace intermediate q terms
        for (size_t j = y_ai; j < oi; j++)
        {
            x -= q_i[j] - q_tilde_i[j];
        }

        if( presence_absence )
        {
            double a = std::max(d_i[i], times[oi]);
            double Ls_plus_a = oi > 0 ? std::min(y_a, times[oi-1]) : y_a;
            double Ls = Ls_plus_a - a;

            // replace H_i
            x += log(1.0 - exp(-Ls*fossil[oi]) ) - H[i];
        }
        else
        {
            // replace q terms at oi
            x += q(oi, o) - q(oi, o, true);
        }

        // compute definite integral
        lnProb += log(-expm1(x));
    }

    // if this is a sampled tree
    // replace extinction events with sampling events
    if( extended == false )
    {
        for(size_t i = 0; i < taxa.size(); i++)
        {
            size_t di = l(d_i[i]);

            // check constraints
            if( presence_absence )
            {
                if( d_i[i] > 0.0 )
                {
                    // yi == di
                    if( youngest_intervals[i] != di )
                    {
                        return RbConstants::Double::neginf;
                    }

                    // if the tip is a sampling event in the past
                    // then replace one unobserved fossil sample with an observed fossil sample
                    // i.e increment the observed fossil count
                    double Ls = times[di-1] - std::max(d_i[i], times[di]);
                    lnProb += fossil[di] - log( 1.0 - exp( - Ls * fossil[di] ) );
                }
            }
            // y == d
            else if ( fabs(d_i[i] - taxa[i].getAgeRange().getMin()) > 1E-5 )
            {
                return RbConstants::Double::neginf;
            }

            // if the tip is a sampling event in the past
            // replace observed extinction time with unobserved extinction time
            if ( d_i[i] > 0.0 )
            {
                lnProb -= death[di];
                lnProb += log( fossil[di] ) + log( p(d_i[i], di) );
            }
        }
    }

    // condition on survival
    if ( condition == "survival" )
    {
        double ps = log( pSurvival( getOriginAge(), 0) );

        lnProb -= use_origin ? ps : 2.0*ps;
    }

    return lnProb;
}


double PiecewiseConstantFossilizedBirthDeathProcess::getMaxTaxonAge( const TopologyNode& node ) const
{
    if( node.isTip() )
    {
        return node.getTaxon().getAgeRange().getMax();
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


double PiecewiseConstantFossilizedBirthDeathProcess::lnProbTreeShape(void) const
{
    // the fossilized birth death divergence times density is derived for an unlabeled oriented tree
    // so we convert to a labeled oriented tree probability by multiplying by 1 / n!
    // where n is the number of extant tips

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
double PiecewiseConstantFossilizedBirthDeathProcess::pSurvival(double start, double end) const
{
    double t = start;

    //std::vector<double> fossil_bak = fossil;

    //std::fill(fossil.begin(), fossil.end(), 0.0);

    double p0 = p(l(t), t);

    //fossil = fossil_bak;

    return 1.0 - p0;
}


/**
 *
 */
void PiecewiseConstantFossilizedBirthDeathProcess::simulateClade(std::vector<TopologyNode *> &n, double age, double present)
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the minimum birth age
    std::vector<double> first_occurrences;

    double current_age = RbConstants::Double::inf;
    double minimum_age = 0.0;

    for (size_t i = 0; i < n.size(); ++i)
    {
        // make sure the tip age is equal to the last occurrence
        if( n[i]->isTip() )
        {
            double min = n[i]->getTaxon().getAgeRange().getMin();
            double max = n[i]->getTaxon().getAgeRange().getMax();

            // in the extended tree, tip ages are extinction times
            if( extended )
            {
                n[i]->setAge( present + rng->uniform01() * ( min - present ) );
            }
            // in the sampled tree, tip ages are sampling times
            else
            {
                n[i]->setAge(min);
            }
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

    // reset the age
    double max_age = getOriginAge();

    if( minimum_age > max_age )
    {
        std::stringstream s;
        s << "Tree age is " << max_age << " but oldest fossil occurrence is " << minimum_age;
        throw(RbException(s.str()));
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
            }
            if( current_age >= first_occurrences[i] )
            {
                active_right_nodes.push_back( n[i] );
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

        // only simulate if there are at least two valid/active nodes and one active left node
        if ( active_nodes.size() <= 2 || active_right_nodes.empty() )
        {
            current_age = next_node_age;
        }
        else
        {
            // now we simulate new ages
            double next_sim_age = simulateNextAge(active_nodes.size()-2, age, present, current_age);

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

        if ( n.size() > 2 && current_age >= age  ) throw RbException("Unexpected number of taxa (remaining #taxa was " + StringUtilities::toString(n.size()) + " and age was " + current_age + " with maximum age of " + age + ") in tree simulation");

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
                throw(RbException("Cannot simulate clade of age " + StringUtilities::toString(age) + ", minimum age is " + StringUtilities::toString(minimum_age) ));
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
        parent->setAge( age );

        // insert the parent to our list
        n.push_back( parent );
    }
    else
    {
        throw RbException("Unexpected number of taxa (" + StringUtilities::toString(n.size()) + ") in tree simulation");
    }


}


/**
 * Simulate new speciation times.
 */
double PiecewiseConstantFossilizedBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant FBDP

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = l(present);

    // get the parameters
    double age = origin - present;
    double b = birth[i];
    double d = death[i];
    //double f = fossil[i];
    double r = homogeneous_rho->getValue();


    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    else
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }

    return present + t;
}


int PiecewiseConstantFossilizedBirthDeathProcess::updateStartEndTimes( const TopologyNode& node ) const
{
    if( node.isTip() )
    {
        return find(taxa.begin(), taxa.end(), node.getTaxon()) - taxa.begin();
    }

    int species = -1;

    std::vector<TopologyNode* > children = node.getChildren();

    bool sa = node.isSampledAncestor(true);

    for(int c = 0; c < children.size(); c++)
    {
        const TopologyNode& child = *children[c];

        int i = updateStartEndTimes(child);

        // if child is a tip, set the species/end time
        if( child.isTip() )
        {
            d_i[i] = child.getAge();
        }

        // is child a new species?
        // set start time at this node
        if( ( sa == false && c > 0 ) || ( sa && !child.isSampledAncestor() ) )
        {
            b_i[i] = node.getAge(); // y_{a(i)}

            if( sa ) I.push_back(i);
        }
        // child is the ancestral species
        else
        {
            // propagate species index
            species = i;

            // if this is the root
            // set the start time to the origin
            if( node.isRoot() )
            {
                b_i[i] = getOriginAge();
            }
        }
    }

    return species;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void PiecewiseConstantFossilizedBirthDeathProcess::updateStartEndTimes( void ) const
{
    I.clear();

    updateStartEndTimes(getValue().getRoot());
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void PiecewiseConstantFossilizedBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    AbstractPiecewiseConstantFossilizedRangeProcess::swapParameterInternal(oldP, newP);
}

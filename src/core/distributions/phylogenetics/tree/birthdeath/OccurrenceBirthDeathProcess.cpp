#include "OccurrenceBirthDeathProcess.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "RbMathMatrix.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

#include "ComputeLikelihoodsLtMt.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    sa             Start age of the process.
 * \param[in]    l              Speciation rate.
 * \param[in]    m              Extinction rate.
 * \param[in]    p              Fossil sampling rate.
 * \param[in]    o              Occurrence sampling rate.
 * \param[in]    rho            Sampling probability at present time.
 * \param[in]    r              Removal probability after sampling.
 * \param[in]    n              Algorithm accuracy (maximal number of hidden lineages).
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Extinct and extant taxa and their range ages.
 * \param[in]    tau            Times for which we want to compute the density.
 * \param[in]    uo             If true t is the origin time otherwise the root age of the process.
 * \param[in]    mt             If true computes densities with the Mt forward traversal algorithm otherwise uses Lt backward one.
 * \param[in]    tr             Initial tree (facultative).
 */
OccurrenceBirthDeathProcess::OccurrenceBirthDeathProcess(   const TypedDagNode<double> *sa,
                                                            const TypedDagNode<double> *l,
                                                            const TypedDagNode<double> *m,
                                                            const TypedDagNode<double> *p,
                                                            const TypedDagNode<double> *o,
                                                            const TypedDagNode<double> *rho,
                                                            const TypedDagNode<double> *r,
                                                            const TypedDagNode<long> *n,

                                                            const std::string& cdt,
                                                            const std::vector<Taxon> &tn,
                                                            const std::vector<double> &tau,
                                                            bool uo,
                                                            bool mt,

                                                            TypedDagNode<Tree> *tr) : AbstractBirthDeathProcess( sa, cdt, tn, uo ),
    start_age( sa ),
    lambda( l ),
    mu( m ),
    psi( p ),
    omega( o ),
    rho( rho ),
    removalPr( r ),
    maxHiddenLin( n ),
    cond (cdt),
    taxa (tn),
    time_points( tau ),
    useOrigin (uo),
    useMt ( mt )

{
    addParameter( start_age );
    addParameter( lambda );
    addParameter( mu );
    addParameter( psi );
    addParameter( omega );
    addParameter( rho );
    addParameter( removalPr );
    addParameter( maxHiddenLin );


    if (tr != NULL)
    {
      delete value;
      value = &(tr->getValue());
    }
    else
    {
      simulateTree();
    }

}


OccurrenceBirthDeathProcess::~OccurrenceBirthDeathProcess( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
OccurrenceBirthDeathProcess* OccurrenceBirthDeathProcess::clone( void ) const
{

    return new OccurrenceBirthDeathProcess( *this );
}



/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return    The log-probability density.
 */
double OccurrenceBirthDeathProcess::computeLnProbabilityDivergenceTimes( void ) const
{

    // prepare the probability computation
    prepareProbComputation();

    // compute the log-likelihood : use ComputeLikelihoodsBackwardsLt (backward traversal of the tree) or ComputeLikelihoodsForwardsMt (forward traversal of the tree)
    const RevBayesCore::Tree tree(*value);

    if (useMt) {
        const std::vector<double> time_points_Mt( 1, start_age->getValue() );
        MatrixReal B_Mt = RevBayesCore::ComputeLikelihoodsForwardsMt(start_age, lambda, mu, psi, omega, rho, removalPr, maxHiddenLin, cond, time_points_Mt, useOrigin, tree);
        
        const double rh = rho->getValue();
        const long N = maxHiddenLin->getValue();
        const size_t k = value->getNumberOfExtantTips();
        
        double likelihood = B_Mt[0][0];
        for(int i = 1; i < N+1; i++){
            likelihood += B_Mt[0][i] * pow(rh,k) * pow(1.0 - rh,i);
        return likelihood;
        }
    }
    const std::vector<double> time_points_Lt(1, 0.0);
    MatrixReal B_Lt = RevBayesCore::ComputeLikelihoodsBackwardsLt(start_age, lambda, mu, psi, omega, rho, removalPr, maxHiddenLin, cond, time_points_Lt, useOrigin, tree);
    
    double likelihood = B_Lt[0][0];
    return likelihood;
}


double OccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
{

    double lnProbTimes = 0.0;
    double process_time = getOriginAge();
    size_t num_initial_lineages = 2;
    const TopologyNode& root = value->getRoot();

    if (use_origin) {
        // If we are conditioning on survival from the origin,
        // then we must divide by 2 the log survival probability computed by AbstractBirthDeathProcess
        num_initial_lineages = 1;
    }

    // if conditioning on root, root node must be a "true" bifurcation event
    else if (root.getChild(0).isSampledAncestor() || root.getChild(1).isSampledAncestor())
    {
        return RbConstants::Double::neginf;
    }

    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double death_rate = mu->getValue();
    double serial_rate = psi->getValue();
    double sampling_prob = rho->getValue();

    // get helper variables
    double a = birth_rate - death_rate - serial_rate;
    double c1 = std::fabs(sqrt(a * a + 4 * birth_rate * serial_rate));
    double c2 = -(a - 2 * birth_rate * sampling_prob) / c1;

    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();

    // classify nodes
    int num_sampled_ancestors = 0;
    int num_extant_taxa = 0;

    std::vector<double> serial_tip_ages = std::vector<double>();
    std::vector<double> internal_node_ages = std::vector<double>();
    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = value->getNode( i );

        if ( n.isFossil() && n.isSampledAncestor() )
        {
            // node is sampled ancestor
            num_sampled_ancestors++;
        }
        else if ( n.isFossil() && !n.isSampledAncestor() )
        {
            // node is serial leaf
            serial_tip_ages.push_back( n.getAge() );
        }
        else if ( n.isTip() && !n.isFossil() )
        {
            // node is extant leaf
            num_extant_taxa++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            if (!n.isRoot() || use_origin)
            {
                // node is bifurcation event (a "true" node)
                internal_node_ages.push_back( n.getAge() );
            }
        }
    }

    // add the log probability for the serial sampling events
    if (serial_rate == 0.0)
    {
        if ( serial_tip_ages.size() + num_sampled_ancestors > 0 )
        {
            return RbConstants::Double::neginf;
            //throw RbException("The serial sampling rate is zero, but the tree has serial sampled tips.");
        }
    }
    else
    {
        lnProbTimes += (serial_tip_ages.size() + num_sampled_ancestors) * log( serial_rate );
    }

    // add the log probability for sampling the extant taxa
    if (num_extant_taxa > 0)
    {
        lnProbTimes += num_extant_taxa * log( 4.0 * sampling_prob );
    }

    // add the log probability of the initial sequences
    lnProbTimes += -lnQ(process_time, c1, c2) * num_initial_lineages;

    // add the log probability for the internal node ages
    lnProbTimes += internal_node_ages.size() * log( birth_rate );
    for (size_t i=0; i<internal_node_ages.size(); i++)
    {
        lnProbTimes -= lnQ(internal_node_ages[i], c1, c2);
    }

    // add the log probability for the serial tip ages
    for (size_t i=0; i < serial_tip_ages.size(); i++)
    {
        double t = serial_tip_ages[i];
        lnProbTimes += log(pZero(t, c1, c2)) + lnQ(t, c1, c2);
    }

    // condition on survival
    if ( condition == "survival")
    {
        lnProbTimes -= num_initial_lineages * log(1.0 - pHatZero(process_time));
    }
    // condition on nTaxa
    else if ( condition == "nTaxa" )
    {
        lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
    }

    std::cout << "Compute lnProbTimes output : " << lnProbTimes << std::endl;
    return lnProbTimes;

}


double OccurrenceBirthDeathProcess::functionP(double t, double z) const
{
    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();

    double gamma = birth + death + ps + om;
    double sqrtDelta = sqrt( pow(gamma, 2) - 4.0 * birth * death );
    double x1 = (gamma - sqrtDelta)/(2*birth);
    double x2 = (gamma + sqrtDelta)/(2*birth);
    
    return pow((sqrtDelta/birth),2) * pow((1/((x2-z) - (x1-z)*exp(-sqrtDelta*t))),2) * exp(-sqrtDelta*t) * (1-z);
}

double OccurrenceBirthDeathProcess::functionU(double t, double z) const
{
    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();

    double gamma = birth + death + ps + om;
    double sqrtDelta = sqrt( pow(gamma, 2) - 4.0 * birth * death );
    double x1 = (gamma - sqrtDelta)/(2*birth);
    double x2 = (gamma + sqrtDelta)/(2*birth);
    
    double numerator = x1*(x2-z) - x2*(x1-z)*exp(-sqrtDelta*t);
    double denominator = (x2-z) - (x1-z)*exp(-sqrtDelta*t);
    return numerator/denominator;
}


double OccurrenceBirthDeathProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / n!
    // where n is the number of extant tips, m is the number of sampled extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    int num_extinct = (int)value->getNumberOfExtinctTips();
    int num_sa = (int)value->getNumberOfSampledAncestors();

    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa - num_extinct);
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Speciation rate at time t.
 */
 double OccurrenceBirthDeathProcess::pSurvival(double start, double end) const
 {
        return 1.0 - pHatZero(end);
 }



/**
 * Simulate new speciation times.
 */
double OccurrenceBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant FBDP

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the parameters
    double age = origin - present;
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho->getValue();

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


double OccurrenceBirthDeathProcess::pZero(double t, double c1, double c2) const
{
//work in progress function u and p
    // get helper variables
    const double birth = lambda->getValue();

    const double death = mu->getValue();

    const double ps = psi->getValue();

    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();

    const double gamma = birth + death + ps + om;

    // double a = birth - death - ps;
    // double c1 = std::fabs(sqrt(a * a + 4 * birth * ps));
    // double c2 = -(a - 2 * birth * rh) / c1;
    // // roots of the polynmial ODE
    // double Delta = pow(gamma,2) - 4*birth*death;
    // double x1 = (gamma - sqrt(Delta))/(2*birth);
    // double x2 = (gamma + sqrt(Delta))/(2*birth);

    double b = lambda->getValue();
    double d = mu->getValue();
    double f = psi->getValue();
    double v1 = exp(-c1*t) * (1.0-c2) - (1.0+c2);
    double v2 = exp(-c1*t) * (1.0-c2) + (1.0+c2);
    double v = (b + d + f + c1 * (v1 / v2)) / (2.0 * b);
    return v;
}



double OccurrenceBirthDeathProcess::lnQ(double t, double c1, double c2) const
{
    //return log( 2*(1-c2*c2) + exp(-c1*t)*(1-c2)*(1-c2) + exp(c1*t)*(1+c2)*(1+c2) );

    // numerically safe code
    return c1*t + 2 * log( exp(-c1*t) * (1-c2) + (1+c2));
}


double OccurrenceBirthDeathProcess::pHatZero(double t) const
 {
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho->getValue();
    double val = r * (b-d) / (b*r + (b*(1-r)-d)*exp(-(b-d)*t));
    return 1.0 - val;
}




/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void OccurrenceBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == psi)
    {
        psi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == omega)
    {
        omega = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == removalPr)
    {
        removalPr = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }

}

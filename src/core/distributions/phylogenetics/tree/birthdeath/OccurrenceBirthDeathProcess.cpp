#include "OccurrenceBirthDeathProcess.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    t              Time of the origin/present/length of the process.
 * \param[in]    l              Speciation rate.
 * \param[in]    m              Extinction rate.
 * \param[in]    p              Fossil sampling rate.
 * \param[in]    o              Occurrence sampling rate.
 * \param[in]    rho            Sampling probability at present time.
 * \param[in]    r              Removal probability after sampling.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tau            Times for which we want to compute the density.
 * \param[in]    uo             If true t is the origin time otherwise the root age of the process.
 * \param[in]    tr             Initial tree (facultative).
 */
OccurrenceBirthDeathProcess::OccurrenceBirthDeathProcess( const TypedDagNode<double> *t,
                                                          const TypedDagNode<double> *l,
                                                          const TypedDagNode<double> *m,
                                                          const TypedDagNode<double> *p,
                                                          const TypedDagNode<double> *o,
                                                          const TypedDagNode<double> *rho,
                                                          const TypedDagNode<double> *r,

                                                          const std::string& cdt,
                                                          const std::vector<Taxon> &tn,
                                                          const TypedDagNode< RbVector<double> > *tau,
                                                          bool uo,

                                                          TypedDagNode<Tree> *tr) : AbstractBirthDeathProcess( t, cdt, tn, uo ),
    tor( t ),
    lambda( l ),
    mu( m ),
    psi( p ),
    omega( o ),
    rho( rho ),
    removalPr( r ),
    dn_time_points ( tau )

{
    addParameter( tor );
    addParameter( lambda );
    addParameter( mu );
    addParameter( psi );
    addParameter( omega );
    addParameter( rho );
    addParameter( removalPr );


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

    // step 1: Create the list of all event times
    poolTimes();
   // std::cout << "Printing the list of events below" << std::endl;
   // for (int i =0; i < events.size(); i++){
   //     std::cout << events[i].time << " : " << events[i].type << std::endl;
   // }

    //std::cout << "Parameter values are the following : lambda, mu, psi, rho = " << lambda -> getValue() << " , " << mu -> getValue() << " , " << psi -> getValue() << " , " << rho -> getValue() << std::endl;

    // compute the log-likelihood : use ComputeLt (backward traversal of the tree) or ComputeMt (forward traversal of the tree)
    // double lnProbTimes_Mt = ComputeMt();
    // double lnProbTimes_Lt = ComputeLt();

    double lnProbTimes = computeLnProbabilityTimes();
    // std::cout << "Difference : "<< lnProbTimes-lnProbTimes_Mt << std::endl;

    // double lnProbTimes2 = computeLnProbabilityTimes2();

    return lnProbTimes;
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

double OccurrenceBirthDeathProcess::computeLnProbabilityTimes2( void ) const
{
    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double death_rate = mu->getValue();
    double serial_rate = psi->getValue();
    double sampling_prob = rho->getValue();

    // get helper variables
    double a = birth_rate - death_rate - serial_rate;
    double c1 = std::fabs(sqrt(a * a + 4 * birth_rate * serial_rate));
    double c2 = -(a - 2 * birth_rate * sampling_prob) / c1;
    
    double lnProbTimes = 0;
    for(int h = 0; h < events.size(); h++){
        // deal with t > tor

        double th = events[h].time;
        std::string type = events[h].type;

        if(type == "sampled ancestor"){
            lnProbTimes += log( serial_rate );
        }
        else if( type == "terminal non-removed"){
            lnProbTimes += log( functionU(th, 1-sampling_prob) ) - log( functionP(th, 1-sampling_prob) ) + log( serial_rate );
        }
        else if( type == "branching time"){
            lnProbTimes += log( functionP(th, 1.0-sampling_prob)) + log( birth_rate );
            //lnProbTimes += log( 4.0 * sampling_prob ) - lnQ(th, c1, c2) + log( birth_rate );
        }
        else if( type == "origin"){
            lnProbTimes +=  log( functionP(th, 1.0 - sampling_prob) );
        }
    }

    std::cout << "Compute lnProbTimes2 output : " << lnProbTimes << std::endl;
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

/**
 * Construct the vector containig all branching and sampling times + time points for which we want to compute the density.
 */
void OccurrenceBirthDeathProcess::poolTimes( void ) const
{
    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();

    extant = 0;

    // classify nodes
    events.clear();
    events.push_back(Event(tor->getValue(), "origin"));

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = value->getNode( i );
        
        //isFossil is an optional condition to obtain sampled ancestor node ages
        /*
        Node labels :
        fl = fossil leaf
        b  = "true" bifurcation
        b' = "false" bifurcation (artefact of the sampled ancestors representation)
        sa = sampled ancestor
        el = extant leaf

         __|___
        |  b   |
        |      |
        fl     |
             b'|___ sa
               |
               |
               el

         1. Pick a fossil among those with brl > 0 (prob = 1/m)
         2. Set brl = 0
         */

        if ( n.isFossil() && n.isSampledAncestor() )  //isFossil is optional (all sampled ancestors are fossils)
        {
            // node is sampled ancestor
            events.push_back(Event(n.getAge(), "sampled ancestor")) ;

        }
        else if ( n.isFossil() && !n.isSampledAncestor() )
        {
            // node is fossil leaf (named terminal non-removed in Lt)
            // events.push_back(Event(n.getAge(),"fossil leaf") ;
            events.push_back(Event(n.getAge(),"terminal non-removed")) ;
        }
        else if ( !n.isFossil() && n.isTip() )
        {
            // node is extant leaf : only their number is necessary to compute Lt and Mt
            // events.push_back(Event(n.getAge(),"extant leaf") ;
            // events.push_back(Event(n.getAge(), "terminal non-removed")) ;
            //std::cout << n.getSpeciesName() << std::endl;
            extant++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            // std::cout << "Is branching node root ? " << n.isRoot() << std::endl;
            
            // node is a "true" bifurcation event
            events.push_back(Event(n.getAge(),"branching time")) ;
        }
        else
        {
            std::cout << "Warning : non-categorized node" << std::endl;
        }   

    }

    events.push_back(Event(0.0,"present time")) ;
}

  //Pool all observations together into

    // // add the log probability for the serial sampling events
    // if (ps == 0.0)
    // {
    //     if ( serial_tip_ages.size() + sampled_ancestors_ages.size() > 0 )
    //     {
    //         return RbConstants::Double::neginf;
    //         //throw RbException("The serial sampling rate is zero, but the tree has serial sampled tips.");
    //     }
    // }
    // else
    // {
    //     lnProbTimes += (serial_tip_ages.size() + sampled_ancestors_ages.size() ) * log( ps );
    // }
    //
    // // add the log probability for sampling the extant taxa
    // if (num_extant_taxa > 0)
    // {
    //     lnProbTimes += num_extant_taxa * log( 4.0 * rh );
    // }
    //
    // // add the log probability of the initial sequences
    // lnProbTimes += -lnQ(process_time, c1, c2) * num_initial_lineages;
    //
    // // add the log probability for the internal node ages
    // lnProbTimes += internal_node_ages.size() * log( birth );
    // for (size_t i=0; i<internal_node_ages.size(); i++)
    // {
    //     lnProbTimes -= lnQ(internal_node_ages[i], c1, c2);
    // }
    //
    // // add the log probability for the serial tip ages
    // for (size_t i=0; i < serial_tip_ages.size(); i++)
    // {
    //     double t = serial_tip_ages[i];
    //     lnProbTimes += log(pZero(t, c1, c2)) + lnQ(t, c1, c2);
    // }
    //
    // // condition on survival
    // if ( condition == "survival")
    // {
    //     lnProbTimes -= num_initial_lineages * log(1.0 - pHatZero(process_time));
    // }
    // // condition on nTaxa
    // else if ( condition == "nTaxa" )
    // {
    //     lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
    // }



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

/**
 * Get the log-transformed probability of the current value under the current parameter values : call Lt or Mt.
 *
 * \return    The log-probability density.
 */
// double OccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
//  {
//  double lnProbTimes = ComputeMt();
//  return lnProbTimes;
//  }

//from Warnock & Manceau computeLt function
 /**
 * Compute the log-transformed probability of the current value under the current parameter values : breadth-first forward traversal algorithm.
 *
 * \return    The log-probability density.
 */
double OccurrenceBirthDeathProcess::ComputeMt( void ) const
{
    // order times oldest to youngest
    std::sort( events.begin(), events.end(), AgeCompareReverse() );

    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();
    const double gamma = birth + death + ps + om;

    size_t N = 20; // accuracy of the algorithm
    size_t S = dn_time_points->getValue().size();

    // Initialize an empty matrix and a cursor to write lines in this matrix
    MatrixReal B(S, (N + 1), 0.0);
    size_t indxJ = S-1;

    // We start at the time of origin, supposedly the first time in the vector of events
    size_t k = 1;
    RbVector<double> Mt(N+1, 0.0);
    Mt[0] = 1;
    double thPlusOne = events[0].time;

    // Then we iterate over the next events
    for(int h = 1; h < events.size(); h++){

        // First, deal with the update on time period (th, thPlusOne)
        double th = events[h].time;

        if(th > tor->getValue()) {
            std::cout << "ERROR : th > tor : " << th << " > " << tor->getValue() << std::endl;
            continue;
        };

        if( th != thPlusOne ){
            MatrixReal A( (N+1), (N+1), 0.0 );
            
            for(int i = 0; i < (N + 1); i++){
                A[i][i] = gamma * (k + i) * (th-thPlusOne);
                if (i < N) A[i][i+1] = -death * (i + 1) * (th-thPlusOne);
                if (i > 0) A[i][i-1] = -birth * (2*k + i - 1) * (th-thPlusOne);
            }
            
            RbMath::expMatrixPade(A, A, 4);
            Mt = A * Mt;
        }

        // Second, deal with the update at punctual event th
        std::string type = events[h].type;

        if(type == "time slice"){
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = Mt[i];
            }
            indxJ -= 1;
        }

        if(type == "terminal removed"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= ps * rp;
            }
            k -= 1;
        }

        if(type == "terminal non-removed"){
            for(int i = N; i > 0; i--){
                Mt[i] = Mt[i-1] * ps * (1-rp);
            }
            Mt[0] = 0;
            k -= 1;
        }

        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= ps * (1-rp);
            }
        }

        if(type == "occurrence removed"){
            for(int i = 0; i < N; i++){
                Mt[i] = Mt[i+1] * (i+1) * om * rp;
            }
         }

        if(type == "occurrence non-removed"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= (k+i) * om * (1-rp);
            }
        }

        if(type == "branching time"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= birth;
            }
            k += 1;
        }

        thPlusOne = th;
    }

    double likelihood = Mt[0];
    for(int i = 1; i < N+1; i++){
        likelihood += Mt[i] * pow(rh,k) * pow(1.0 - rh,i);
    }

    std::cout << "Compute Mt output : " << log(likelihood) << std::endl;
    return log(likelihood);
}


//from Warnock & Manceau computeLt function
/**
* Compute the log-transformed probability of the current value under the current parameter values : breadth-first backward traversal algorithm.
*
* \return    The log-probability density.
*/
double OccurrenceBirthDeathProcess::ComputeLt( void ) const
{
    // order times youngest to oldest
    std::sort( events.begin(), events.end(), AgeCompare() );

    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();
    const double gamma = birth + death + ps + om;

    size_t N = 20; // accuracy of the algorithm

    // Initialize an empty matrix and a cursor indxJ to write lines in this matrix
    size_t S = dn_time_points->getValue().size();
    MatrixReal B(S, (N + 1), 0.0);
    size_t indxJ = 0;

    // We start at time 0 with type "present" in the vector of events
    size_t k = extant;
    double thMinusOne = events[0].time;
    RbVector<double> Lt(N+1, pow(rh, k));
    for ( int i = 0; i < (N + 1); i++){
        Lt[i] *= pow( 1.0-rh, i );
    }

    // We then iterate over the following events until finding the time of origin
    for(int h = 1; h < events.size(); h++){

        // First deal with the update along (thMinusOne, th)
        double th = events[h].time;

        if( th != thMinusOne ){

            MatrixReal A( (N+1), (N+1), 0.0 );
            for(int i = 0; i < (N + 1); i++){
              A[i][i] = -gamma * (k + i) * (th - thMinusOne);
              if (i < N) A[i][i+1] = birth * ( (2 * k) + i ) * (th - thMinusOne);
              if (i > 0) A[i][i-1] = death * i * (th - thMinusOne);
            }
            RbMath::expMatrixPade(A, A, 4);
            Lt = A * Lt;
        }

        // Second, deal with the update at time th
        std::string type = events[h].type;

        if(type == "time slice"){
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = Lt[i];
            }
            indxJ += 1;
        }

        if(type == "terminal removed"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * ps * rp;
            }
            k += 1;
        }

        if(type == "terminal non-removed"){
            for(int i = 0; i < N; i++){
                Lt[i] = Lt[i+1] * ps * (1.0-rp);
            }
            k += 1;
        }

        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * ps * (1.0-rp);
            }
        }

        if(type == "occurrence removed"){
            for(int i = N; i > 0; i--){
                Lt[i] = Lt[i-1] * i * om * rp;
            }
            Lt[0] = 0;
         }

        if(type == "occurrence non-removed"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * (k+i) * om * (1-rp);
            }
        }

        if(type == "branching time"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * birth;
            }
            k -= 1;
        }
        
        thMinusOne = th;
    }

    std::cout << "Compute Lt output : " << log(Lt[0]) << std::endl;
    return log(Lt[0]);
}

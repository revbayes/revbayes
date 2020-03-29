#include "ComputeLikelihoodsLtMt.h"

#include <vector>
#include <iostream>
#include <cmath>

#include "RbMathMatrix.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "Tree.h"
#include "TopologyNode.h"

namespace RevBayesCore {
    class DagNode;
    class MatrixReal; }

using namespace RevBayesCore;

/**
 * Compute the joint probability density of the observations made up to any time t and the population
 * size at that time, as time decreases towards present : breadth-first forward traversal algorithm.
 *
 * \param[in]    start_age              Start age of the process.
 * \param[in]    lambda                 Speciation rate.
 * \param[in]    mu                     Extinction rate.
 * \param[in]    psi                    Extinction sampling rate.
 * \param[in]    omega                  Occurrence sampling rate.
 * \param[in]    rho                    Sampling probability at present time.
 * \param[in]    removalPr              Removal probability after sampling.
 * \param[in]    maxHiddenLin           Algorithm accuracy (maximal number of hidden lineages).
 * \param[in]    cond                   Condition of the process (none/survival/#Taxa).
 * \param[in]    time_points         Times for which we want to compute the density.
 * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
 * \param[in]    timeTree               Tree for ancestral populations size inference.
 *
 * \return    The matrix of Mt values through time.
*/
MatrixReal RevBayesCore::ComputeLikelihoodsForwardsMt(    const TypedDagNode<double> *start_age,
                                                          const TypedDagNode<double> *lambda,
                                                          const TypedDagNode<double> *mu,
                                                          const TypedDagNode<double> *psi,
                                                          const TypedDagNode<double> *omega,
                                                          const TypedDagNode<double> *rho,
                                                          const TypedDagNode<double> *removalPr,
                                                          const TypedDagNode<long> *maxHiddenLin,

                                                          const std::string &cond,
                                                          const std::vector<double> &time_points,
                                                          bool useOrigin,
                                                          const std::vector<double> &occurrence_ages,
                                                          const Tree &timeTree)
{

    // Construct the vector containig all branching and sampling times + time points for which we want to compute the density.
    struct Event {
            Event(double d, std::string s) : time(d), type(s) {};

            std::string type;
            double time;

            std::string getEventType(void){ return type; }
            double getEventTime(void){ return time; }

        };

    // get node/time variables
    const size_t num_nodes = timeTree.getNumberOfNodes();

    // number of extant lineages
    size_t k = 0;

    // vector of events
    std::vector<Event>         events;

    // classify nodes
    events.clear();
    events.push_back(Event(start_age->getValue(), "origin"));

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = timeTree.getNode( i );

        //isFossil is an optional condition to obtain sampled ancestor node ages
        /*
        Node labels :
        fl = fossil leaf
        b  = "true" bifurcation
        b' = "false" bifurcation (artefact of the sampled ancestors representation)
        sa = sampled ancestor
        el = extant leaf
        . = time points at which density is computed

         __|___             .
        |  b   |
        |      |            .
        fl     |
             b'|___ sa      .
               |
               |            .
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
            events.push_back(Event(n.getAge(),"fossil leaf")) ;
            // events.push_back(Event(n.getAge(),"terminal non-removed")) ;
        }
        else if ( !n.isFossil() && n.isTip() )
        {
            // node is extant leaf : only their number is necessary to compute Lt and Mt
            // events.push_back(Event(n.getAge(),"extant leaf") ;
            // events.push_back(Event(n.getAge(), "terminal non-removed")) ;
            // std::cout << n.getSpeciesName() << std::endl;
            k++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            // node is a "true" bifurcation event
            events.push_back(Event(n.getAge(),"branching time")) ;
        }
    }

    const RbVector<double> tau = time_points;

    for (size_t i = 0; i < tau.size(); i++)
    {
        events.push_back(Event(tau[i],"time point")) ;
    }

    for ( int i=0; i < occurrence_ages.size(); ++i)
    {
        events.push_back(Event(occurrence_ages[i],"occurrence non-removed")) ;
    }

    events.push_back(Event(0.0,"present time")) ;

    // mutable std::vector<Event> events = RevBayesCore::poolTimes(*start_age, *time_points, timeTree);

    // order times oldest to youngest
    struct AgeCompareReverse {
        bool operator()(const Event first, const Event second) {
            return first.time > second.time;
        }
    };
    std::sort( events.begin(), events.end(), AgeCompareReverse() );



    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();
    const long N = maxHiddenLin->getValue();
    // const RbVector<double> tau = time_points->getValue();

    const size_t S = tau.size();
    const double gamma = birth + death + ps + om;

    // Initialize an empty matrix and a cursor to write lines in this matrix
    MatrixReal B(S, (N + 1), 0.0);
    size_t indxJ = S-1;

    // Count the number of extant lineages
    // mutable size_t k = 0;
    // for(int h = 1; h < events.size(); h++){
    //     if(type == "extant leaf"){
    //         k += 1;
    //     }
    // }

    // Number of extant lineages
    k = 1;

    // We start at the time of origin, supposedly the first time in the vector of events
    RbVector<double> Mt(N+1, 0.0);
    Mt[0] = 1;
    double thPlusOne = events[0].time;
    std::cout << "Event time : " << events[0].time << " - Event type : " << events[0].type << std::endl;
    std::cout << k << std::endl;
    
    if(thPlusOne != start_age->getValue()) {
        std::cout << "WARNING : thPlusOne != start_age : " << thPlusOne << " != " << start_age->getValue() << std::endl;
    };

    // Then we iterate over the next events
    for(int h = 1; h < events.size(); h++){

        // First, deal with the update on time period (th, thPlusOne)
        double th = events[h].time;

        if(th > start_age->getValue()) {
            std::cout << "WARNING : th > start_age : " << th << " > " << start_age->getValue() << std::endl;
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

        if(type == "time point"){
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

        if(type == "fossil leaf"){
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

        std::cout << "Event time : " << th << " - Event type : " << type << " -> Mt[0] : " << Mt[0] << " / Mt[1] : " << Mt[1] << std::endl;

        thPlusOne = th;
    }

    // double likelihood = Mt[0];
    // for(int i = 1; i < N+1; i++){
    //     likelihood += Mt[i] * pow(rh,k) * pow(1.0 - rh,i);
    // }

    return B;
}

/**
 * Compute the probability density of observations made down any time t, conditioned on the population
 * size at that time, as time increases towards the past : breadth-first backward traversal algorithm.
 *
 * \param[in]    start_age              Start age of the process.
 * \param[in]    lambda                 Speciation rate.
 * \param[in]    mu                     Extinction rate.
 * \param[in]    psi                    Extinction sampling rate.
 * \param[in]    omega                  Occurrence sampling rate.
 * \param[in]    rho                    Sampling probability at present time.
 * \param[in]    removalPr              Removal probability after sampling.
 * \param[in]    maxHiddenLin           Algorithm accuracy (maximal number of hidden lineages).
 * \param[in]    cond                   Condition of the process (none/survival/#Taxa).
 * \param[in]    time_points         Times for which we want to compute the density.
 * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
 * \param[in]    timeTree               Tree for ancestral populations size inference.
 *
 * \return    The matrix of Lt values through time.
 */
MatrixReal RevBayesCore::ComputeLikelihoodsBackwardsLt(   const TypedDagNode<double> *start_age,
                                                          const TypedDagNode<double> *lambda,
                                                          const TypedDagNode<double> *mu,
                                                          const TypedDagNode<double> *psi,
                                                          const TypedDagNode<double> *omega,
                                                          const TypedDagNode<double> *rho,
                                                          const TypedDagNode<double> *removalPr,
                                                          const TypedDagNode<long> *maxHiddenLin,

                                                          const std::string& cond,
                                                          const std::vector<double> &time_points,
                                                          bool useOrigin,
                                                          const std::vector<double> &occurrence_ages,
                                                          const Tree &timeTree)
{
    // Construct the vector containig all branching and sampling times + time points for which we want to compute the density.
    struct Event {
            Event(double d, std::string s) : time(d), type(s) {};

            std::string type;
            double time;

            std::string getEventType(void){ return type; }
            double getEventTime(void){ return time; }

        };

    // get node/time variables
    const size_t num_nodes = timeTree.getNumberOfNodes();

    // number of extant lineages
    size_t k = 0;

    // vector of events
    std::vector<Event>         events;

    // classify nodes
    events.clear();
    events.push_back(Event(start_age->getValue(), "origin"));

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = timeTree.getNode( i );

        //isFossil is an optional condition to obtain sampled ancestor node ages
        /*
        Node labels :
        fl = fossil leaf
        b  = "true" bifurcation
        b' = "false" bifurcation (artefact of the sampled ancestors representation)
        sa = sampled ancestor
        el = extant leaf
        . = time points at which density is computed

         __|___             .
        |  b   |
        |      |            .
        fl     |
             b'|___ sa      .
               |
               |            .
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
            events.push_back(Event(n.getAge(),"fossil leaf")) ;
            // events.push_back(Event(n.getAge(),"terminal non-removed")) ;
        }
        else if ( !n.isFossil() && n.isTip() )
        {
            // node is extant leaf : only their number is necessary to compute Lt and Mt
            // events.push_back(Event(n.getAge(),"extant leaf") ;
            // events.push_back(Event(n.getAge(), "terminal non-removed")) ;
            // std::cout << n.getSpeciesName() << std::endl;
            k++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            // node is a "true" bifurcation event
            events.push_back(Event(n.getAge(),"branching time")) ;
        }
        else
        {
            std::cout << "Warning : non-categorized node" << std::endl;
        }
    }

    const RbVector<double> tau = time_points;

    for (size_t i = 0; i < tau.size(); i++)
    {
        events.push_back(Event(tau[i],"time point")) ;
    }

    for ( int i=0; i < occurrence_ages.size(); ++i)
    {
    events.push_back(Event(occurrence_ages[i],"occurrence non-removed")) ;
    }
    
    events.push_back(Event(0.0,"present time")) ;

    // mutable std::vector<Event> events = RevBayesCore::poolTimes(*start_age, *time_points, timeTree);

    // order times youngest to oldest
    struct AgeCompare {
        bool operator()(const Event first, const Event second) {
            return first.time < second.time;
        }
    };
    std::sort( events.begin(), events.end(), AgeCompare() );





    const double birth = lambda->getValue();
    const double death = mu->getValue();
    const double ps = psi->getValue();
    const double om = omega->getValue();
    const double rh = rho->getValue();
    const double rp = removalPr->getValue();
    const long N = maxHiddenLin->getValue();
    // const RbVector<double> tau = time_points->getValue();

    const size_t S = tau.size();
    const double gamma = birth + death + ps + om;

    // Initialize an empty matrix and a cursor indxJ to write lines in this matrix
    MatrixReal B(S, (N + 1), 0.0);
    size_t indxJ = 0;

    // We start at time 0 with type "present" in the vector of events
    // size_t k = extant;
    double thMinusOne = events[0].time;
    RbVector<double> Lt(N+1, pow(rh, k));
    for ( int i = 0; i < (N + 1); i++){
        Lt[i] *= pow( 1.0-rh, i );
    }

    std::cout << "Event time : " << events[0].time << " - Event type : " << events[0].type << std::endl;

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

        if(type == "time point"){
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

        if(type == "fossil leaf"){
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

        std::cout << "Event time : " << th << " - Event type : " << type << " -> Lt[0] : " << Lt[0] << " / Lt[1] : " << Lt[1] << std::endl;

        thMinusOne = th;
    }

    return B;
}

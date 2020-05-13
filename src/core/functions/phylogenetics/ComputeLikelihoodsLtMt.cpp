#include "ComputeLikelihoodsLtMt.h"

#include <vector>
#include <iostream>
#include <cmath>

#include "RbConstants.h"
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
 * Compute the joint log-probability density of the observations made up to any time t and the population
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
 * \param[in]    time_points            Times for which we want to compute the density.
 * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
 * \param[in]    returnLogLikelihood    If true the function returns the log likelihood instead of the full B_Mt matrix.
 * \param[in]    verbose                If true displays warnings and information messages.
 * \param[in]    occurrence_ages        Vector of occurrence ages.
 * \param[in]    timeTree               Tree for ancestral populations size inference.
 *
 * \return    The matrix of log-Mt values through time.
*/
MatrixReal RevBayesCore::ComputeLnProbabilitiesForwardsMt(  const TypedDagNode<double> *start_age,
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
                                                            bool returnLogLikelihood,
                                                            bool verbose,
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
    size_t nb_extant = 0;

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
            // node is fossil leaf
            // For now, we assume there is no information/labels on removal
            // @todo add fossil-removed/non-removed
            events.push_back(Event(n.getAge(),"fossil leaf")) ;
            // events.push_back(Event(n.getAge(),"terminal non-removed")) ;
        }
        else if ( !n.isFossil() && n.isTip() )
        {
            // node is extant leaf : only their number is necessary to compute Lt and Mt
            // events.push_back(Event(n.getAge(),"extant leaf") ;
            // events.push_back(Event(n.getAge(), "terminal non-removed")) ;
            // if (verbose){std::cout << n.getSpeciesName() << std::endl;}
            nb_extant++;
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
        events.push_back(Event(occurrence_ages[i],"occurrence")) ;
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
    MatrixReal B(S, (N + 1), RbConstants::Double::neginf);
    size_t indxJ = S-1;

    // Count the number of extant lineages
    // mutable size_t k = 0;
    // for(int h = 1; h < events.size(); h++){
    //     if(type == "extant leaf"){
    //         k++;
    //     }
    // }

    // Count event types :
    size_t k = 1;                       // Number of lineages
    size_t nb_fossil_leafs = 0;         // Number of fossil leafs
    size_t nb_sampled_ancestors = 0;    // Number of sampled ancestors
    size_t nb_occurrences = 0;          // Number of fossil occurrences
    size_t nb_branching_times = 0;      // Number of branching times

    // Recording the correction terms c to avoid divergence towards extreme values (outside of double precision)
    double log_correction = 0;
    double c;

    // Estimating the smallest sufficient N value for getting most of the probability density
    size_t N_optimal = 0;
    size_t N_optimal_tmp;

    // We start at the time of origin, supposedly the first time in the vector of events
    RbVector<double> Mt(N+1, 0.0);
    Mt[0] = 1;
    double thPlusOne = events[0].time;
    // if (verbose){std::cout << k << std::endl;}

    if(thPlusOne != start_age->getValue()) {
        if (verbose){std::cout << "WARNING : thPlusOne != start_age : " << thPlusOne << " != " << start_age->getValue() << " - type : " << events[0].type << std::endl;}
    };

    // Then we iterate over the next events
    for(int h = 0; h < events.size(); h++){

        // First, deal with the update on time period (th, thPlusOne)
        double th = events[h].time;
        std::string type = events[h].type;

        // Events older than the start age : accepted
        if(th > start_age->getValue()) {
            if (verbose){std::cout << "WARNING : th > start_age : " << th << " > " << start_age->getValue() << " -> type : " << type << std::endl;}

            // Time points before the start age
            if(type == "time point"){
                for(int i = 0; i < N+1; i++){
                    B[indxJ][i] = log(Mt[i]) + log_correction;
                }
                indxJ--;
            }
            // Other events before the start age : rejected
            else{
                if(returnLogLikelihood){
                    MatrixReal LogLikelihood_reject(1, 1, RbConstants::Double::neginf);
                    return LogLikelihood_reject;
                }
                else{
                    MatrixReal B_reject(S, (N + 1), RbConstants::Double::neginf);
                    return B_reject;
                }
            }

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
            // if (verbose){std::cout << "\nA : " << A << std::endl;}
            // if (verbose){std::cout << "\nMt : " << Mt << std::endl;}
            // if (verbose){std::cout << "\nlog(Mt[0]) : " << log(Mt[0]) << " - log(Mt[1]) : " << log(Mt[1]) << " - log(Mt[2]) : " << log(Mt[2]) << " - log(Mt[3]) : " << log(Mt[3]) << std::endl;}
            Mt = A * Mt;
        }

        // Second, deal with the update at punctual event th

        if(type == "time point"){
            // Combine the scaling factors for all the events until present
            double events_factor_log = log_correction;
            // if (verbose){std::cout << "log_correction : " << log_correction << std::endl;}
            if (ps != 0){
                events_factor_log += log(ps) * (nb_fossil_leafs + nb_sampled_ancestors);
                // if (verbose){std::cout << "log(ps) * (nb_fossil_leafs + nb_sampled_ancestors) : " << log(ps) * (nb_fossil_leafs + nb_sampled_ancestors) << std::endl;}
            }
            if (om != 0){
                events_factor_log += log(om) * nb_occurrences;
                // if (verbose){std::cout << "log(om) * nb_occurrences : " << log(om) * nb_occurrences << std::endl;}
            }
            if (birth != 0){
                events_factor_log += log(birth) * nb_branching_times;
                // if (verbose){std::cout << "log(birth) * nb_branching_times : " << log(birth) * nb_branching_times << std::endl;}
            }
            // if (verbose){std::cout << "events_factor_log : " << events_factor_log << std::endl;}
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = log(Mt[i]) + events_factor_log;
            }
            // if (verbose){std::cout << "Time point scaled log : log(Mt[0]) + events_factor_log : " << log(Mt[0]) + events_factor_log << " / log(Mt[N]) + events_factor_log : " << log(Mt[N]) + events_factor_log << std::endl;}
            indxJ--;
        }

        // if(type == "terminal removed"){
        //     for(int i = 0; i < N+1; i++){
        //         Mt[i] *= ps * rp;
        //     }
        //     k--;
        // }
        //
        // if(type == "fossil leaf"){
        //     for(int i = N; i > 0; i--){
        //         Mt[i] = Mt[i-1] * ps * (1-rp);
        //     }
        //     Mt[0] = 0;
        //     k--;
        // }

        // For now, we assume there is no information/labels on removal
        // @todo add fossil-removed/non-removed
        if(type == "fossil leaf"){
            for(int i = N; i > 0 ; i--){
                // Mt[i] = Mt[i-1] * ps * (1-rp) + ps * rp * Mt[i];
                Mt[i] = Mt[i-1] * (1-rp) + rp * Mt[i];
            }
            // Mt[0] *= ps * rp;
            Mt[0] *= rp;
            k--;
            nb_fossil_leafs++;
        }


        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                // Mt[i] *= ps * (1-rp);
                Mt[i] *= (1-rp);
            }
            nb_sampled_ancestors++;
        }

        // if(type == "occurrence removed"){
        //     for(int i = 0; i < N; i++){
        //         Mt[i] = Mt[i+1] * (i+1) * om * rp;
        //     }
        //  }
        //
        // if(type == "occurrence"){
        //     for(int i = 0; i < N+1; i++){
        //         Mt[i] *= (k+i) * om * (1-rp);
        //     }
        // }

        // For now, we assume there is no information/labels on removal
        // @todo add occurrence-removed/non-removed
        if(type == "occurrence"){
            for(int i = 0; i < N; i++){
                // Mt[i] = Mt[i+1] * (i+1) * om * rp + (k+i) * om * (1-rp) * Mt[i];
                Mt[i] = Mt[i+1] * (i+1) * rp + (k+i) * (1-rp) * Mt[i];
            }
            // Mt[N] *= (k+N) * om * (1-rp);
            Mt[N] *= (k+N) * (1-rp);
            nb_occurrences++;
         }
        // if(type == "occurrence"){
        //     // Correct probability vector to avoid divergence towards extreme values (outside of double precision) by scaling to a unit maximum value
        //     double max_Mt = *std::max_element(Mt.begin(), Mt.end());
        //     c = 1/max_Mt;
        //     // c=exp(1);
        //     for(int i = 0; i < N; i++){
        //         // Mt[i] = Mt[i+1] * (i+1) * om * rp + (k+i) * om * (1-rp) * Mt[i];
        //         Mt[i] = (Mt[i+1] * (i+1) * rp + (k+i) * (1-rp) * Mt[i]) * c;
        //    }
        //     // Mt[N] *= (k+N) * om * (1-rp);
        //     Mt[N] *= (k+N) * (1-rp) * c;
        //     nb_occurrences++;
            
        //     log_correction -= log(c);
        //     if (verbose){std::cout << "log(c) : " << log(c) << std::endl;}
        //  }


        if(type == "branching time"){
            // for(int i = 0; i < N+1; i++){
            //     Mt[i] *= birth;
            // }
            k++;
            nb_branching_times++;
        }
        
        // Correct probability vector to avoid divergence towards extreme values (outside of double precision) by scaling to a unit maximum value
        double max_Mt = *std::max_element(Mt.begin(), Mt.end());
        c = 1/max_Mt;
        for(int i = 0; i < N+1; i++){
            Mt[i] = Mt[i] * c;
        }
        log_correction -= log(c);
        
        // if (verbose){std::cout << "Event time : " << th << " - Event type : " << type << " -> log(c) : " << log(c) << " -> Mt[0] : " << Mt[0] << " / Mt[1] : " << Mt[1] << " / Mt[N] : " << Mt[N] << std::endl;}
        
        // Check that N is big enough
        if (Mt[N]>0.001){
            if (verbose){std::cout << "\nWARNING : Mt[N] contains a non-negligeable probability ( Mt[N] > max(Mt)/1000 ) -> you should increase N\n" << std::endl;}

            // In that case the optimal N value cannot be estimated because it is greater than the chosen one
            N_optimal = N+1;
        }
        else{
            N_optimal_tmp = N;
            while (Mt[N_optimal_tmp-1]<0.001){
                N_optimal_tmp--;
            }
            N_optimal = std::max(N_optimal, N_optimal_tmp);
        }
        // if (verbose){std::cout << "\nThe smallest sufficient N value ( such as Mt[N_optimal] < max(Mt)/1000 for all t ) is " << N_optimal << "\n" << std::endl;}

        thPlusOne = th;
    }
    
    // Give the estimated optimal N value
    if (N_optimal != N+1){
        if (verbose){std::cout << "\nThe smallest sufficient N value ( such as Mt[N_optimal] < max(Mt)/1000 for all t ) is " << N_optimal << "\n" << std::endl;}
    }
    else{
        if (verbose){std::cout << "\nWARNING : There is at least one time t such as Mt[N] contains a non-negligeable probability ( Mt[N] > max(Mt)/1000 ) -> you should increase N\n" << std::endl;}
    }

    if(returnLogLikelihood){
        double events_factor_log = log_correction;
        if (ps != 0){events_factor_log += log(ps) * (nb_fossil_leafs + nb_sampled_ancestors); }
        if (om != 0){events_factor_log += log(om) * nb_occurrences; }
        if (birth != 0) {events_factor_log += log(birth) * nb_branching_times; }
        
        double likelihood = 0;
        for(int i = 0; i < N+1; i++){
            likelihood += Mt[i] * pow(rh, k) * pow(1.0 - rh, i);
        }
        MatrixReal LogLikelihood(1, 1, log(likelihood) + events_factor_log);
        return LogLikelihood;
    }

    return B;
}



/**
 * Compute the joint log probability density of observations made down any time t, conditioned on the population
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
 * \param[in]    time_points            Times for which we want to compute the density.
 * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
 * \param[in]    verbose                If true displays warnings and information messages.
 * \param[in]    occurrence_ages        Vector of occurrence ages.
 * \param[in]    timeTree               Tree for ancestral populations size inference.
 *
 * \return    The matrix of log-Lt values through time.
 */
MatrixReal RevBayesCore::ComputeLnProbabilitiesBackwardsLt( const TypedDagNode<double> *start_age,
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
                                                            bool verbose,
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
          // node is a fossil leaf
          // For now, we assume there is no information/labels on removal
          // @todo add fossil-removed/non-removed
            events.push_back(Event(n.getAge(),"fossil leaf")) ;
            // events.push_back(Event(n.getAge(),"terminal non-removed")) ;
        }
        else if ( !n.isFossil() && n.isTip() )
        {
            // node is extant leaf : only their number is necessary to compute Lt and Mt
            // events.push_back(Event(n.getAge(),"extant leaf") ;
            // events.push_back(Event(n.getAge(), "terminal non-removed")) ;
            // if (verbose){std::cout << n.getSpeciesName() << std::endl;}
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
    // For now, we assume there is no information/labels on removal
    // @todo add occurrence-removed/non-removed
    events.push_back(Event(occurrence_ages[i],"occurrence")) ;
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
    MatrixReal B(S, (N + 1), RbConstants::Double::neginf);
    size_t indxJ = 0;

    // Count event types :
    size_t nb_fossil_leafs = 0;         // Number of fossil leafs
    size_t nb_sampled_ancestors = 0;    // Number of sampled ancestors
    size_t nb_occurrences = 0;          // Number of fossil occurrences
    size_t nb_branching_times = 0;      // Number of branching times

    // Recording the correction terms c to avoid divergence towards extreme values (outside of double precision)
    double log_correction = 0;
    double c;

    // Estimating the smallest sufficient N value for getting most of the probability density
    size_t N_optimal = 0;
    size_t N_optimal_tmp;

    // We start at time 0 with type "present" in the vector of events
    // size_t k = extant;
    double thMinusOne = events[0].time;
    RbVector<double> Lt(N+1, pow(rh, k));
    for ( int i = 0; i < (N + 1); i++){
        Lt[i] *= pow( 1.0-rh, i );
    }

    // We then iterate over the following events until finding the time of origin
    for(int h = 0; h < events.size(); h++){

        double th = events[h].time;
        std::string type = events[h].type;

        // Events older than the start age
        if(th > start_age->getValue()) {
            if (verbose){std::cout << "WARNING : th > start_age : " << th << " > " << start_age->getValue() << " -> type : " << type << std::endl;}

            // Time points older than the start age : accepted
            if(type == "time point"){
                double events_factor_log = log_correction;
                if (ps != 0){events_factor_log    += log(ps) * (nb_fossil_leafs + nb_sampled_ancestors);}
                if (om != 0){events_factor_log    += log(om) * nb_occurrences;}
                if (birth != 0){events_factor_log += log(birth) * nb_branching_times;}
                
                for(int i = 0; i < N+1; i++){
                    B[indxJ][i] = log(Lt[i]) + events_factor_log;
                }
                indxJ++;
            }
            // Other events older than the start age : rejected
            else{
                MatrixReal B_reject(S, (N + 1), RbConstants::Double::neginf);
                return B_reject;
            }

            continue;
        };

        // First deal with the update along (thMinusOne, th)
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
        if(type == "time point"){
            // Combine the scaling factors for all the events until present
            double events_factor_log = log_correction;
            // if (verbose){std::cout << "log_correction : " << log_correction << std::endl;}
            if (ps != 0){
                events_factor_log += log(ps) * (nb_fossil_leafs + nb_sampled_ancestors);
                // if (verbose){std::cout << "log(ps) * (nb_fossil_leafs + nb_sampled_ancestors) : " << log(ps) * (nb_fossil_leafs + nb_sampled_ancestors) << std::endl;}
            }
            
            if (om != 0){
                events_factor_log += log(om) * nb_occurrences;
                // if (verbose){std::cout << "log(om) * nb_occurrences : " << log(om) * nb_occurrences << std::endl;}
            }
            
            if (birth != 0){
                events_factor_log += log(birth) * nb_branching_times;
                // if (verbose){std::cout << "log(birth) * nb_branching_times : " << log(birth) * nb_branching_times << std::endl;}
            }

            // if (verbose){std::cout << "events_factor_log : " << events_factor_log  << " / Lt : " << Lt << std::endl;}
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = log(Lt[i]) + events_factor_log;
            }
            // if (verbose){std::cout << "Time point scaled log : log(Lt[0]) + events_factor_log : " << log(Lt[0]) + events_factor_log << " / log(Lt[N]) + events_factor_log : " << log(Lt[N]) + events_factor_log << std::endl;}
            indxJ++;
        }

        // if(type == "terminal removed"){
        //     for(int i = 0; i < N+1; i++){
        //         Lt[i] = Lt[i] * ps * rp;
        //     }
        //     k++;
        // }
        //
        // if(type == "fossil leaf"){
        //     for(int i = 0; i < N; i++){
        //         Lt[i] = Lt[i+1] * ps * (1.0-rp);
        //     }
        //     k++;
        // }

       // For now, we assume there is no information/labels on removal
       // @todo add fossil-removed/non-removed
       if(type == "fossil leaf"){
           for(int i = 0; i < N; i++){
               // Lt[i] = Lt[i] * ps * rp + Lt[i+1] * ps * (1.0-rp) ;
               Lt[i] = Lt[i] * rp + Lt[i+1] * (1.0-rp) ;
           }
           // Lt[N] = Lt[N] * ps * rp;
           Lt[N] = Lt[N] * rp;
           k++;
           nb_fossil_leafs++;
       }

        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                // Lt[i] = Lt[i] * ps * (1.0-rp);
                Lt[i] = Lt[i] * (1.0-rp);
            }
            nb_sampled_ancestors++;
        }

        // if(type == "occurrence removed"){
        //     for(int i = N; i > 0; i--){
        //         Lt[i] = Lt[i-1] * i * om * rp;
        //     }
        //     Lt[0] = 0;
        //  }

        // if(type == "occurrence"){
        //     for(int i = 0; i < N+1; i++){
        //         Lt[i] = Lt[i] * (k+i) * om * (1-rp);
        //     }
        // }

        // For now, we assume there is no information/labels on removal
       // @todo add occurrence-removed/non-removed
       if(type == "occurrence"){
           for(int i = N; i > 0; i--){
               // Lt[i] = Lt[i-1] * i * om * rp + Lt[i] * (k+i) * om * (1-rp);
               Lt[i] = Lt[i-1] * i * rp + Lt[i] * (k+i) * (1-rp);
           }
           // Lt[0] = Lt[0] * k * om * (1-rp);
           Lt[0] = Lt[0] * k * (1-rp);
           nb_occurrences++;
        }

        if(type == "branching time"){
            // for(int i = 0; i < N+1; i++){
            //     Lt[i] = Lt[i] * birth;
            // }
            k--;
            nb_branching_times++;
        }

        // Correct probability vector to avoid divergence towards extreme values (outside of double precision) by scaling to a unit maximum value
        double max_Lt = *std::max_element(Lt.begin(), Lt.end());
        c = 1/max_Lt;
        for(int i = 0; i < N+1; i++){
            Lt[i] = Lt[i] * c;
        }
        log_correction -= log(c);
        
        // if (verbose){std::cout << "Event time : " << th << " - Event type : " << type << " -> log(c) : " << log(c) << " -> Lt[0] : " << Lt[0] << " / Lt[1] : " << Lt[1] << " / Lt[N] : " << Lt[N] << std::endl;}
        
        // Check that N is big enough
        if (Lt[N]>0.001){
            if (verbose){std::cout << "\nWARNING : Lt[N] contains a non-negligeable probability ( Lt[N] > max(Lt)/1000 ) -> you should increase N\n" << std::endl;}

            // In that case the optimal N value cannot be estimated because it is greater than the chosen one
            N_optimal = N+1;
        }
        else{
            N_optimal_tmp = N;
            while (Lt[N_optimal_tmp-1]<0.001){
                N_optimal_tmp--;
            }
            N_optimal = std::max(N_optimal, N_optimal_tmp);
        }
        // if (verbose){std::cout << "\nThe smallest sufficient N value ( such as Lt[N_optimal] < max(Lt)/1000 for all t ) is " << N_optimal << "\n" << std::endl;}

        // if (verbose){std::cout << "Event time : " << th << " - Event type : " << type << " -> k : " << k << " -> Lt[0] : " << Lt[0] << " / Lt[1] : " << Lt[1] << " / Lt[N] : " << Lt[N] << std::endl;}

        thMinusOne = th;
    }

    // Give the estimated optimal N value
    if (N_optimal != N+1){
        if (verbose){std::cout << "\nThe smallest sufficient N value ( such as Lt[N_optimal] < max(Lt)/1000 for all t ) is " << N_optimal << "\n" << std::endl;}
    }
    else{
        if (verbose){std::cout << "\nWARNING : There is at least one time t such as Lt[N] contains a non-negligeable probability ( Lt[N] > max(Lt)/1000 ) -> you should increase N\n" << std::endl;}
    }

    return B;
}



// /**
//  * Compute the joint probability density of the observations made up to any time t and the population
//  * size at that time, as time decreases towards present : breadth-first forward traversal algorithm.
//  *
//  * \param[in]    start_age              Start age of the process.
//  * \param[in]    lambda                 Speciation rate.
//  * \param[in]    mu                     Extinction rate.
//  * \param[in]    psi                    Extinction sampling rate.
//  * \param[in]    omega                  Occurrence sampling rate.
//  * \param[in]    rho                    Sampling probability at present time.
//  * \param[in]    removalPr              Removal probability after sampling.
//  * \param[in]    maxHiddenLin           Algorithm accuracy (maximal number of hidden lineages).
//  * \param[in]    cond                   Condition of the process (none/survival/#Taxa).
//  * \param[in]    time_points         Times for which we want to compute the density.
//  * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
//  * \param[in]    timeTree               Tree for ancestral populations size inference.
//  *
//  * \return    The matrix of Mt values through time.
// */
MatrixReal RevBayesCore::ComputeLikelihoodsForwardsMtPiecewise(   const TypedDagNode<double> *start_age,
                                                                  const std::vector<double> &timeline,
                                                                  const std::vector< double > &lambda,
                                                                  const std::vector< double > &mu,
                                                                  const std::vector< double > &psi,
                                                                  const std::vector< double > &omega,
                                                                  const TypedDagNode<double> *rho,
                                                                  const std::vector< double > &removalPr,
                                                                  const TypedDagNode<long> *maxHiddenLin,
                                                                  const std::string& cond,
                                                                  const std::vector<double> &time_points,
                                                                  bool useOrigin,
                                                                  const std::vector<double> &occurrence_ages,
                                                                  const Tree &timeTree  )
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
            // node is fossil leaf
            // For now, we assume there is no information/labels on removal
            // @todo add fossil-removed/non-removed
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
        events.push_back(Event(occurrence_ages[i],"occurrence")) ;
    }
    //adding rate shift events to the events vector
    const std::vector<double> d = timeline;
    for ( int i=1; i < d.size(); i++)
    {
    events.push_back(Event(timeline[i],"rate shift")) ;
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


    const std::vector<double> birth = lambda;
    std::cout << "birth vector" << birth << std::endl;

    const std::vector<double> death = mu;
    const std::vector<double> ps = psi;
    const std::vector<double> om = omega;
    const double rh = rho->getValue();
    const std::vector<double> rp = removalPr;

    const long N = maxHiddenLin->getValue();
    // const RbVector<double> tau = time_points->getValue();

    const size_t S = tau.size();
    const std::vector<double> gamma = birth + death + ps + om;

    // Initialize an empty matrix and a cursor to write lines in this matrix
    MatrixReal B(S, (N + 1), 0.0);

    double birth_current = birth.back();
    double death_current = death.back();
    double ps_current = ps.back();
    double om_current = om.back();
    double rp_current = rp.back();
    double gamma_current = gamma.back();
    size_t indxJ = S-1;

    // Count the number of extant lineages
    // mutable size_t k = 0;
    // for(int h = 1; h < events.size(); h++){
    //     if(type == "extant leaf"){
    //         k++;
    //     }
    // }

    // Number of extant lineages
    k = 1;

    // We start at the time of origin, supposedly the first time in the vector of events
    RbVector<double> Mt(N+1, 0.0);
    Mt[0] = 1;
    double thPlusOne = events[0].time;
    // std::cout << "Event time : " << events[0].time << " - Event type : " << events[0].type << std::endl;
    // std::cout << k << std::endl;

    if(thPlusOne != start_age->getValue()) {
        std::cout << "WARNING : thPlusOne != start_age : " << thPlusOne << " != " << start_age->getValue() << std::endl;
        return B;
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
                A[i][i] = gamma_current * (k + i) * (th-thPlusOne);
                if (i < N) A[i][i+1] = -death_current * (i + 1) * (th-thPlusOne);
                if (i > 0) A[i][i-1] = -birth_current * (2*k + i - 1) * (th-thPlusOne);
            }

            RbMath::expMatrixPade(A, A, 4);
            Mt = A * Mt;
        }

        // Second, deal with the update at punctual event th
        std::string type = events[h].type;
        if(type == "rate shift"){
        //change rates appropriately, rate changes provided in timeline are in an ascending order
        size_t where = LocateTimeSliceIndex(th,timeline) - 1;

        std::cout <<"birth rate was " << birth_current << std::endl;
        std::cout <<"birth rate shifts to "<< birth[where] <<std::endl;
        birth_current = birth[where];
        death_current = death[where];
        ps_current = ps[where];
        om_current = om[where];
        rp_current = rp[where];
        gamma_current = gamma[where];

        }

        if(type == "time point"){
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = Mt[i];
            }
            indxJ--;
        }

        if(type == "terminal removed"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= ps_current * rp_current;
            }
            k--;
        }

        if(type == "fossil leaf"){
            for(int i = N; i > 0; i--){
                Mt[i] = Mt[i-1] * ps_current * (1-rp_current);
            }
            Mt[0] = 0;
            k--;
        }

        // For now, we assume there is no information/labels on removal
        // @todo add fossil-removed/non-removed
        // if(type == "fossil leaf"){
        //     for(int i = N; i > 0 ; i--){
        //         Mt[i] = Mt[i-1] * ps * (1-rp) + ps * rp * Mt[i];
        //     }
        //     Mt[0] *= ps * rp;
        //     k--;
        // }


        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= ps_current * (1-rp_current);
            }
        }

        if(type == "occurrence removed"){
            for(int i = 0; i < N; i++){
                Mt[i] = Mt[i+1] * (i+1) * om_current * rp_current;
            }
         }

        if(type == "occurrence"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= (k+i) * om_current * (1-rp_current);
            }
        }

        // For now, we assume there is no information/labels on removal
        // @todo add occurrence-removed/non-removed
        // if(type == "occurrence"){
        //     for(int i = 0; i < N; i++){
        //         Mt[i] = Mt[i+1] * (i+1) * om * rp + (k+i) * om * (1-rp) * Mt[i];
        //     }
        //     Mt[N] *= (k+N) * om * (1-rp);
        //  }


        if(type == "branching time"){
            for(int i = 0; i < N+1; i++){
                Mt[i] *= birth_current;
            }
            k++;
        }

        // std::cout << "Event time : " << th << " - Event type : " << type << " -> Mt[0] : " << Mt[0] << " / Mt[1] : " << Mt[1] << std::endl;

        thPlusOne = th;
    }

    // double likelihood = Mt[0];
    // for(int i = 1; i < N+1; i++){
    //     likelihood += Mt[i] * pow(rh,k) * pow(1.0 - rh,i);
    // }

    return B;
}
//
// /**
//  * Compute the probability density of observations made down any time t, conditioned on the population
//  * size at that time, as time increases towards the past : breadth-first backward traversal algorithm.
//  *
//  * \param[in]    start_age              Start age of the process.
//  * \param[in]    lambda                 Speciation rate.
//  * \param[in]    mu                     Extinction rate.
//  * \param[in]    psi                    Extinction sampling rate.
//  * \param[in]    omega                  Occurrence sampling rate.
//  * \param[in]    rho                    Sampling probability at present time.
//  * \param[in]    removalPr              Removal probability after sampling.
//  * \param[in]    maxHiddenLin           Algorithm accuracy (maximal number of hidden lineages).
//  * \param[in]    cond                   Condition of the process (none/survival/#Taxa).
//  * \param[in]    time_points            Times for which we want to compute the density.
//  * \param[in]    useOrigin              If true the start age is the origin time otherwise the root age of the process.
//  * \param[in]    timeTree               Tree for ancestral populations size inference.
//  *
//  * \return    The matrix of Lt values through time.
//  */
MatrixReal RevBayesCore::ComputeLikelihoodsBackwardsLtPiecewise(  const TypedDagNode<double> *start_age,
                                                    const std::vector< double > &timeline,
                                                    const std::vector< double > &lambda,
                                                    const std::vector< double > &mu,
                                                    const std::vector< double > &psi,
                                                    const std::vector< double > &omega,
                                                    const TypedDagNode<double> *rho,
                                                    const std::vector< double > &removalPr,
                                                    const TypedDagNode<long> *maxHiddenLin,
                                                    const std::string& cond,
                                                    const std::vector<double> &time_points,
                                                    bool useOrigin,
                                                    const std::vector<double> &occurrence_ages,
                                                    const Tree &timeTree  )
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
          // node is a fossil leaf
          // For now, we assume there is no information/labels on removal
          // @todo add fossil-removed/non-removed
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
    // For now, we assume there is no information/labels on removal
    // @todo add occurrence-removed/non-removed
    events.push_back(Event(occurrence_ages[i],"occurrence")) ;
    }
    //adding rate shift events to the events vector
    const std::vector<double> d = timeline;
    for ( int i=1; i < d.size(); i++)
    {
    events.push_back(Event(timeline[i],"rate shift")) ;
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





    const std::vector<double> birth = lambda;
    const std::vector<double> death = mu;
    const std::vector<double> ps = psi;
    const std::vector<double> om = omega;
    const double rh = rho->getValue();
    const std::vector<double> rp = removalPr;
    const long N = maxHiddenLin->getValue();
    // const RbVector<double> tau = time_points->getValue();

    const size_t S = tau.size();
    const std::vector<double> gamma = birth + death + ps + om;

    // Initialize an empty matrix and a cursor indxJ to write lines in this matrix as well as the parameter sets.
    MatrixReal B(S, (N + 1), 0.0);
    double birth_current = birth[0];
    double death_current = death[0];
    double ps_current = ps[0];
    double om_current = om[0];
    double rp_current = rp[0];
    double gamma_current = gamma[0];

    size_t indxJ = 0;


    // We start at time 0 with type "present" in the vector of events
    // size_t k = extant;
    double thMinusOne = events[0].time;
    RbVector<double> Lt(N+1, pow(rh, k));
    for ( int i = 0; i < (N + 1); i++){
        Lt[i] *= pow( 1.0-rh, i );
    }

    // std::cout << "Event time : " << events[0].time << " - Event type : " << events[0].type << std::endl;

    // We then iterate over the following events until finding the time of origin
    for(int h = 1; h < events.size(); h++){

        // First deal with the update along (thMinusOne, th)
        double th = events[h].time;

        if( th != thMinusOne ){

            MatrixReal A( (N+1), (N+1), 0.0 );
            for(int i = 0; i < (N + 1); i++){
              A[i][i] = -gamma_current * (k + i) * (th - thMinusOne);
              if (i < N) A[i][i+1] = birth_current * ( (2 * k) + i ) * (th - thMinusOne);
              if (i > 0) A[i][i-1] = death_current * i * (th - thMinusOne);
            }
            RbMath::expMatrixPade(A, A, 4);
            Lt = A * Lt;
        }

        // Second, deal with the update at time th
        std::string type = events[h].type;
        if(type == "rate shift"){
        //change rates appropriately
        size_t where = LocateTimeSliceIndex(th,timeline) ;
        std::cout <<"birth rate was " << birth_current << std::endl;
        std::cout <<"birth rate shifts to "<< birth[where] <<std::endl;

        birth_current = birth[where];
        death_current = death[where];
        ps_current = ps[where];
        om_current = om[where];
        rp_current = rp[where];
        gamma_current = gamma[where];
        }

        if(type == "time point"){
            for(int i = 0; i < N+1; i++){
                B[indxJ][i] = Lt[i];
            }
            indxJ++;
        }

        if(type == "terminal removed"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * ps_current * rp_current;
            }
            k++;
        }

        if(type == "fossil leaf"){
            for(int i = 0; i < N; i++){
                Lt[i] = Lt[i+1] * ps_current * (1.0-rp_current);
            }
            k++;
        }

       // For now, we assume there is no information/labels on removal
       // @todo add fossil-removed/non-removed
       // if(type == "fossil leaf"){
       //     for(int i = 0; i < N; i++){
       //         Lt[i] = Lt[i] * ps * rp + Lt[i+1] * ps * (1.0-rp) ;
       //     }
       //     Lt[N] = Lt[N] * ps * rp;
       //     k++;
       // }

        if(type == "sampled ancestor"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * ps_current * (1.0-rp_current);
            }
        }

        if(type == "occurrence removed"){
            for(int i = N; i > 0; i--){
                Lt[i] = Lt[i-1] * i * om_current * rp_current;
            }
            Lt[0] = 0;
         }

        if(type == "occurrence"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * (k+i) * om_current * (1-rp_current);
            }
        }

        // For now, we assume there is no information/labels on removal
       // @todo add occurrence-removed/non-removed
       // if(type == "occurrence"){
       //     for(int i = N; i > 0; i--){
       //         Lt[i] = Lt[i-1] * i * om * rp + Lt[i] * (k+i) * om * (1-rp);
       //     }
       //     Lt[0] = Lt[0] * k * om * (1-rp);
       //  }

        if(type == "branching time"){
            for(int i = 0; i < N+1; i++){
                Lt[i] = Lt[i] * birth_current;
            }
            k--;
        }

        // std::cout << "Event time : " << th << " - Event type : " << type << " -> Lt[0] : " << Lt[0] << " / Lt[1] : " << Lt[1] << std::endl;

        thMinusOne = th;
    }

    return B;
}


size_t RevBayesCore::LocateTimeSliceIndex(const double &t, const std::vector<double> &timeline)
{
    return std::distance(timeline.begin(),std::find(timeline.begin(),timeline.end(),t));
}

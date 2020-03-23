#include "ComputeLikelihoodsLtMtFunction.h"

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


/** ComputeLikelihoodsLtMtFunction of a MatrixReal Constructor
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
 * \param[in]    tau            Times for which we want to compute the density.
 * \param[in]    uo             If true sa is the origin time otherwise the root age of the process.
 * \param[in]    mt             If true computes densities with the Mt forward traversal algorithm otherwise uses Lt backward one.
 * \param[in]    tr             Tree for ancestral populations size inference.
*/
ComputeLikelihoodsLtMtFunction::ComputeLikelihoodsLtMtFunction( const TypedDagNode<double> *sa,
                                                                const TypedDagNode<double> *l,
                                                                const TypedDagNode<double> *m,
                                                                const TypedDagNode<double> *p,
                                                                const TypedDagNode<double> *o,
                                                                const TypedDagNode<double> *rho,
                                                                const TypedDagNode<double> *r,
                                                                const TypedDagNode<long> *n,

                                                                const std::string& cdt,
                                                                const std::vector<double> &tau,
                                                                bool uo,
                                                                bool mt,
                                                                const std::vector<double> &O,
                                                                const TypedDagNode<Tree> *tr) : TypedFunction<MatrixReal>( new MatrixReal(tau.size(), (n->getValue() + 1), 0.0) ),
    start_age( sa ),
    lambda( l ),
    mu( m ),
    psi( p ),
    omega( o ),
    rho( rho ),
    removalPr( r ),
    maxHiddenLin( n ),
    cond (cdt),
    time_points ( tau ),
    useOrigin (uo),
    useMt ( mt ),
    occurrence_ages( O ),
    timeTree (tr)
{
    this->addParameter( start_age );
    this->addParameter( lambda );
    this->addParameter( mu );
    this->addParameter( psi );
    this->addParameter( omega );
    this->addParameter( rho );
    this->addParameter( removalPr );
    this->addParameter( maxHiddenLin );
    this->addParameter( timeTree );

    poolTimes();
    update();
}



ComputeLikelihoodsLtMtFunction::~ComputeLikelihoodsLtMtFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


ComputeLikelihoodsLtMtFunction* ComputeLikelihoodsLtMtFunction::clone( void ) const
{
    return new ComputeLikelihoodsLtMtFunction( *this );
}



/**
 * Compute B, where you do the work, required
 *
 * \return A matrix of the joint probability densities of the tree with occurrences at each event
 */
void ComputeLikelihoodsLtMtFunction::update( void )
{
    if (useMt) {
        MatrixReal B = ComputeMt();
    }
    else {
        MatrixReal B = ComputeLt();
    }

    *this->value = B;
}

/**
 * Construct the vector containig all branching and sampling times + time points for which we want to compute the density.
 */
void ComputeLikelihoodsLtMtFunction::poolTimes( void ) const
{
    // get node/time variables
    size_t num_nodes = timeTree->getValue().getNumberOfNodes();

    extant = 0;

    // classify nodes
    events.clear();
    events.push_back(Event(start_age->getValue(), "origin"));

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = timeTree->getValue().getNode( i );

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
    for ( int i=0; i < occurrence_ages.size(); ++i)
    {
    events.push_back(Event(occurrence_ages[i],"occurrence non-removed")) ;
    }


    events.push_back(Event(0.0,"present time")) ;
}

/**
 * Compute the joint probability density of the observations made up to any time t and the population
 * size at that time, as time decreases towards present : breadth-first forward traversal algorithm.
 *
 * \return    The matrix of Mt values through time.
 */
MatrixReal ComputeLikelihoodsLtMtFunction::ComputeMt( void ) const
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

    size_t S = time_points.size();
    long N = maxHiddenLin->getValue();

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

        if(th > start_age->getValue()) {
            std::cout << "ERROR : th > start_age : " << th << " > " << start_age->getValue() << std::endl;
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
 * \return    The matrix of Lt values through time.
 */
MatrixReal ComputeLikelihoodsLtMtFunction::ComputeLt( void ) const
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

    // Initialize an empty matrix and a cursor indxJ to write lines in this matrix
    size_t S = time_points.size();
    long N = maxHiddenLin->getValue();
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

    return B;
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void ComputeLikelihoodsLtMtFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == start_age)
    {
        start_age = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == psi)
    {
        psi = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == omega)
    {
        omega = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == removalPr)
    {
        removalPr = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == maxHiddenLin)
    {
        maxHiddenLin = static_cast<const TypedDagNode< long >* >( newP );
    }
    else if (oldP == timeTree)
    {
        timeTree = static_cast<const TypedDagNode< Tree >* >( newP );
    }
}

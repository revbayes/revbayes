//
//  InferAncestralPopSizeFunction.cpp
//
//  Initiated by Rachel Warnock, Marc Manceau 30.01.2020.
//  Completed by Jérémy Andréoletti, Antoine Swaans 03.2020.
//
#include "InferAncestralPopSizeFunction.h"

// #include "ComputeLikelihoodsLtMtFunction.h"
#include "ComputeLikelihoodsLtMt.h"

#include <vector>
#include <iostream>
#include <cmath>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "Tree.h"
#include "TopologyNode.h"
// other includes

namespace RevBayesCore {
	class DagNode;
	class MatrixReal; } //tt

using namespace RevBayesCore;

/** InferAncestralPopSizeFunction
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)

 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    sa             Time of the origin/present/length of the process.
 * \param[in]    l              Speciation rate.
 * \param[in]    m              Extinction rate.
 * \param[in]    p              Fossil sampling rate.
 * \param[in]    o              Occurrence sampling rate.
 * \param[in]    rho            Sampling probability at present time.
 * \param[in]    r              Removal probability after sampling.
 * \param[in]    n              Algorithm accuracy (maximal number of hidden lineages).
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tau            Times for which we want to compute the density.
 * \param[in]    uo             If true t is the origin time otherwise the root age of the process.
 * \param[in]    tr             Tree for ancestral populations size inference.
 */
InferAncestralPopSizeFunction::InferAncestralPopSizeFunction( 	const TypedDagNode<double> *sa,
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

	                                                          	TypedDagNode<Tree> *tr) : TypedFunction<MatrixReal>( new MatrixReal(tau.size(), (n->getValue() + 1), 0.0) ),
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

	// poolTimes();
	update();
}

InferAncestralPopSizeFunction::~InferAncestralPopSizeFunction( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
InferAncestralPopSizeFunction* InferAncestralPopSizeFunction::clone( void ) const
{
    return new InferAncestralPopSizeFunction( *this );
}

/**
 * Compute Kt, where you do the work, required
 *
 * \return A density matrix of the number of hidden lineages through time
 */
void InferAncestralPopSizeFunction::update( void )
{
    size_t S = time_points.size();
    long N = maxHiddenLin->getValue();
    const Tree tree = timeTree->getValue();

    //bool useMt = false;
	MatrixReal B_Lt = RevBayesCore::ComputeLikelihoodsBackwardsLt(start_age, lambda, mu, psi, omega, rho, removalPr, maxHiddenLin, cond, time_points, useOrigin, tree);
	
	//useMt = true;
	MatrixReal B_Mt = RevBayesCore::ComputeLikelihoodsForwardsMt(start_age, lambda, mu, psi, omega, rho, removalPr, maxHiddenLin, cond, time_points, useOrigin, tree);

	// Realize the Hadamar Product of B_Lt and B_Mt
	MatrixReal D_Kt(S, (N + 1), 0.0);
	for (size_t i = 0; i < S; i++)
        {
            for (size_t j = 0; j < (N+1); j++)
            {
                D_Kt[i][j] = B_Lt[i][j] * B_Mt[i][j];
            }
        }

    // Transpose Lt
    MatrixReal B_Lt_T((N + 1), S, 0.0);
    RbMath::transposeMatrix(B_Lt, B_Lt_T);

    // Normalize
	D_Kt *= 1/(B_Lt_T*B_Mt)[0][0];

	*this->value = D_Kt; // this will eventually be the lk of the tree + occurrences

}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void InferAncestralPopSizeFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
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


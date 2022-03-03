//
//  InferAncestralPopSizeFunction.cpp
//
//  Initiated by Rachel Warnock, Marc Manceau 30.01.2020.
//  Completed by Jérémy Andréoletti, Antoine Zwaans 03.2020.
//
#include "InferAncestralPopSizeFunction.h"
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
 * \param[in]    vb             If true displays warnings and information messages.
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
                                                                const TypedDagNode< RevBayesCore::RbVector<double> > *O,
	                                                          		const std::vector<double> &tau,
                                                                bool uo,
                                                                bool vb,

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
    occurrences( O ),
    time_points ( tau ),
    useOrigin (uo),
    verbose (vb),
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
    this->addParameter( occurrences );
    this->addParameter( timeTree );

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
    const std::vector<double> occurrence_ages = occurrences->getValue();

    bool useMt;
		std::vector<double> timeline{0.0};
		std::vector<double> lambd{lambda->getValue()};
		std::vector<double> m{mu->getValue()};
		std::vector<double> ps{psi->getValue()};
		std::vector<double> omeg{omega->getValue()};
		std::vector<double> removalP{removalPr->getValue()};




	MatrixReal B_Lt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age->getValue(), timeline, lambd, m, ps, omeg, rho, removalP, maxHiddenLin, cond, time_points, useOrigin, useMt = false, verbose, occurrence_ages, tree);
	MatrixReal B_Mt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age->getValue(), timeline, lambd, m, ps, omeg, rho, removalP, maxHiddenLin, cond, time_points, useOrigin, useMt = true, verbose, occurrence_ages, tree);

	// Realize the normalized Hadamar Product of B_Lt and B_Mt
	MatrixReal D_Kt(S, (N + 1), 0.0);
    for (size_t j = 0; j < S; j++){
        // Compute the mean log-probability (later substracted in order to avoid extreme values, out of the double precision)
        double B_Lt_Mt_log_mean = 0.0;
        for (size_t i = 0; i < (N+1); i++)
        {
            B_Lt_Mt_log_mean += B_Lt_log[j][i] + B_Mt_log[j][i];
        }
        B_Lt_Mt_log_mean /= N;

        // Compute the normalization factor
        double norm = 0.0;
        for (size_t i = 0; i < (N+1); i++)
        {
            norm += exp(B_Lt_log[j][i] + B_Mt_log[j][i] - B_Lt_Mt_log_mean);
        }

        for (size_t i = 0; i < (N+1); i++)
        {
            D_Kt[j][i] = exp(B_Lt_log[j][i] + B_Mt_log[j][i] - B_Lt_Mt_log_mean)/norm;
        }
    }

	*this->value = D_Kt;

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

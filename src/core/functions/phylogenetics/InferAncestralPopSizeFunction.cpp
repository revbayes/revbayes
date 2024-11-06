#include "InferAncestralPopSizeFunction.h"
#include "ComputeLikelihoodsLtMt.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "Tree.h"
#include "TopologyNode.h"
// other includes

namespace RevBayesCore {
	class DagNode;
	class MatrixReal; }

using namespace RevBayesCore;

/** InferAncestralPopSizeFunctionPiecewise
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * We delegate most parameters to the base class and initialize the members.
 * \param[in]    sa                        Start age of the process.
 * \param[in]    inspeciation              Speciation/birth rate(s).
 * \param[in]    inextinction              Extinction/death rate(s).
 * \param[in]    inserialsampling          Serial sampling rate(s).
 * \param[in]    inoccurrence              Occurrence sampling rate(s).
 * \param[in]    ineventsampling           Sampling probability at present time.
 * \param[in]    intreatment               Probabilit(y|ies) of death upon sampling (treatment).
 * \param[in]    n                         Maximum number of hidden lineages (algorithm accuracy).
 * \param[in]    cdt                       Condition of the process (survival/survival2).
 * \param[in]    O                         Vector of occurrence ages.
 * \param[in]    tau            		   Time points at which we compute the density.
 * \param[in]    vb                        Verbose
 * \param[in]    tr                        Tree for which ancestral pop. size has to be computed.
 */

InferAncestralPopSizeFunction::InferAncestralPopSizeFunction( 	const TypedDagNode<double> *sa,
                                      							const TypedDagNode<double> *inspeciation,
	                                                          	const TypedDagNode<double> *inextinction,
	                                                          	const TypedDagNode<double> *inserialsampling,
	                                                          	const TypedDagNode<double> *inoccurrence,
	                                                          	const TypedDagNode<double> *ineventsampling,
	                                                          	const TypedDagNode<double> *intreatment,
		                                                        const TypedDagNode<std::int64_t> *n,
	                                                          	const std::string& cdt,
                                                                const TypedDagNode< RevBayesCore::RbVector<double> > *O,
	                                                          	const std::vector<double> &tau,
                                                                bool vb,
	                                                          	TypedDagNode<Tree> *tr) :

    TypedFunction<MatrixReal>( new MatrixReal(tau.size(), (n->getValue() + 1), 0.0) ),

    start_age( sa ),
    lambda( inspeciation ),
    mu( inextinction ),
    psi( inserialsampling ),
    omega( inoccurrence ),
    rho( ineventsampling ),
    removalPr( intreatment ),
    maxHiddenLin( n ),
    cond (cdt),
    occurrences( O ),
    time_points ( tau ),
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





InferAncestralPopSizeFunction::~InferAncestralPopSizeFunction( void )
{}





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
 * Compute the density matrix of the number of hidden lineages through time
 *
 * \return A density matrix of the number of hidden lineages through time
 */
void InferAncestralPopSizeFunction::update( void )
{
    size_t S = time_points.size();
    std::int64_t N = maxHiddenLin->getValue();

    const Tree& tree = timeTree->getValue();
    const double sa = start_age->getValue();
    const std::vector<double>& occurrence_ages = occurrences->getValue();

    bool useMt;
	std::vector<double> timeline{0.0};
	std::vector<double> lambd{lambda->getValue()};
	std::vector<double> m{mu->getValue()};
	std::vector<double> ps{psi->getValue()};
	std::vector<double> omeg{omega->getValue()};
	std::vector<double> removalP{removalPr->getValue()};

	MatrixReal B_Lt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(sa, timeline, lambd, m, ps, omeg, rho, removalP, maxHiddenLin, cond,
                                                                          time_points, useMt = false, verbose, occurrence_ages, tree);
    MatrixReal B_Mt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(sa, timeline, lambd, m, ps, omeg, rho, removalP, maxHiddenLin, cond,
                                                                          time_points, useMt = true, verbose, occurrence_ages, tree);

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
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void InferAncestralPopSizeFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == start_age)
    {
        start_age    = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == lambda)
    {
        lambda       = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == mu)
    {
        mu           = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == psi)
    {
        psi          = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == omega)
    {
        omega        = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == rho)
    {
        rho          = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == removalPr)
    {
        removalPr    = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == maxHiddenLin)
    {
        maxHiddenLin = static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }
    else if (oldP == timeTree)
    {
        timeTree     = static_cast<const TypedDagNode< Tree >* >( newP );
    }
}

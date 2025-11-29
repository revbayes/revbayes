#include "InferAncestralPopSizeFunctionPiecewise.h"
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
 * \param[in]    tau            	       Time points at which we compute the density.
 * \param[in]    t                         Tree for which ancestral pop. size has to be computed.
 * \param[in]    ht                        Rate interval change times of the piecewise constant process.
 */



InferAncestralPopSizeFunctionPiecewise::InferAncestralPopSizeFunctionPiecewise( 	const TypedDagNode<double> *sa,
                      																														const DagNode *inspeciation,
                                          																				const DagNode *inextinction,
                                          																				const DagNode *inserialsampling,
                                          																				const DagNode *inoccurrence,
                                          																				const DagNode *ineventsampling,
                                          																				const DagNode *intreatment,
                                          																				const TypedDagNode<std::int64_t> *n,
                                          																				const std::string& cdt,
								                                          										    const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                          																				const std::vector<double> &tau,
                                          																				TypedDagNode<Tree> *tr,
																																									const TypedDagNode< RbVector<double> > *ht) :

    TypedFunction<MatrixReal>( new MatrixReal(tau.size(), (n->getValue() + 1), 0.0) ),

    start_age( sa ),
    maxHiddenLin( n ),
    cond (cdt),
    occurrences( O ),
    time_points ( tau ),
    timeTree (tr),
		interval_times(ht)

{
	homogeneous_lambda   = NULL;
	heterogeneous_lambda = NULL;

	homogeneous_mu       = NULL;
	heterogeneous_mu     = NULL;

	homogeneous_psi      = NULL;
	heterogeneous_psi    = NULL;

	homogeneous_o        = NULL;
	heterogeneous_o      = NULL;

	homogeneous_rho      = NULL;

	homogeneous_r        = NULL;
	heterogeneous_r      = NULL;


	std::vector<double> times = timeline;
	std::vector<double> times_sorted_ascending = times;

	sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );

	if ( times != times_sorted_ascending )
	{
			throw(RbException("Rate change times must be provided in ascending order."));
	}

	addParameter( interval_times );

	heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
	homogeneous_lambda   = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

	addParameter( homogeneous_lambda );
	addParameter( heterogeneous_lambda );

	heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
	homogeneous_mu   = dynamic_cast<const TypedDagNode<double >*>(inextinction);

	addParameter( homogeneous_mu );
	addParameter( heterogeneous_mu );

	heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inserialsampling);
	homogeneous_psi   = dynamic_cast<const TypedDagNode<double >*>(inserialsampling);

	addParameter( homogeneous_psi );
	addParameter( heterogeneous_psi );

	heterogeneous_r = dynamic_cast<const TypedDagNode<RbVector<double> >*>(intreatment);
	homogeneous_r   = dynamic_cast<const TypedDagNode<double >*>(intreatment);

	addParameter( homogeneous_r );
	addParameter( heterogeneous_r );

	heterogeneous_o = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inoccurrence);
	homogeneous_o   = dynamic_cast<const TypedDagNode<double >*>(inoccurrence);

	addParameter( homogeneous_o );
	addParameter( heterogeneous_o );

	homogeneous_rho = dynamic_cast<const TypedDagNode<double >*>(ineventsampling);

	addParameter( homogeneous_rho );


	//check that lengths of vector arguments are sane
	if ( heterogeneous_lambda != NULL && !(interval_times->getValue().size() == heterogeneous_lambda->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument lambda must have one more element than timeline."));
	}

	if ( heterogeneous_mu != NULL && !(interval_times->getValue().size() == heterogeneous_mu->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument mu must have one more element than timeline."));
	}

	if ( heterogeneous_psi != NULL && !(interval_times->getValue().size() == heterogeneous_psi->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument psi must have one more element than timeline."));
	}

	if ( heterogeneous_r != NULL && !(interval_times->getValue().size() == heterogeneous_r->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument r must have one more element than timeline."));
	}

	if ( heterogeneous_o != NULL && !(interval_times->getValue().size() == heterogeneous_o->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument o must have one more element than timeline."));
	}


	update();
}





InferAncestralPopSizeFunctionPiecewise::~InferAncestralPopSizeFunctionPiecewise( void )
{}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
InferAncestralPopSizeFunctionPiecewise* InferAncestralPopSizeFunctionPiecewise::clone( void ) const
{
    return new InferAncestralPopSizeFunctionPiecewise( *this );
}





/**
 * Compute the density matrix of the number of hidden lineages through time.
 *
 * \return A density matrix of the Maximal number of hidden lineages through time
 */
void InferAncestralPopSizeFunctionPiecewise::update( void )
{
		updateVectorParameters();
    size_t S = time_points.size();
    std::int64_t N = maxHiddenLin->getValue();

    const Tree& tree = timeTree->getValue();
    const std::vector<double>& occurrence_ages = occurrences->getValue(); 

		bool useMt;
		bool verbose;
		MatrixReal B_Lt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age->getValue(), timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, time_points, useMt = false, verbose = true, occurrence_ages, tree);
		MatrixReal B_Mt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age->getValue(), timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, time_points, useMt = true, verbose = true, occurrence_ages, tree);
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





void InferAncestralPopSizeFunctionPiecewise::updateVectorParameters( void ) const
{
    // clean and get timeline
    timeline.clear();
    timeline = interval_times->getValue();


    timeline.insert(timeline.begin(),0.0);

    // clean all the sets
    lambda.clear();
    mu.clear();
    psi.clear();
    r.clear();
    omega.clear();
    psi_event.clear();

    // Get vector of birth rates
    if ( heterogeneous_lambda != NULL )
    {
      lambda = heterogeneous_lambda->getValue();
    }
    else
    {
      lambda = std::vector<double>(timeline.size(),homogeneous_lambda->getValue());
    }

    // Get vector of death rates
    if ( heterogeneous_mu != NULL )
    {
      mu = heterogeneous_mu->getValue();
    }
    else
    {
      mu = std::vector<double>(timeline.size(),homogeneous_mu->getValue());
    }

    // Get vector of serial sampling rates
    if ( heterogeneous_psi != NULL )
    {
      psi = heterogeneous_psi->getValue();
    }
    else
    {
      psi = std::vector<double>(timeline.size(),homogeneous_psi->getValue());
    }

    // Get vector of conditional death upon sampling probabilities
    if ( heterogeneous_r != NULL )
    {
      r = heterogeneous_r->getValue();
    }
    else
    {
      r = std::vector<double>(timeline.size(),homogeneous_r->getValue());
    }

    // Get vector of occcurrence rates.
    if ( heterogeneous_o != NULL )
    {
      omega = heterogeneous_o->getValue();
    }
    else
    {
      omega = std::vector<double>(timeline.size(),homogeneous_o->getValue());
    }

		psi_event = std::vector<double>(timeline.size(),0.0);
    if ( homogeneous_rho != NULL )
    {
        // User specified the sampling fraction at the present
        psi_event[0] = homogeneous_rho->getValue();
    }
    else
    {
        // set the final sampling to one (for sampling at the present)
        psi_event[0] = 1.0;
  	}
}





/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void InferAncestralPopSizeFunctionPiecewise::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == start_age)
    {
        start_age            = static_cast<const TypedDagNode< double >* >( newP );
    }
	else if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
	else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda   = static_cast<const TypedDagNode<double>* >( newP );
    }
	else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu     = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
	else if (oldP == homogeneous_mu)
    {
        homogeneous_mu       = static_cast<const TypedDagNode<double>* >( newP );
    }
	else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi    = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
	else if (oldP == homogeneous_psi)
    {
        homogeneous_psi      = static_cast<const TypedDagNode<double>* >( newP );
    }
	else if (oldP == heterogeneous_o)
    {
        heterogeneous_o      = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_o)
    {
        homogeneous_o        = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho      = static_cast<const TypedDagNode<double>* >( newP );
    }
	else if (oldP == heterogeneous_r)
    {
        heterogeneous_r      = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_r)
    {
        homogeneous_r        = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == maxHiddenLin)
    {
        maxHiddenLin         = static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }
    else if (oldP == timeTree)
    {
        timeTree             = static_cast<const TypedDagNode< Tree >* >( newP );
    }
}

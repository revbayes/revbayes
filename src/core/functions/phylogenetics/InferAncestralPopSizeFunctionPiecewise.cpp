//
//  InferAncestralPopSizeFunctionPiecewise.cpp
//
//  Initiated by Rachel Warnock, Marc Manceau 30.01.2020.
//  Completed by Jérémy Andréoletti, Antoine Zwaans 03.2020.
//
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
// other includes

namespace RevBayesCore {
	class DagNode;
	class MatrixReal; } //tt

using namespace RevBayesCore;

/** InferAncestralPopSizeFunctionPiecewise
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
InferAncestralPopSizeFunctionPiecewise::InferAncestralPopSizeFunctionPiecewise( 	const TypedDagNode<double> *sa,
                                                  						const DagNode *inspeciation,
	                                                          	const DagNode *inextinction,
	                                                          	const DagNode *inserialsampling,
	                                                          	const DagNode *inoccurrence,
	                                                          	const DagNode *ineventsampling,
	                                                          	const DagNode *intreatment,
	                                                          	const TypedDagNode<long> *n,
	                                                          	const std::string& cdt,
                                                              const TypedDagNode< RevBayesCore::RbVector<double> > *O,
	                                                          	const std::vector<double> &tau,
	                                                          	bool uo,
	                                                          	TypedDagNode<Tree> *tr,
																															const TypedDagNode< RbVector<double> > *ht) : TypedFunction<MatrixReal>( new MatrixReal(tau.size(), (n->getValue() + 1), 0.0) ),
    start_age( sa ),
    maxHiddenLin( n ),
    cond (cdt),
    occurrences( O ),
    time_points ( tau ),
    useOrigin (uo),
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
	heterogeneous_rho    = NULL;

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
	homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

	addParameter( homogeneous_lambda );
	addParameter( heterogeneous_lambda );

	heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
	homogeneous_mu = dynamic_cast<const TypedDagNode<double >*>(inextinction);

	addParameter( homogeneous_mu );
	addParameter( heterogeneous_mu );

	heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inserialsampling);
	homogeneous_psi = dynamic_cast<const TypedDagNode<double >*>(inserialsampling);

	addParameter( homogeneous_psi );
	addParameter( heterogeneous_psi );

	heterogeneous_r = dynamic_cast<const TypedDagNode<RbVector<double> >*>(intreatment);
	homogeneous_r = dynamic_cast<const TypedDagNode<double >*>(intreatment);

	addParameter( homogeneous_r );
	addParameter( heterogeneous_r );

	heterogeneous_o = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inoccurrence);
	homogeneous_o = dynamic_cast<const TypedDagNode<double >*>(inoccurrence);

	addParameter( homogeneous_o );
	addParameter( heterogeneous_o );

	heterogeneous_rho = dynamic_cast<const TypedDagNode<RbVector<double> >*>(ineventsampling);
	homogeneous_rho = dynamic_cast<const TypedDagNode<double >*>(ineventsampling);

	addParameter( homogeneous_rho );
	addParameter( heterogeneous_rho );


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

	if ( heterogeneous_rho != NULL && !(interval_times->getValue().size() == heterogeneous_rho->getValue().size() - 1) )
	{
		throw(RbException("If provided as a vector, argument rho must have one more element than timeline."));
	}

	update();
}

InferAncestralPopSizeFunctionPiecewise::~InferAncestralPopSizeFunctionPiecewise( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

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
 * Compute Kt, where you do the work, required
 *
 * \return A density matrix of the number of hidden lineages through time
 */
void InferAncestralPopSizeFunctionPiecewise::update( void )
{
		updateVectorParameters();
    size_t S = time_points.size();
    long N = maxHiddenLin->getValue();

    const Tree tree = timeTree->getValue();
    const std::vector<double> occurrence_ages = occurrences->getValue();

		bool useMt;
		bool verbose;
	MatrixReal B_Lt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age, timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, time_points, useOrigin, useMt = false, verbose = true, occurrence_ages, tree);
	std::cout << "LT is ok, go to Mt" << std::endl;
	MatrixReal B_Mt_log = RevBayesCore::ComputeLnProbabilityDensitiesOBDP(start_age, timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, time_points, useOrigin, useMt = true, verbose = true, occurrence_ages, tree);
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

    // Get vector of event sampling probabilities
    if ( heterogeneous_rho != NULL )
    {
      // User has specified psi_event_0,...,psi_event_{l-1}
      psi_event = heterogeneous_rho->getValue();
    }
    else
    {
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
        start_age = static_cast<const TypedDagNode< double >* >( newP );
    }
		if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
		else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
		else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
		else if (oldP == homogeneous_mu)
    {
        homogeneous_mu = static_cast<const TypedDagNode<double>* >( newP );
    }
		else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
		else if (oldP == homogeneous_psi)
    {
        homogeneous_psi = static_cast<const TypedDagNode<double>* >( newP );
    }
		else if (oldP == heterogeneous_o)
    {
        heterogeneous_o = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_o)
    {
        homogeneous_o = static_cast<const TypedDagNode<double>* >( newP );
    }
		else if (oldP == heterogeneous_rho)
    {
        heterogeneous_rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho = static_cast<const TypedDagNode<double>* >( newP );
    }
		else if (oldP == heterogeneous_r)
    {
        heterogeneous_r = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_r)
    {
        homogeneous_r = static_cast<const TypedDagNode<double>* >( newP );
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

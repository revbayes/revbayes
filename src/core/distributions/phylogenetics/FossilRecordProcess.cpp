#include "FossilRecordProcess.h"

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "StochasticNode.h"
#include "MatrixReal.h"
#include "RbException.h"

using namespace RevBayesCore;


/**
 * Constructor. The observation half of the factored FBD-range model: conditions on a
 * birth-death-range SKELETON (matrix-first prototype: a MatrixReal-valued stochastic node
 * whose distribution derives from AbstractFossilizedBirthDeathRangeProcess) and contributes
 * the fossil-occurrence (reporting) log-density. Data (this node's value) is the taxon vector.
 */
FossilRecordProcess::FossilRecordProcess(const DagNode *sk, const std::string &rep, const std::vector<Taxon> &t) :
    TypedDistribution< RbVector<Taxon> >( new RbVector<Taxon>(t) ),
    skeleton_node( sk ),
    skeleton( NULL ),
    reporting( rep ),
    taxa( t )
{
    addParameter( skeleton_node );
    resolveSkeleton();

    if ( skeleton == NULL )
    {
        throw RbException("dnFossilRecord requires a fossilized birth-death range skeleton (e.g. dnFBDRP) as its first argument.");
    }

    // push the reporting model onto the skeleton: one source of truth for the tau1 support + reporting term
    skeleton->setReportingModel( reporting );
}


FossilRecordProcess* FossilRecordProcess::clone( void ) const
{
    return new FossilRecordProcess( *this );
}


/** Downcast the skeleton stochastic node's distribution to the shared range base (the type-bridge). */
void FossilRecordProcess::resolveSkeleton( void )
{
    skeleton = NULL;
    const StochasticNode<MatrixReal> *sn = dynamic_cast<const StochasticNode<MatrixReal> *>( skeleton_node );
    if ( sn != NULL )
    {
        skeleton = const_cast<AbstractFossilizedBirthDeathRangeProcess *>(
                       dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &sn->getDistribution() ) );
    }
}


double FossilRecordProcess::computeLnProbability( void )
{
    // psi, timeline, b/d, tau1/tau_last and the reporting model all live on the skeleton;
    // this is exactly the term computeLnProbabilityRanges omits when report_internally is false.
    return skeleton->computeLnFossilTotal();
}


void FossilRecordProcess::redrawValue( void )
{
    // Clamped-data distribution for now. Forward Poisson-thin simulator is TODO (needs the
    // dating timeline + the Inf-encoding reader fix; see project tasks).
    (*this->value) = taxa;
}


void FossilRecordProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == skeleton_node )
    {
        skeleton_node = newP;
        resolveSkeleton();
    }
}

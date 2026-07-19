#include "FossilRecordProcess.h"

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "StochasticNode.h"
#include "MatrixReal.h"
#include "Tree.h"
#include "TimeInterval.h"
#include "RbException.h"

using namespace RevBayesCore;


/** Downcast a ranges stochastic node (MatrixReal-valued dnFBDRP or Tree-valued dnFBDSP) to the shared range base. */
AbstractFossilizedBirthDeathRangeProcess* FossilRecordProcess::rangesOf(const DagNode *n)
{
    const StochasticNode<MatrixReal> *mn = dynamic_cast<const StochasticNode<MatrixReal> *>( n );
    if ( mn != NULL )
    {
        return const_cast<AbstractFossilizedBirthDeathRangeProcess *>(
                   dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &mn->getDistribution() ) );
    }
    const StochasticNode<Tree> *tn = dynamic_cast<const StochasticNode<Tree> *>( n );
    if ( tn != NULL )
    {
        return const_cast<AbstractFossilizedBirthDeathRangeProcess *>(
                   dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &tn->getDistribution() ) );
    }
    return NULL;
}


/**
 * Constructor. The observation half of the factored FBD-range model: conditions on a birth-death
 * range process (a dnFBDRP or dnFBDSP node) and contributes the fossil-occurrence (reporting)
 * log-density. The occurrences are read from the ranges' taxa, so they need not be supplied again
 * here; this node's value is that same taxon vector.
 */
FossilRecordProcess::FossilRecordProcess(const DagNode *rn, const std::string &rep) :
    TypedDistribution< RbVector<Taxon> >( new RbVector<Taxon>() ),
    ranges_node( rn ),
    ranges( rangesOf( rn ) ),
    reporting( rep ),
    taxa()
{
    if ( ranges == NULL )
    {
        throw RbException("dnFossilRecord requires a fossilized birth-death range process (e.g. dnFBDRP) as its first argument.");
    }

    addParameter( ranges_node );

    // the reporting term reads psi and the timeline off the range process, not through its value
    const std::vector<const DagNode*> &rp = ranges->getRangeParameters();
    for (std::vector<const DagNode*>::const_iterator it = rp.begin(); it != rp.end(); ++it)
    {
        addParameter( *it );
    }

    // the occurrences are the range process's taxa
    taxa = ranges->getTaxa();
    *this->value = RbVector<Taxon>( taxa );

    // push the reporting model onto the ranges: one source of truth for the tau1 support + reporting term
    ranges->setReportingModel( reporting );
}


FossilRecordProcess* FossilRecordProcess::clone( void ) const
{
    return new FossilRecordProcess( *this );
}


/** Downcast the ranges stochastic node's distribution to the shared range base (the type-bridge). */
void FossilRecordProcess::resolveRanges( void )
{
    ranges = rangesOf( ranges_node );
}


/**
 * Clamp: the likelihood reads the occurrences from the range process, not from this value, so a
 * clamped record must match that range process. Validate here (matched by name) and reject a
 * mismatch rather than let the wrong data pass silently. TODO: once redrawValue simulates a record,
 * the clamped value becomes the data read directly.
 */
void FossilRecordProcess::setValue(RbVector<Taxon> *v, bool force)
{
    const std::vector<Taxon> &ref = ranges->getTaxa();

    if ( v->size() != ref.size() )
    {
        throw RbException("dnFossilRecord: the clamped record has a different number of taxa than the range process.");
    }

    for (size_t i = 0; i < v->size(); ++i)
    {
        const Taxon &c = (*v)[i];

        const Taxon *r = NULL;
        for (size_t j = 0; j < ref.size(); ++j)
        {
            if ( ref[j].getName() == c.getName() ) { r = &ref[j]; break; }
        }
        if ( r == NULL )
        {
            throw RbException("dnFossilRecord: taxon '" + c.getName() + "' is not in the range process.");
        }
        if ( c.getOccurrences() != r->getOccurrences() )
        {
            throw RbException("dnFossilRecord: the occurrences for taxon '" + c.getName() + "' do not match the range process.");
        }
    }

    TypedDistribution< RbVector<Taxon> >::setValue( v, force );
}


double FossilRecordProcess::computeLnProbability( void )
{
    // psi, timeline, b/d, tau1/tau_last and the reporting model all live on the range process;
    // this is exactly the term computeLnProbabilityRanges omits when report_internally is false.
    return ranges->computeLnFossilTotal();
}


void FossilRecordProcess::redrawValue( void )
{
    // Clamped-data distribution for now. Forward Poisson-thin simulator is TODO (needs the
    // dating timeline + the Inf-encoding reader fix; see project tasks).
    (*this->value) = taxa;
}


void FossilRecordProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == ranges_node )
    {
        ranges_node = newP;
        resolveRanges();
    }
}

#include "ComputeEmpiricalWithinSpeciesVariancesFunction.h"

#include <cmath>
#include <string>

#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ComputeEmpiricalWithinSpeciesVariancesFunction::ComputeEmpiricalWithinSpeciesVariancesFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *s, const std::vector<Taxon> &ta, bool lt, RevBayesCore::TypedDagNode<double>* dv ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    tau( t ),
    data( d ),
    site( s ),
    log_transformed( lt ),
    taxa( ta ),
    default_var( dv )
{
    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( data );
    addParameter( s );

    resetWithinSpeciesVariances();
    update();
}


ComputeEmpiricalWithinSpeciesVariancesFunction::~ComputeEmpiricalWithinSpeciesVariancesFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ComputeEmpiricalWithinSpeciesVariancesFunction* ComputeEmpiricalWithinSpeciesVariancesFunction::clone( void ) const
{
    return new ComputeEmpiricalWithinSpeciesVariancesFunction( *this );
}


double ComputeEmpiricalWithinSpeciesVariancesFunction::computeMeanForSpecies(const std::string &name, size_t index)
{

    double mean = 0.0;
    double num_samples = getNumberOfSamplesForSpecies(name);

    const ContinuousCharacterData &d = data->getValue();

    for (size_t i=0; i<taxa.size(); ++i)
    {

        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            const ContinuousTaxonData& taxon = d.getTaxonData( t.getName() );
            mean += taxon.getCharacter(index);

        }

    }

    // normalize
    mean /= num_samples;


    return mean;
}


double ComputeEmpiricalWithinSpeciesVariancesFunction::computeWithinSpeciesVariance(const std::string &name, size_t index, bool log_tr)
{

    double num_samples = getNumberOfSamplesForSpecies(name);
    double var = 0.0;

    if ( num_samples > 1 )
    {
        double mean = computeMeanForSpecies(name, index);

        const ContinuousCharacterData &d = data->getValue();

        for (size_t i=0; i<taxa.size(); ++i)
        {

            const Taxon &t = taxa[i];
            if ( name == t.getSpeciesName() )
            {
                const ContinuousTaxonData& taxon = d.getTaxonData( t.getName() );
                var += (taxon.getCharacter(index) - mean) * (taxon.getCharacter(index) - mean);

            }

        }

        // normalize
        var /= num_samples;

        if ( log_tr == true )
        {
            var = log(var);
        }
    }
    else
    {
        var = default_var->getValue();
        //static_cast< double >( default_var );
    }

    return var;
}


double ComputeEmpiricalWithinSpeciesVariancesFunction::getNumberOfSamplesForSpecies(const std::string &name)
{

    double num_samples = 0.0;

    for (size_t i=0; i<taxa.size(); ++i)
    {

        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ++num_samples;
        }

    }

    return num_samples;
}


void ComputeEmpiricalWithinSpeciesVariancesFunction::resetWithinSpeciesVariances( void )
{
    size_t num_tips = tau->getValue().getNumberOfTips();

    // check if the vectors need to be resized
    within_species_variance     = std::vector<double>(num_tips, 0);

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    size_t site_index = site->getValue()-1;

    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {

            const std::string &name = (*it)->getName();
            within_species_variance[(*it)->getIndex()] = computeWithinSpeciesVariance(name,site_index,log_transformed);

        }

    }

}

void ComputeEmpiricalWithinSpeciesVariancesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == data)
    {
        data = static_cast<const TypedDagNode<ContinuousCharacterData>* >( newP );
    }
    else if (oldP == site)
    {
        site = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
}


void ComputeEmpiricalWithinSpeciesVariancesFunction::update( void )
{
    RbVector<double> &v = *value;
    const Tree &tree = tau->getValue();

    size_t num_tips = tree.getNumberOfTips();

    if ( v.size() != num_tips )
    {
        v.resize( num_tips );
    }

    for (size_t i=0; i<num_tips; ++i)
    {
        v[i] = within_species_variance[i];
    }
}

#include "ComputeWithinSpeciesVarianceFromCharacterDataFunction.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbException.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ComputeWithinSpeciesVarianceFromCharacterDataFunction::ComputeWithinSpeciesVarianceFromCharacterDataFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *vs, const TypedDagNode<std::int64_t> *ns, MISSING_TREATMENT mtr ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    data( d ),
    variance_site( vs ),
    num_sample_site( ns )
{
    missing_var_treatment = mtr;

    // add the lambda parameter as a parent
    addParameter( data );
    addParameter( variance_site );
    addParameter( num_sample_site );

    resetWithinSpeciesVariances();
    update();
}


ComputeWithinSpeciesVarianceFromCharacterDataFunction::~ComputeWithinSpeciesVarianceFromCharacterDataFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ComputeWithinSpeciesVarianceFromCharacterDataFunction* ComputeWithinSpeciesVarianceFromCharacterDataFunction::clone( void ) const
{
    return new ComputeWithinSpeciesVarianceFromCharacterDataFunction( *this );
}


double ComputeWithinSpeciesVarianceFromCharacterDataFunction::computeWithinSpeciesVariance(const std::string &name, size_t v_site_index, size_t n_site_index)
{

    double num_samples = getNumberOfSamplesForSpecies(name, n_site_index);
    double var = 0.0;

    if ( num_samples > 1 )
    {
        const ContinuousCharacterData &d = data->getValue();

        const ContinuousTaxonData& taxon = d.getTaxonData( name );
        var = taxon.getCharacter(v_site_index);

        // normalize
        var /= num_samples;

    }
    else
    {
        // change here with options MISSING_TREATMENT
        if ( missing_var_treatment == MEAN )
        {
            var = computeMeanWithinSpeciesVariance();
        }
        else if ( missing_var_treatment == MEDIAN )
        {
            var = computeMedianWithinSpeciesVariance();
        }
        else if ( missing_var_treatment == NONE )
        {
            var = -1.0;
        }
        else
        {
            throw RbException( "Argument missingVarianceTreatment must be one of \"mean\", \"median\" or \"none\"" );
        }
        //static_cast< double >( default_var );
    }

    return var;
}


double ComputeWithinSpeciesVarianceFromCharacterDataFunction::getNumberOfSamplesForSpecies(const std::string &name, size_t n_site_index)
{

    const ContinuousCharacterData &d = data->getValue();
    const ContinuousTaxonData& taxon = d.getTaxonData( name );

    double num_samples = taxon.getCharacter(n_site_index);

    return num_samples;
}


std::vector<std::string> ComputeWithinSpeciesVarianceFromCharacterDataFunction::getAlphabeticalSpeciesNames(void)
{
    const ContinuousCharacterData &d = data->getValue();
    const std::vector<Taxon> &taxa = d.getTaxa();

    std::vector<std::string> species_names;

    for (size_t i=0; i<taxa.size(); ++i)
    {

      const std::string &name = taxa[i].getSpeciesName();
      species_names.push_back(name);

    }

    species_names.erase(std::unique(species_names.begin(), species_names.end()), species_names.end());
    sort( species_names.begin(), species_names.end() );

    return species_names;
}


double ComputeWithinSpeciesVarianceFromCharacterDataFunction::computeMeanWithinSpeciesVariance( void )
{

    // some of the sites may have been excluded
    size_t v_site_index = variance_site->getValue()-1;
    size_t n_site_index = num_sample_site->getValue()-1;

    double mean_var              = 0.0;
    size_t num_taxa_multi_sample = 0.0;

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        double num_samples = getNumberOfSamplesForSpecies(name, n_site_index);

        if ( num_samples > 1 )
        {
            mean_var += computeWithinSpeciesVariance(name, v_site_index, n_site_index);
            num_taxa_multi_sample++;
        }

    }

    mean_var /= num_taxa_multi_sample;
    return mean_var;
}

double ComputeWithinSpeciesVarianceFromCharacterDataFunction::computeMedianWithinSpeciesVariance( void )
{
    size_t v_site_index = variance_site->getValue()-1;
    size_t n_site_index = num_sample_site->getValue()-1;

    std::vector<double> vars = std::vector<double>(0, 0);

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        double num_samples = getNumberOfSamplesForSpecies(name, n_site_index);

        if ( num_samples > 1 )
        {
            double var = computeWithinSpeciesVariance(name, v_site_index, n_site_index);
            vars.push_back(var);
        }

    }

    sort( vars.begin(), vars.end() );

    double med_var = 0.0;
    if (vars.size() % 2 != 0) // if the number of elements is odd
    {
        med_var = vars[vars.size() / 2];
    }
    else                      // if the number of elements is odd
    {
        med_var = (vars[(vars.size() - 1) / 2] + vars[vars.size() / 2]) / 2.0;
    }

    return med_var;
}


void ComputeWithinSpeciesVarianceFromCharacterDataFunction::resetWithinSpeciesVariances( void )
{

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    // check if the vectors need to be resized
    within_species_variance     = std::vector<double>(num_species, 0);

    // some of the sites may have been excluded
    size_t v_site_index = variance_site->getValue()-1;
    size_t n_site_index = num_sample_site->getValue()-1;

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        within_species_variance[i] = computeWithinSpeciesVariance(name, v_site_index, n_site_index);

    }

}

void ComputeWithinSpeciesVarianceFromCharacterDataFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == data)
    {
        data = static_cast<const TypedDagNode<ContinuousCharacterData>* >( newP );
    }
    else if (oldP == variance_site)
    {
        variance_site = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else if (oldP == num_sample_site)
    {
        num_sample_site = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}


void ComputeWithinSpeciesVarianceFromCharacterDataFunction::update( void )
{
    RbVector<double> &v = *value;


    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    if ( v.size() != num_species )
    {
        v.resize( num_species );
    }

    for (size_t i=0; i<num_species; ++i)
    {
        v[i] = within_species_variance[i];
    }
}

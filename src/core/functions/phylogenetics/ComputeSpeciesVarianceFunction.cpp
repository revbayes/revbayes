#include "ComputeSpeciesVarianceFunction.h"

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

ComputeSpeciesVarianceFunction::ComputeSpeciesVarianceFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *s, const std::vector<Taxon> &ta, MISSING_TREATMENT mtr, bool bool_vom ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    data( d ),
    site( s ),
    taxa( ta ),
    compute_VarOfMean( bool_vom )
{
    missing_var_treatment = mtr;

    // add the lambda parameter as a parent
    addParameter( data );
    addParameter( site );

    reset();
    update();
}


ComputeSpeciesVarianceFunction::~ComputeSpeciesVarianceFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ComputeSpeciesVarianceFunction* ComputeSpeciesVarianceFunction::clone( void ) const
{
    return new ComputeSpeciesVarianceFunction( *this );
}


double ComputeSpeciesVarianceFunction::computeMeanForSpecies(const std::string &name, size_t index)
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


double ComputeSpeciesVarianceFunction::computeTipErrorOrVarianceForSpecies(const std::string &name, size_t index)
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
        var /= ( num_samples - 1 );

        // if standard error of mean trait is desired
        if ( compute_VarOfMean )
        {
            var /= num_samples;
        }


    }
    else
    {
        // change here with options MISSING_TREATMENT
        if ( missing_var_treatment == MEAN )
        {
            var = computeMeanErrorOrVarianceAcrossSpecies();
        }
        else if ( missing_var_treatment == MEDIAN )
        {
            var = computeMedianErrorOrVarianceAcrossSpecies();
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


double ComputeSpeciesVarianceFunction::getNumberOfSamplesForSpecies(const std::string &name)
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


std::vector<std::string> ComputeSpeciesVarianceFunction::getAlphabeticalSpeciesNames(void)
{

    std::vector<std::string> species_names;

    for (size_t i=0; i<taxa.size(); ++i)
    {

      const std::string &name = taxa[i].getSpeciesName();
      species_names.push_back(name);

    }

    species_names.erase(unique(species_names.begin(), species_names.end()), species_names.end());
    sort( species_names.begin(), species_names.end() );

    return species_names;
}


double ComputeSpeciesVarianceFunction::computeMeanErrorOrVarianceAcrossSpecies( void )
{

    // some of the sites may have been excluded
    size_t site_index = site->getValue()-1;

    double mean_var              = 0.0;
    double num_species_multi_sample = 0.0;

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        double num_samples = getNumberOfSamplesForSpecies(name);

        if ( num_samples > 1 )
        {
            mean_var += computeTipErrorOrVarianceForSpecies(name,site_index);
            num_species_multi_sample++;
        }

    }

    mean_var /= num_species_multi_sample;
    return mean_var;
}

double ComputeSpeciesVarianceFunction::computeMedianErrorOrVarianceAcrossSpecies( void )
{
    size_t site_index = site->getValue()-1;

    std::vector<double> vars = std::vector<double>(0, 0);

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        double num_samples = getNumberOfSamplesForSpecies(name);

        if ( num_samples > 1 )
        {
            double var = computeTipErrorOrVarianceForSpecies(name,site_index);
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


void ComputeSpeciesVarianceFunction::reset( void )
{

    std::vector<std::string> species_names = getAlphabeticalSpeciesNames();
    size_t num_species = species_names.size();

    // check if the vectors need to be resized
    within_species_variance     = std::vector<double>(num_species, 0);

    // some of the sites may have been excluded
    size_t site_index = site->getValue()-1;

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        within_species_variance[i] = computeTipErrorOrVarianceForSpecies(name,site_index);

    }

}

void ComputeSpeciesVarianceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == data)
    {
        data = static_cast<const TypedDagNode<ContinuousCharacterData>* >( newP );
    }
    else if (oldP == site)
    {
        site = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}


void ComputeSpeciesVarianceFunction::update( void )
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

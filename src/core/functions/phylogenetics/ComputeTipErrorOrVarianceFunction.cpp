#include "ComputeTipErrorOrVarianceFunction.h"

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

ComputeTipErrorOrVarianceFunction::ComputeTipErrorOrVarianceFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *s, const std::vector<Taxon> &ta, MISSING_TREATMENT mtr, bool err ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    data( d ),
    site( s ),
    taxa( ta ),
    compute_SEM( err )
{
    missing_var_treatment = mtr;

    // add the lambda parameter as a parent
    addParameter( data );
    addParameter( site );

    reset();
    update();
}


ComputeTipErrorOrVarianceFunction::~ComputeTipErrorOrVarianceFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ComputeTipErrorOrVarianceFunction* ComputeTipErrorOrVarianceFunction::clone( void ) const
{
    return new ComputeTipErrorOrVarianceFunction( *this );
}


double ComputeTipErrorOrVarianceFunction::computeMeanForSpecies(const std::string &name, size_t index)
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


double ComputeTipErrorOrVarianceFunction::computeTipErrorOrVarianceForSpecies(const std::string &name, size_t index)
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

        // if standard error of mean trait is desired
        if ( compute_SEM )
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


double ComputeTipErrorOrVarianceFunction::getNumberOfSamplesForSpecies(const std::string &name)
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


std::vector<std::string> ComputeTipErrorOrVarianceFunction::getAlphabeticalSpeciesNames(void)
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


double ComputeTipErrorOrVarianceFunction::computeMeanErrorOrVarianceAcrossSpecies( void )
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

double ComputeTipErrorOrVarianceFunction::computeMedianErrorOrVarianceAcrossSpecies( void )
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


void ComputeTipErrorOrVarianceFunction::reset( void )
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

void ComputeTipErrorOrVarianceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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


void ComputeTipErrorOrVarianceFunction::update( void )
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

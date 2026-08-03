#include "GetContinuousCharacterAsVectorFunction.h"

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

GetContinuousCharacterAsVectorFunction::GetContinuousCharacterAsVectorFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *s, VECTOR_ORDER ord ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    data( d ),
    site_index( s )
{
    order_by = ord;

    // add the lambda parameter as a parent
    addParameter( data );
    addParameter( site_index );

    reset();
    update();
}


GetContinuousCharacterAsVectorFunction::~GetContinuousCharacterAsVectorFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



GetContinuousCharacterAsVectorFunction* GetContinuousCharacterAsVectorFunction::clone( void ) const
{
    return new GetContinuousCharacterAsVectorFunction( *this );
}


double GetContinuousCharacterAsVectorFunction::getContinuousCharacter(const std::string &name, size_t site_index)
{

    double cont_char = 0.0;


    const ContinuousCharacterData &d = data->getValue();

    const ContinuousTaxonData& taxon = d.getTaxonData( name );
    cont_char = taxon.getCharacter(site_index);



    return cont_char;
}


std::vector<std::string> GetContinuousCharacterAsVectorFunction::getAlphabeticalSpeciesNames(void)
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

void GetContinuousCharacterAsVectorFunction::reset( void )
{

    std::vector<std::string> species_names;
    if ( order_by == ALPHABETICAL )
    {
        species_names = getAlphabeticalSpeciesNames();

    }
    else
    {
        throw RbException( "Currently, we only support ordering by species names alphabetically." );
    }

    size_t num_species = species_names.size();

    // check if the vectors need to be resized
    continuous_character = std::vector<double>(num_species, 0);

    // some of the sites may have been excluded
    size_t s = site_index->getValue()-1;

    for (size_t i=0; i<species_names.size(); ++i)
    {

        std::string name = species_names[i];
        continuous_character[i] = getContinuousCharacter(name, s);

    }

}

void GetContinuousCharacterAsVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == data)
    {
        data = static_cast<const TypedDagNode<ContinuousCharacterData>* >( newP );
    }
    else if (oldP == site_index)
    {
        site_index = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}


void GetContinuousCharacterAsVectorFunction::update( void )
{
    RbVector<double> &v = *value;


    std::vector<std::string> species_names;
    if ( order_by == ALPHABETICAL )
    {
        species_names = getAlphabeticalSpeciesNames();

    }
    else
    {
        throw RbException( "Currently, we only support ordering by species names alphabetically." );
    }

    size_t num_species = species_names.size();

    if ( v.size() != num_species )
    {
        v.resize( num_species );
    }

    for (size_t i=0; i<num_species; ++i)
    {
        v[i] = continuous_character[i];
    }
}

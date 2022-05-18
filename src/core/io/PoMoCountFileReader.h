#ifndef PoMoCountFileReader_H
#define PoMoCountFileReader_H

#include "DelimitedDataReader.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousDiscreteCharacterData.h"
#include "PoMoState.h"

#include <string>
#include <vector>

namespace RevBayesCore {


    /**
     * Reader for allele count filesfor use with PoMo matrices.
     * Same format as read by IQ-Tree.
     * Description: http://www.iqtree.org/doc/Polymorphism-Aware-Models/
     * The first line contains COUNTSFILE  then the number of populations and the number of sites.
     * The following line contains the header: chromosome, site position, population names.
     * The following lines contain the states.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Bastien Boussau)
     * @since 2015-03-03, version 1.0
     *
     */
    class PoMoCountFileReader : public DelimitedDataReader {

    public:
        
        enum FORMAT { PoMo, NaturalNumbers };

        PoMoCountFileReader(const std::string& fn, const size_t vps = 10, FORMAT f=FORMAT::PoMo);
        PoMoCountFileReader(const PoMoCountFileReader& r);
        virtual                                                 ~PoMoCountFileReader();
        
        PoMoCountFileReader&                                    operator=(const PoMoCountFileReader& n);                                         //!< Assignment operator

        
        const size_t                                            getNumberOfPopulations( void );
        const size_t                                            getNumberOfSites( void );
        const AbstractHomologousDiscreteCharacterData&          getMatrix( void );
        const size_t                                            getVirtualPopulationSize( void );

    protected:

        size_t                                                  number_of_populations;
        size_t                                                  number_of_sites;
        size_t                                                  virtual_population_size;
        std::vector<string>                                     names;
        AbstractHomologousDiscreteCharacterData*                matrix;
        FORMAT                                                  data_format;
    };

}

#endif

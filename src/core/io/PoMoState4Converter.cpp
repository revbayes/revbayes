#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "PoMoState4Converter.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "PoMoState4.h"
#include "RbException.h"
#include "Taxon.h"


using namespace RevBayesCore;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
PoMoState4Converter::PoMoState4Converter( void )
{

}


/**
 * Data converter from DNA into PoMoState4.
 *
 * This function concverts a DNA matrix into a PoMoState4 matrix of given virtualPopulationSize,
 * using the given mapping between sequence name and species name.
 */
HomologousDiscreteCharacterData<PoMoState4>* PoMoState4Converter::convertData(
                                                                                const AbstractHomologousDiscreteCharacterData &d,
                                                                                const unsigned int virtualPopulationSize,
                                                                                const std::map<std::string, std::string> sequenceNameToSpeciesName) {
    HomologousDiscreteCharacterData<PoMoState4>* data = new HomologousDiscreteCharacterData<PoMoState4> ();
    //First, build a vector of frequencies according to the PoMo model
    std::vector<double> tempFreq (5, 0.0);
    std::vector< std::vector<double> > frequencies ( 4+ (virtualPopulationSize-1) * 6, tempFreq);
    double unit = 1.0/(double)virtualPopulationSize;
    frequencies[0][0]=1.0;
    frequencies[1][1]=1.0;
    frequencies[2][2]=1.0;
    frequencies[3][3]=1.0;
    //AC
    size_t start = 3;
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][0]=(i-start)*unit;
        frequencies[i][1]=(virtualPopulationSize-(i-start))*unit;
       // std::cout << "i: "<< i << " ; "<<(i-start)*unit << " ; "<<(virtualPopulationSize-(i-start))*unit <<std::endl;
    }
    //AG
    start = 3 + virtualPopulationSize-1;
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][0]=(i-start)*unit;
        frequencies[i][2]=(virtualPopulationSize-(i-start))*unit;
    }
    //AT
    start = 3 + 2*(virtualPopulationSize-1);
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][0]=(i-start)*unit;
        frequencies[i][3]=(virtualPopulationSize-(i-start))*unit;
    }
    //CG
    start = 3 + 3*(virtualPopulationSize-1);
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][1]=(i-start)*unit;
        frequencies[i][2]=(virtualPopulationSize-(i-start))*unit;
    }
    //CT
    start = 3 + 4*(virtualPopulationSize-1);
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][1]=(i-start)*unit;
        frequencies[i][3]=(virtualPopulationSize-(i-start))*unit;
    }
    //GT
    start = 3 + 5*(virtualPopulationSize-1);
    for (size_t i = start+1; i< start + 1+ (virtualPopulationSize-1) ; ++i) {
        frequencies[i][2]=(i-start)*unit;
        frequencies[i][3]=(virtualPopulationSize-(i-start))*unit;
    }

   /* std::cout << "Frequencies: "<<std::endl;
    for (size_t j = 0; j < 30; j++) {
    for (size_t i = 0 ; i < frequencies[j].size(); ++i) {
        std::cout << frequencies[j][i] << ",";
    }
        PoMoState4* pol = new PoMoState4(virtualPopulationSize);
        pol->setState((size_t) (j+1) );
        std::cout << "  "<<pol->getStringValue() <<std::endl;
    }*/

    //Second, go through the map to find all the species present in the data
    std::map<std::string, std::vector<std::string> > speciesNameToSequenceNames;
    for (std::map<std::string, std::string>::const_iterator it = sequenceNameToSpeciesName.begin(); it != sequenceNameToSpeciesName.end(); it++) {
        if ( speciesNameToSequenceNames.find(it->second) != speciesNameToSequenceNames.end() ) {
            speciesNameToSequenceNames[it->second].push_back(it->first);
        }
        else {
            std::vector<std::string> temp (1, it->first);
            speciesNameToSequenceNames[it->second] = temp;
        }
    }
    //std::cout << "Found "<<speciesNameToSequenceNames.size() << " species." <<std::endl;
    std::vector<double> counts (5, 0.0); //A C G T -
    std::vector<double> countsBackup (5, 0.0); //A C G T -
    for (std::map<std::string, std::vector<std::string> >::const_iterator it = speciesNameToSequenceNames.begin(); it != speciesNameToSequenceNames.end(); it++) {
        DiscreteTaxonData<PoMoState4> tax (it->first);
        for (size_t c = 0; c < d.getNumberOfCharacters(); ++c) {
            for (std::vector<std::string>::const_iterator seq = it->second.begin(); seq != it->second.end(); seq++) {
                size_t index = d.getIndexOfTaxon(*seq);
                const DnaState* ch = static_cast< const DnaState*>(&(d.getCharacter(index, c)));
                size_t chIndex = 0;
                if (ch->getStringValue()=="A") {
                    chIndex = 0;
                }
                else if (ch->getStringValue()=="C") {
                    chIndex = 1;
                }
                else if (ch->getStringValue()=="G") {
                    chIndex = 2;
                }
                else if (ch->getStringValue()=="T") {
                    chIndex = 3;
                }
                else if (ch->getStringValue()=="-") {
                    chIndex = 4;
                }
                //std::cout << "ch->getStateIndex(): " << ch->getStateIndex() << " aka "<< ch->getStringValue()<< " ch->getState(): " << ch->getState() << std::endl;
                counts[chIndex]++;
            }
            //Now we have all the counts for this species,
            //we need to use these counts to build a PoMoState4
            tax.addCharacter(*convertCounts(counts, virtualPopulationSize, frequencies) );
            //Resetting counts
            counts = countsBackup;
        }
        data->addTaxonData(tax);
    }




    return data;
}


PoMoState4* PoMoState4Converter::convertCounts(std::vector<double> &counts,
                                                           const unsigned int virtualPopulationSize,
                                                           std::vector< std::vector<double> > &frequencies) {

    //First, normalize the counts vector
    double sum = 0.0;
    for (size_t i = 0 ; i < counts.size(); ++i)
    {
      //  std::cout << counts[i] << ",";
        sum += counts[i];
    }
    //std::cout << "  ";
    for (size_t i = 0 ; i < counts.size(); ++i)
    {
        counts[i] /= sum;
       // std::cout << counts[i] << ",";
    }
    //If the site is all gaps
    if (counts[4] == 1.0)
    {
        PoMoState4* pol = new PoMoState4(4);
        pol->setStateByIndex((size_t)0);
        return pol;
    }
    //Now, compare the counts vector to the frequencies vector, to find the most alike frequency.
    size_t index=0;
    double minDiff = 10000000000;
    double diff = 0.0;
    for (size_t i = 0 ; i < frequencies.size(); ++i)
    {
        diff=0.0;
        for (size_t j = 0 ; j < counts.size(); ++j)
        {
            diff += (frequencies[i][j] -counts[j]) * (frequencies[i][j] -counts[j]);
        }
        if (diff < minDiff) {
            minDiff = diff;
            index = i;
           /* for (size_t j = 0 ; j < counts.size(); ++j) {
                std::cout << frequencies[i][j]<< " , ";
            }*/
        }
    }
    PoMoState4* pol = new PoMoState4(virtualPopulationSize);
    pol->setStateByIndex((size_t) (index+1) );
    std::cout << "polgetstringvalue:  "<<pol->getStringValue() <<std::endl;;
    return pol;
}

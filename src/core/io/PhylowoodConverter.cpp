//
//  PhylowoodConverter.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "DelimitedDataReader.h"
#include "NclReader.h"
#include "PhylowoodConverter.h"
#include "RbFileManager.h"
#include "Tree.h"
#include "TimeAtlas.h"
#include "TimeAtlasDataReader.h"
#include "GeographicArea.h"
#include "TopologyNode.h"

using namespace RevBayesCore;

RevBayesCore::PhylowoodConverter::PhylowoodConverter(const path &sfn, const path &tfn, const path &gfn, const path &pfn, double b, const std::string& ct, const std::string& bt) :
    stateFilename(sfn),
    treeFilename(tfn),
    geoFilename(gfn),
    phwFilename(pfn),
    chartype(ct),
    bgtype(bt),
    burn(b)
{
    
    convert();
}

PhylowoodConverter::~PhylowoodConverter(void)
{
    delete atlas;
    delete dat;
}

void PhylowoodConverter::convert(void) {
    
    // get atlas info
    TimeAtlasDataReader* tadr = new TimeAtlasDataReader(geoFilename);
    atlas = new TimeAtlas(tadr);
    numAreas = atlas->getNumAreas();
    if (bgtype == "Range")
    {
        num_states = std::pow(2, (double)atlas->getNumAreas());
        makeBits();
    }
    else if (bgtype == "Area") {
        num_states = atlas->getNumAreas();
    }

    numEpochs = atlas->getNumEpochs();
    
    // get tree info
    NclReader reader = NclReader();
    std::vector<Tree> trees;
    std::vector<Tree*>* tmp = reader.readTimeTrees( treeFilename );
    if (tmp->size() > 0)
    {
        tree = (*tmp)[0];
    }
    delete tmp;
    
    num_nodes = tree->getNumberOfNodes();
    
    // get state info
    dat = new DelimitedDataReader(stateFilename);
    makeMarginalAreaProbs();
    
    std::string s = buildPhylowoodString();
    
    createDirectoryForFile( phwFilename );

    // open the stream to the file
    std::fstream outStream;
    outStream.open( phwFilename.string(), std::fstream::out );
    outStream << s;
    outStream.close();
}



std::string PhylowoodConverter::buildPhylowoodString(void)
{
    std::stringstream nhxStrm;
    
    // begin nexus file
    nhxStrm << "#NEXUS" << "\n\n";
    
    // phylowood settings block
    nhxStrm << "Begin phylowood;\n";
    nhxStrm << "\tdrawtype pie\n";
    nhxStrm << "\tmodeltype biogeography\n";
    nhxStrm << "\tareatype discrete\n";
    nhxStrm << "\tmaptype clean\n";
    nhxStrm << "\tpieslicestyle even\n";
    nhxStrm << "\tpiefillstyle outwards\n";
    nhxStrm << "\ttimestart -" << tree->getRoot().getAge() << "\n";
    nhxStrm << "\tmarkerradius " << 200 << "\n";
    nhxStrm << "\tminareaval " << 0.1 << "\n";
    nhxStrm << "End;\n\n";
    
    // bayarea-fig block
    nhxStrm << "Begin bayarea-fig;\n";
    nhxStrm << "\tmapheight\t100\n";
    nhxStrm << "\tmapwidth\t150\n";
    nhxStrm << "\tcanvasheight\t2000\n";
    nhxStrm << "\tcanvaswidth\t1000\n";
    nhxStrm << "\tminareaval\t0.15\n";
    nhxStrm << "\tareacolors black\n";
    nhxStrm << "\tareatypes";
    for (unsigned i = 0; i < numAreas; i++)
        nhxStrm << " 1";
    nhxStrm << "\n";
    nhxStrm << "\tareanames Default\n";
    nhxStrm << "End;\n\n";
    
    // taxa block
    nhxStrm << "Begin taxa;\n";
    nhxStrm << "\tDimensions ntax=" << tree->getNumberOfTips() << ";\n";
    nhxStrm << "\tTaxlabels\n";
    for (unsigned i = 0; i < tree->getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getName() << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // geo block
    nhxStrm << "Begin geo;\n";
    nhxStrm << "\tDimensions ngeo=" << numAreas << ";\n";
    nhxStrm << "\tCoords\n";
    
    for (unsigned i = 0; i < numAreas; i++)
    {
        double lat = atlas->getAreas()[numEpochs-1][i]->getLatitude();
        double lon = atlas->getAreas()[numEpochs-1][i]->getLongitude();
        nhxStrm << "\t\t" << i << "\t" << lat << "\t" << lon;
        if (i != 0)
            nhxStrm << ",";
        nhxStrm << "\n";

    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // tree block
    nhxStrm << "Begin trees;\n";
    nhxStrm << "\tTranslate\n";
    for (unsigned i = 0; i < tree->getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getIndex() << "\t" << p->getName();
            if (i < (tree->getNumberOfNodes() - 1))
                nhxStrm << ",";
            nhxStrm << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    
    // write tree string
    std::string treeStr = "";
    treeStr = buildExtendedNewick(&tree->getRoot());
//    std::cout << treeStr << "\n";
//    std::cout << "nhxStr\n" << treeStr << "\n";
    nhxStrm << "tree TREE1 = " << treeStr << ";\n";
    nhxStrm << "End;\n";
    
//    std::cout << "[";
//    for (size_t i = 0; i < numAreas; i++)
//    {
//        if (i != 0)
//            std::cout << ",";
//        std::cout << (double)childCharacterCounts[0][i] / numSamples;
//    }
//    std::cout << "]\n";
    
    return nhxStrm.str();
}


void PhylowoodConverter::makeMarginalAreaProbs(void) {

    // read data
    std::vector<std::vector<std::string> > stateChars = dat->getChars();
    size_t numSamples = stateChars.size() - 1;
    std::string standardStates = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

    // get state counts/node
    marginalStartProbs.resize(num_nodes);
    marginalEndProbs.resize(num_nodes);
    for (size_t i = 0; i < num_nodes; i++)
    {
        marginalStartProbs[i].resize(numAreas, 0.0);
        marginalEndProbs[i].resize(numAreas, 0.0);
    }

    // get data
    int skip = 1;
    std::vector<std::vector<unsigned> > idx;
    for (size_t i = 0; i < stateChars.size(); i++)
    {
        for (size_t j = 0; j < stateChars[i].size(); j++)
        {
            
            if (j < skip)
            {
                idx.push_back(std::vector<unsigned>(2,-1));
                continue;
            }
            
            // header
            if (i == 0)
            {
                std::vector<unsigned> tmp(2);
                const std::string& s = stateChars[0][j];
                std::string tok;
                std::stringstream ss(s.c_str());

                // branch end
                std::getline(ss, tok, '_');
                tmp[0] = (tok == "start") ? 0 : 1;
                
                // node
                std::getline(ss, tok, '_');
                tmp[1] = std::atoi(tok.c_str());
                
                idx.push_back(tmp);
            }
            
            // range-encoded states (e.g. 3 -> 1100, 1 -> 1000, etc.)
            else if (bgtype == "Range")
            {
                // get node, state, and bits for i,j cell in tab-delimited file
                size_t nide_idx = idx[j][1];
                unsigned stateIdx = 0;
                if (chartype == "NaturalNumbers")
                    stateIdx = std::atoi(stateChars[i][j].c_str());
                else if (chartype == "Standard")
                    stateIdx = (unsigned)standardStates.find( stateChars[i][j] );
                
                const std::vector<unsigned>& b = bits[stateIdx];
                
                // start
                if (idx[j][0] == 0) {
                    for (size_t k = 0; k < b.size(); k++)
                    {
                        marginalStartProbs[nide_idx][k] += b[k];
                    }
                }
                // end
                else if (idx[j][0] == 1) {
                    for (size_t k = 0; k < b.size(); k++)
                    {
                        marginalEndProbs[nide_idx][k] += b[k];
                    }
                }
            }
            
            // area-encoded states (e.g. 3 -> 0001, 1 -> 1000, etc.)
            else if (bgtype == "Area")
            {
                bool ambig = false;
                
                // get node, state, and bits for i,j cell in tab-delimited file
                size_t nide_idx = idx[j][1];
                unsigned stateIdx = 0;
                if (chartype == "NaturalNumbers")
                    stateIdx = std::atoi(stateChars[i][j].c_str());
                else if (chartype == "Standard")
                {
                    if (stateChars[i][j] == "-")
                        ambig = true;
                    stateIdx = (unsigned)standardStates.find( stateChars[i][j] );
                }

                // start
                if (idx[j][0] == 0) {
                    marginalStartProbs[nide_idx][stateIdx] += 1;
                }
                // end
                else if (idx[j][0] == 1) {
                    if (ambig) {
                        for (size_t jdx = 0; jdx < marginalEndProbs[nide_idx].size(); jdx++)
                        {
                            marginalEndProbs[nide_idx][jdx] += 1;
                        }
                    }
                    else {
                        marginalEndProbs[nide_idx][stateIdx] += 1;
                    }
                }
            }
        }
    }
    
    // normalize
    for (size_t i = 0; i < marginalEndProbs.size(); i++)
    {
        for (size_t j = 0; j < marginalEndProbs[i].size(); j++)
        {
            marginalStartProbs[i][j] = (double)(marginalStartProbs[i][j]) / numSamples;
            marginalEndProbs[i][j]   = (double)(marginalEndProbs[i][j])   / numSamples;
        }
    }
}

void PhylowoodConverter::makeBits(void)
{
    bits = std::vector<std::vector<unsigned> >(num_states, std::vector<unsigned>(numAreas, 0));
    for (size_t i = 1; i < num_states; i++)
    {
        size_t n = i;
        for (size_t j = 0; j < numAreas; j++)
        {
            bits[i][j] = n % 2;
            n /= 2;
            if (n == 0)
                break;
        }
    }
}

std::string PhylowoodConverter::buildCharacterHistoryString(const TopologyNode* n, unsigned end)
{
    std::stringstream ss;
    
    const size_t& nide_idx = n->getIndex();
    const std::vector<double>& characters = ( end == 0 ? marginalStartProbs[nide_idx] : marginalEndProbs[nide_idx] );
    for (size_t i = 0; i < characters.size(); i++)
    {
        if (i != 0)
            ss << ",";
        ss << characters[i];
    }
    
    return ss.str();
}


/* Build newick string */
std::string PhylowoodConverter::buildExtendedNewick(const TopologyNode* n ) {
    // create the newick string
    std::stringstream o;
    
    // extended data is only found on admixture nodes
    std::string additionalInfo = "";
    
    // loop over admixture nodes per branch
    std::stringstream characterStream;
    
    double br = 1.0;
    
    characterStream << "[";
    
    // character history
    characterStream << "&ch={" << buildCharacterHistoryString(n, 1) << "}";
    characterStream << ",&pa={" << buildCharacterHistoryString(n, 0) << "}";
    
    // ... whatever else
    characterStream << "]";
    
    additionalInfo = characterStream.str();
    
    // test whether this is a internal or external node
    if (n->isTip()) {
        // this is a tip so we just return the name of the node
        o << n->getIndex() << additionalInfo << ":" << tree->getBranchLengthForNode(*n);
    }
    
    else {
        o << "(";
        for (size_t i=0; i<(n->getNumberOfChildren()-1); i++) {
            o << buildExtendedNewick( &n->getChild(i) ) << ",";
        }
        o << buildExtendedNewick( &n->getChild(n->getNumberOfChildren()-1) ) << ")" << additionalInfo << ":" << tree->getBranchLengthForNode(*n) * br;
    }
    
    return o.str();
}

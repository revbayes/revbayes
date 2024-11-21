#include <cstdio>
#include <algorithm>
#include <string>
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "AminoAcidState.h"
#include "HomologousDiscreteCharacterData.h"
#include "DnaState.h"
#include "NewickTreeReader.h"
#include "NclReader.h"
#include "NonHomologousDiscreteCharacterData.h"
#include "RbFileManager.h"
#include "RnaState.h"
#include "StandardState.h"
#include "StringUtilities.h"
#include "TaxonMap.h"
#include "Tree.h"
#include "TopologyNode.h"
#include "TreeUtilities.h"
#include "RlUserInterface.h"
#include "AbstractCharacterData.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "DiscreteTaxonData.h"
#include "HomologousCharacterData.h"
#include "NaturalNumbersState.h"
#include "RbException.h"
#include "Taxon.h"
#include "nxsassumptionsblock.h"
#include "nxscharactersblock.h"
#include "nxsdefs.h"
#include "nxsexception.h"
#include "nxsmultiformat.h"
#include "nxsreader.h"
#include "nxsstring.h"
#include "nxstaxablock.h"
#include "nxstreesblock.h"
#include "nxsunalignedblock.h"

using namespace RevBayesCore;


NclReader::NclReader( void ) :
    nexusReader(-1, NxsReader::IGNORE_WARNINGS)
{
    
}


NclReader::NclReader( const NclReader &r ) :
    nexusReader(-1, NxsReader::IGNORE_WARNINGS)
{
    
}


//!< Check and change the indices of the trees to match
void NclReader::checkTreeTaxonIndices( std::vector<Tree* >* trees )
{
    
    // check for NULL pointer
    if ( trees != NULL && trees->size() > 1 )
    {
        
        Tree *reference_tree = (*trees)[0];
        TaxonMap tm = TaxonMap( *reference_tree );
        for (size_t i=1; i<trees->size(); ++i)
        {
            Tree *t = (*trees)[i];
            t->setTaxonIndices( tm );
        }
        
    }
    
}


/** Constructs a tree from NCL */
void NclReader::constructBranchLengthTreefromNclRecursively(TopologyNode* tn, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens, bool& reindex, const NxsSimpleNode* tnNcl, const NxsTaxaBlock *tb)
{
    
    // add the current node to the vector of nodes
    nodes.push_back( tn );
    
    // add the branch length for the node
    double edgeLength = tnNcl->GetEdgeToParent().GetDblEdgeLen();
    edgeLength = (edgeLength < 0 ? 0.0 : edgeLength);
    
    const std::vector<NxsComment>& br_comms = tnNcl->GetEdgeToParent().GetUnprocessedComments();
    int node_idx = -1;
    for (size_t i=0; i<br_comms.size(); ++i)
    {
        const NxsComment& curr_comm = br_comms[i];
        const std::string& comm_text = curr_comm.GetText();
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(comm_text, "=", tokens);
        if ( tokens[0] == "&index" )
        {
            node_idx = StringUtilities::asIntegerNumber( tokens[1] ) - 1;
            tn->setIndex( node_idx );
            break;
        }
        
        
    }
        
    // add the current node to the vector of nodes
    reindex = (reindex || node_idx == -1);
    
    
    //    tn->setBranchLength( edgeLength );
    // We only store the branch length because we can only set it later once the entire tree is build
    brlens.push_back( edgeLength );
    
    // get the children
    std::vector<NxsSimpleNode*> children = tnNcl->GetChildren();
    
    // iterate over all children
    for (std::vector<NxsSimpleNode*>::iterator it = children.begin(); it!=children.end(); it++)
    {
        // create a new tree node with given name
        std::string name = (*it)->GetName();
        if ( (*it)->GetTaxonIndex() < tb->GetNumTaxonLabels() )
        {
            name = tb->GetTaxonLabel( (*it)->GetTaxonIndex() ).BlanksToUnderscores();
        }
        TopologyNode* child = new TopologyNode(name);
        
        // add new node as child
        tn->addChild(child);
        
        // add parent to new node
        child->setParent(tn);
        
        // recursive call for the child to parse the tree
        constructBranchLengthTreefromNclRecursively(child, nodes, brlens, reindex, *it, tb);
    }
}


/** Reads the blocks stored by NCL and converts them to RevBayes character matrices */
std::vector<AbstractCharacterData* > NclReader::convertFromNcl(const path& file_name)
{
    
    std::vector<AbstractCharacterData* > cmv;
    
    size_t num_taxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (unsigned tBlck=0; tBlck<num_taxaBlocks; ++tBlck)
    {
        NxsTaxaBlock* taxaBlock = nexusReader.GetTaxaBlock(tBlck);
        std::string taxaBlockTitle          = taxaBlock->GetTitle();
        const unsigned nCharBlocks          = nexusReader.GetNumCharactersBlocks(taxaBlock);
        const unsigned nUnalignedCharBlocks = nexusReader.GetNumUnalignedBlocks(taxaBlock);
        
        // make alignment objects
        for (unsigned cBlck=0; cBlck<nCharBlocks; cBlck++)
        {
            AbstractCharacterData* m = NULL;
            NxsCharactersBlock* charBlock = nexusReader.GetCharactersBlock(taxaBlock, cBlck);
            std::string charBlockTitle = taxaBlock->GetTitle();
            int dt = charBlock->GetDataType();
            if (dt == NxsCharactersBlock::dna || dt == NxsCharactersBlock::nucleotide)
            {
                m = createDnaMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::rna)
            {
                m = createRnaMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::protein)
            {
                m = createAminoAcidMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::standard)
            {
                m = createStandardMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::continuous)
            {
                m = createContinuousMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::mixed)
            {
                addWarning("Mixed data types are not allowed");
            }
            else
            {
                addWarning("Unknown data type");
            }
            
            if (m != NULL)
            {
                
                
                
                m->setFilename( file_name );
                
                unsigned int nAssumptions = nexusReader.GetNumAssumptionsBlocks(charBlock);
                if ( nAssumptions > 0 )
                {
                    for (unsigned int i = 0; i < nAssumptions; ++i)
                    {
                        NxsAssumptionsBlock *assumption = nexusReader.GetAssumptionsBlock(charBlock,i);
                        size_t nSets = assumption->GetNumCharSets();
                        NxsStringVector names;
                        assumption->GetCharSetNames(names);
                        for (size_t j = 0; j < nSets; ++j)
                        {
                            const NxsUnsignedSet *set = assumption->GetCharSet(names[j]);
                            HomologousCharacterData *m_tmp = dynamic_cast<HomologousCharacterData *>(m)->clone();
                            m_tmp->excludeAllCharacters();
                            for (std::set<unsigned>::iterator k = set->begin(); k != set->end(); k++)
                            {
                                m_tmp->includeCharacter( *k );
                            }
                            m_tmp->removeExcludedCharacters();
                            cmv.push_back( m_tmp );
                            
                        }
                        
                    }
                }
                else
                {
                    cmv.push_back( m );
                }
                
            }
        }
        
        // create unaligned data objects
        for (unsigned cBlck=0; cBlck<nUnalignedCharBlocks; cBlck++)
        {
            AbstractCharacterData* m = NULL;
            NxsUnalignedBlock* charBlock = nexusReader.GetUnalignedBlock(taxaBlock, cBlck);
            std::string charBlockTitle = taxaBlock->GetTitle();
            int dt = charBlock->GetDataType();
            if (dt == NxsCharactersBlock::dna || dt == NxsCharactersBlock::nucleotide)
            {
                m = createUnalignedDnaMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::rna)
            {
                m = createUnalignedRnaMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::protein)
            {
                m = createUnalignedAminoAcidMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::standard)
            {
                //m = createStandardMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::continuous)
            {
                //m = createContinuousMatrix(charBlock);
            }
            else if (dt == NxsCharactersBlock::mixed)
            {
                addWarning("Mixed data types are not allowed");
            }
            else
            {
                addWarning("Unknown data type");
            }
            
            if (m != NULL)
            {
                m->setFilename( file_name );
                cmv.push_back( m );
            }
        }
        
    }
    
    return cmv;
}


/** Converts trees stored by NCL into RevBayes formatted trees */
std::vector<Tree*>* NclReader::convertTreesFromNcl(void)
{
        
    const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    std::vector<Tree*>* rbTreesFromFile = new std::vector<Tree*>();
    for (unsigned t = 0; t < nTaxaBlocks; ++t)
    {
        const NxsTaxaBlock *tb = nexusReader.GetTaxaBlock(t);
        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
        if (nTreesBlocks == 0)
        {
            continue;
        }
        
        for (unsigned i = 0; i < nTreesBlocks; ++i)
        {
            const NxsTreesBlock * trb = nexusReader.GetTreesBlock(tb, i);
            trb->ProcessAllTrees();
            for (unsigned j = 0; j < trb->GetNumTrees(); ++j)
            {
                const NxsFullTreeDescription & ftd = trb->GetFullTreeDescription(j);
                NxsSimpleTree tree(ftd, -1, -1.0);

                Tree* rbTree = translateNclSimpleTreeToBranchLengthTree(tree,tb,ftd.IsRooted());


                rbTreesFromFile->push_back( rbTree );
            }
        }
    }
    
    return rbTreesFromFile;
}

/** Create an object to hold aligned amino acid data */
HomologousDiscreteCharacterData<AminoAcidState>* NclReader::createAminoAcidMatrix(NxsCharactersBlock* charblock)
{
    
    if ( charblock == NULL )
    {
        throw RbException("Trying to create an AA matrix from a NULL pointer.");
    }
    
    // check that the character block is of the correct type
    if (charblock->GetDataType() != NxsCharactersBlock::protein)
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    NxsUnsignedSet charset;
    for (unsigned int i=0; i<charblock->GetNumChar(); ++i)
    {
        charset.insert(i);
    }
    size_t numOrigTaxa = charblock->GetNTax();
    
    // get the set of excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    
    // instantiate the character matrix
    HomologousDiscreteCharacterData<AminoAcidState>* cMat = new HomologousDiscreteCharacterData<AminoAcidState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // get the taxon name
        NxsString   tLabel = charblock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of Standard states
        DiscreteTaxonData<AminoAcidState> dataVec = DiscreteTaxonData<AminoAcidState>(tokens[0]);
        
        for (NxsUnsignedSet::const_iterator cit = charset.begin(); cit != charset.end();cit++)
        {
            AminoAcidState aaState;
            if (charblock->IsGapState(origTaxIndex, *cit) == true)
            {
                aaState.setGapState(true);
//                aaState.setState("-");
            }
            else if ( charblock->IsMissingState(origTaxIndex, *cit) == true || charblock->GetNumStates(origTaxIndex, *cit) >= 20 )
            {
//                aaState.setState("?");
                aaState.setMissingState(true);
            }
            else
            {
                size_t nStates = charblock->GetNumStates(origTaxIndex, *cit);
                aaState.setState( std::string(1, charblock->GetState(origTaxIndex, *cit, 0) ) );
                for (unsigned s=1; s<nStates; s++)
                {
                    aaState.addState( std::string(1, charblock->GetState(origTaxIndex, *cit, s) ) );
                }
            }
            dataVec.addCharacter( aaState );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    setExcluded( charblock, cMat );
    
    return cMat;
}

/* Create an object to hold unaligned aminoacid data */
NonHomologousDiscreteCharacterData<AminoAcidState>* NclReader::createUnalignedAminoAcidMatrix(NxsUnalignedBlock* charblock)
{
    
    // check that the character block is of the correct type
    if ( charblock->GetDataType() != NxsCharactersBlock::protein )
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    unsigned numOrigTaxa = charblock->GetNTax();
    NxsTaxaBlockAPI* taxonBlock = charblock->GetTaxaBlockPtr();
    
    // instantiate the character matrix
    NonHomologousDiscreteCharacterData<AminoAcidState>* cMat = new NonHomologousDiscreteCharacterData<AminoAcidState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = taxonBlock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of amino acid states
        DiscreteTaxonData<AminoAcidState> dataVec = DiscreteTaxonData<AminoAcidState>(tokens[0]);
        
        // add the sequence information for the sequence associated with the taxon
        std::string rowDataAsString = charblock->GetMatrixRowAsStr(origTaxIndex);
        for (size_t i=0; i<rowDataAsString.size(); i++)
        {
            AminoAcidState aaState;
            aaState.setState( std::string(1, rowDataAsString[i] ) );
            dataVec.addCharacter( aaState );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    return cMat;
}

/** Create an object to hold continuous data */
ContinuousCharacterData* NclReader::createContinuousMatrix(NxsCharactersBlock* charblock)
{
    
    // check that the character block is of the correct type
    if (charblock->GetDataType() != NxsCharactersBlock::continuous)
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    NxsUnsignedSet charset;
    for (unsigned int i=0; i<charblock->GetNumChar(); i++)
        charset.insert(i);
    unsigned numOrigTaxa = charblock->GetNTax();
    
    // get the set of excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    
    // instantiate the character matrix
    ContinuousCharacterData* cMat = new ContinuousCharacterData();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = charblock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of Standard states
        ContinuousTaxonData dataVec = ContinuousTaxonData(tokens[0]);
        
        // add the real-valued observation
        for (NxsUnsignedSet::const_iterator cit = charset.begin(); cit != charset.end();cit++)
        {
            double contObs = RbConstants::Double::nan;
            bool is_resolved = false;
//            if ( charblock->IsMissingState(origTaxIndex, *cit) == false )
//            {
                const std::vector<double>& x = charblock->GetContinuousValues( origTaxIndex, *cit, "AVERAGE" );
                if ( x.size() > 0 )
                {
                    contObs = x[0];
                    is_resolved = true;
                }
//            }
            dataVec.addCharacter( contObs, is_resolved );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    setExcluded( charblock, cMat );
    
    return cMat;
}


/** Create an object to hold aligned DNA data */
HomologousDiscreteCharacterData<DnaState>* NclReader::createDnaMatrix(NxsCharactersBlock* charblock)
{
    
    if ( charblock == NULL )
    {
        throw RbException("Trying to create an DNA matrix from a NULL pointer.");
    }
    
    // check that the character block is of the correct type
    if ( charblock->GetDataType() != NxsCharactersBlock::dna )
    {
        std::cout << "Could not read in data matrix of type DNA because the nexus files says the type is:" << std::endl;
        switch ( charblock->GetDataType() )
        {
        case 1:
            std::cout << "Standard" << std::endl;
            break;
                
        case 2:
            std::cout << "DNA" << std::endl;
            break;
                
        case 3:
            std::cout << "RNA" << std::endl;
            break;
                
        case 4:
            std::cout << "Nucleotide" << std::endl;
            break;
                
        case 5:
            std::cout << "Protein" << std::endl;
            break;
                
        case 6:
            std::cout << "Continuous" << std::endl;
            break;
                
        case 7:
            std::cout << "Codon" << std::endl;
            break;
                
        case 8:
            std::cout << "Mixed" << std::endl;
            break;
                
        default:
            std::cout << "Unknown" << std::endl;
            break;
        }
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    NxsUnsignedSet charset;
    for (unsigned int i=0; i<charblock->GetNumChar(); i++)
    {
        charset.insert(i);
    }
    
    unsigned numOrigTaxa = charblock->GetNTax();
    
    // get the set of excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    
    // instantiate the character matrix
    HomologousDiscreteCharacterData<DnaState>* cMat = new HomologousDiscreteCharacterData<DnaState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = charblock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of DNA states
        DiscreteTaxonData<DnaState> dataVec = DiscreteTaxonData<DnaState>( tokens[0] );
        
        // add the sequence information for the sequence associated with the taxon
        for (NxsUnsignedSet::iterator cit = charset.begin(); cit != charset.end(); cit++)
        {
            // add the character state to the matrix
            DnaState dnaState;
            bool isResolved = true;
            if ( charblock->IsGapState(origTaxIndex, *cit) == true )
            {
                dnaState.setState("-");
                dnaState.setGapState(true);
                isResolved = false;
            }
            else if (charblock->IsMissingState(origTaxIndex, *cit) == true)
            {
                dnaState.setState("?");
                dnaState.setMissingState(true);
                isResolved = false;
            }
            else
            {
                dnaState.setState( std::string(1, charblock->GetState(origTaxIndex, *cit, 0) ) );
                for (unsigned int s=1; s<charblock->GetNumStates(origTaxIndex, *cit); s++)
                {
                    dnaState.addState( std::string(1, charblock->GetState(origTaxIndex, *cit, s) ) );
                }
            }
            dataVec.addCharacter(dnaState, isResolved);
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    setExcluded( charblock, cMat );
    
    return cMat;
}

/** Create an object to hold unaligned DNA data */
NonHomologousDiscreteCharacterData<DnaState>* NclReader::createUnalignedDnaMatrix(NxsUnalignedBlock* charblock)
{
    
    // check that the character block is of the correct type
    if ( charblock->GetDataType() != NxsCharactersBlock::dna )
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    unsigned numOrigTaxa = charblock->GetNTax();
    NxsTaxaBlockAPI* taxonBlock = charblock->GetTaxaBlockPtr();
    
    // instantiate the character matrix
    NonHomologousDiscreteCharacterData<DnaState>* cMat = new NonHomologousDiscreteCharacterData<DnaState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = taxonBlock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of DNA states
        DiscreteTaxonData<DnaState> dataVec = DiscreteTaxonData<DnaState>(tokens[0]);
        
        // add the sequence information for the sequence associated with the taxon
        std::string rowDataAsString = charblock->GetMatrixRowAsStr(origTaxIndex);
        for (size_t i=0; i<rowDataAsString.size(); i++)
        {
            DnaState dnaState;
            dnaState.setState( std::string(1, rowDataAsString[i] ) );
            dataVec.addCharacter( dnaState );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    return cMat;
}

/** Create an object to hold aligned RNA data */
HomologousDiscreteCharacterData<RnaState>* NclReader::createRnaMatrix(NxsCharactersBlock* charblock) {
    
    if ( charblock == NULL )
    {
        throw RbException("Trying to create an RNA matrix from a NULL pointer.");
    }
    
    // check that the character block is of the correct type
    if ( charblock->GetDataType() != NxsCharactersBlock::rna )
    {
        std::cout << "Could not read in data matrix of type RNA because the nexus files says the type is:" << std::endl;
        switch ( charblock->GetDataType() )
        {
        case 1:
            std::cout << "Standard" << std::endl;
            break;
                
        case 2:
            std::cout << "DNA" << std::endl;
            break;
                
        case 3:
            std::cout << "RNA" << std::endl;
            break;
                
        case 4:
            std::cout << "Nucleotide" << std::endl;
            break;
                
        case 5:
            std::cout << "Protein" << std::endl;
            break;
                
        case 6:
            std::cout << "Continuous" << std::endl;
            break;
                
        case 7:
            std::cout << "Codon" << std::endl;
            break;
                
        case 8:
            std::cout << "Mixed" << std::endl;
            break;
                
        default:
            std::cout << "Unknown" << std::endl;
            break;
        }
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    NxsUnsignedSet charset;
    for (unsigned int i=0; i<charblock->GetNumChar(); i++)
    {
        charset.insert(i);
    }
    
    unsigned numOrigTaxa = charblock->GetNTax();
    
    // get the set of excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    
    // instantiate the character matrix
    HomologousDiscreteCharacterData<RnaState>* cMat = new HomologousDiscreteCharacterData<RnaState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = charblock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of DNA states
        DiscreteTaxonData<RnaState> dataVec = DiscreteTaxonData<RnaState>( tokens[0] );
        
        // add the sequence information for the sequence associated with the taxon
        for (NxsUnsignedSet::iterator cit = charset.begin(); cit != charset.end(); cit++)
        {
            // add the character state to the matrix
            RnaState rnaState;
            bool isResolved = true;
            if ( charblock->IsGapState(origTaxIndex, *cit) == true )
            {
                rnaState.setState("-");
                rnaState.setGapState(true);
                isResolved = false;
            }
            else if (charblock->IsMissingState(origTaxIndex, *cit) == true)
            {
                rnaState.setState("?");
                rnaState.setMissingState(true);
                isResolved = false;
            }
            else
            {
                rnaState.setState( std::string(1, charblock->GetState(origTaxIndex, *cit, 0) ) );
                for (unsigned int s=1; s<charblock->GetNumStates(origTaxIndex, *cit); s++)
                {
                    rnaState.addState( std::string(1, charblock->GetState(origTaxIndex, *cit, s) ) );
                }
            }
            dataVec.addCharacter(rnaState, isResolved);
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    setExcluded( charblock, cMat );
    
    return cMat;
}

/** Create an object to hold unaligned RNA data */
NonHomologousDiscreteCharacterData<RnaState>* NclReader::createUnalignedRnaMatrix(NxsUnalignedBlock* charblock)
{
    
    // check that the character block is of the correct type
    if ( charblock->GetDataType() != NxsCharactersBlock::rna )
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    unsigned numOrigTaxa = charblock->GetNTax();
    NxsTaxaBlockAPI* taxonBlock = charblock->GetTaxaBlockPtr();
    
    // instantiate the character matrix
    NonHomologousDiscreteCharacterData<RnaState>* cMat = new NonHomologousDiscreteCharacterData<RnaState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = taxonBlock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of DNA states
        DiscreteTaxonData<RnaState> dataVec = DiscreteTaxonData<RnaState>(tokens[0]);
        
        // add the sequence information for the sequence associated with the taxon
        std::string rowDataAsString = charblock->GetMatrixRowAsStr(origTaxIndex);
        for (size_t i=0; i<rowDataAsString.size(); i++)
        {
            RnaState rnaState;
            rnaState.setState( std::string(1, rowDataAsString[i] ) );
            dataVec.addCharacter( rnaState );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    return cMat;
}

/** Create an object to hold standard data */
HomologousDiscreteCharacterData<StandardState>* NclReader::createStandardMatrix(NxsCharactersBlock* charblock)
{
    
    // check that the character block is of the correct type
    if (charblock->GetDataType() != NxsCharactersBlock::standard)
    {
        return NULL;
    }
    
    // get the set of characters (and the number of taxa)
    NxsUnsignedSet charset;
    for (unsigned int i=0; i<charblock->GetNumChar(); i++)
    {
        charset.insert(i);
    }
    unsigned numOrigTaxa = charblock->GetNTax();
    
    // get the set of excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    
    // get the number of states
    const NxsDiscreteDatatypeMapper* mapper = charblock->GetDatatypeMapperForChar(0);
    std::string sym = charblock->GetSymbols();
    size_t nStates = mapper->GetNumStates();
    if (nStates > 100)
    {
        return NULL;
    }
    
    // instantiate the character matrix
    HomologousDiscreteCharacterData<StandardState>* cMat = new HomologousDiscreteCharacterData<StandardState>();
    
    // read in the data, including taxon names
    for (unsigned origTaxIndex=0; origTaxIndex<numOrigTaxa; origTaxIndex++)
    {
        // add the taxon name
        NxsString   tLabel = charblock->GetTaxonLabel(origTaxIndex);
        std::string tName  = NxsString::GetEscaped(tLabel).c_str();
        
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(tName, "|", tokens);
        
        // allocate a vector of Standard states
        DiscreteTaxonData<StandardState> dataVec = DiscreteTaxonData<StandardState>(tokens[0]);

        // add the character information for the data associated with the taxon
        for (NxsUnsignedSet::iterator cit = charset.begin(); cit != charset.end(); cit++)
        {
            // add the character state to the matrix
            StandardState stdState = StandardState(sym.substr(0,1),sym);
            if ( charblock->IsGapState(origTaxIndex, *cit) == true )
            {
                stdState.setGapState(true);
            }
            else if (charblock->IsMissingState(origTaxIndex, *cit) == true)
            {
                stdState.setMissingState(true);
            }
            else
            {
                for (unsigned int s=0; s<charblock->GetNumStates(origTaxIndex, *cit); s++)
                {
                    stdState.setState( std::string(1, charblock->GetState(origTaxIndex, *cit, 0) ) );
                    for (unsigned int s=1; s<charblock->GetNumStates(origTaxIndex, *cit); s++)
                    {
                        stdState.addState( std::string(1, charblock->GetState(origTaxIndex, *cit, s) ) );
                    }
                }
            }
            dataVec.addCharacter( stdState );
        }
        
        // add sequence to character matrix
        cMat->addTaxonData( dataVec );
    }
    
    setExcluded( charblock, cMat );
    
    return cMat;
}


void NclReader::getTranslateTables(std::vector<std::map<int,std::string> >& translateTables)
{

    const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (unsigned t=0; t<nTaxaBlocks; t++)
    {
        const NxsTaxaBlock* tb = nexusReader.GetTaxaBlock(t);
        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
        if (nTreesBlocks == 0)
            return;
        
        std::map<int,std::string> newTable;
        for (unsigned i=0; i<tb->GetNTax(); i++)
        {
            NxsString ts = tb->GetTaxonLabel(i);
            newTable.insert( make_pair(i, ts) );
        }
        
        translateTables.push_back(newTable);
    }
    
}

/** Attempt to determine the type of data this is being read */
std::string NclReader::intuitDataType(std::string& s)
{
    
    // set up some strings containing the valid states for the different data types we are interested
    // in distiguishing
    static std::string dnaStates = "acgtmgrsvwyhkdbn-.?";
    static std::string rnaStates = "acgumgrsvwyhkdbn-.?";
    static std::string aaStates  = "arndcqeghilkmfpstwxyv-.?";
    static std::string stdStates = "0123456789n-.?abcdefghijklmnopqrstuvwxyz()";
    static std::string nucStates = "acgtu";
    
    // and intialize a few variables we'll be needing
    size_t nucCount = 0, nMissing = 0;
    bool notDna = false, notRna = false, notAa = false, notStd = false;
    
    // loop over the string (s) that contains the raw data we look at the state and try to determine if the
    // state rules out certain data types
    StringUtilities::toLower( s );
    s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );
    
    for (size_t i=0; i<s.size(); i++)
    {
        char c = s[i];
        
        if (c == 'n' || c == '-' || c == '.' || c == '?')
            nMissing++;
        
        if (notDna == false)
        {
            bool foundState = false;
            for (size_t j=0; j<19; j++)
            {
                if ( c == dnaStates[j] )
                {
                    foundState = true;
                    break;
                }
            }
            if (foundState == false)
            {
                notDna = true;
            }
            
        }
        
        if (notRna == false)
        {
            bool foundState = false;
            for (size_t j=0; j<19; j++)
            {
                if ( c == rnaStates[j] )
                {
                    foundState = true;
                    break;
                }
            }
            if (foundState == false)
            {
                notRna = true;
            }
            
        }
        
        if (notAa == false)
        {
            bool foundState = false;
            for (size_t j=0; j<23; j++)
            {
                if ( c == aaStates[j] )
                {
                    foundState = true;
                    break;
                }
            }
            if (foundState == false)
            {
                notAa = true;
            }
            
        }
        
        if (notStd == false)
        {
            bool foundState = false;
            for (size_t j=0; j<14; j++)
            {
                if ( c == stdStates[j] )
                {
                    foundState = true;
                    break;
                }
            }
            if (foundState == false)
            {
                notStd = true;
            }
            
        }
        
        for (size_t j=0; j<5; j++)
        {
            if ( c == nucStates[j] )
                nucCount++;
        }
    }
    
    // Try to determine which type of data it is. Note that using only the sequence information in the
    // string s that it is difficult or impossible to rule out certain data types. For example, if the
    // data is DNA, then the DNA states with the ambiguity codes do not rule out amino acid data. However,
    // we would expect a lot of A, C, G, and Ts in our data and not much else if the data is DNA, so we
    // use a bit of "fuzzy" logic to determine if the data is highly likely to be nucleotide, or not.
    if (notDna == false && notRna == true && notAa == true && notStd == true)
        return "dna";
    else if (notDna == true && notRna == false && notAa == true && notStd == true)
        return "rna";
    else if (notDna == true && notRna == true && notAa == false && notStd == true)
        return "protein";
    else if (notDna == true && notRna == true && notAa == true && notStd == false)
        return "standard";
    else if (notDna == false && notRna == true && notAa == false && notStd == true)
    {
        if ( (double)nucCount / (s.size()-nMissing) > 0.8 )
            return "dna";
        else
            return "protein";
    }
    else if (notDna == true && notRna == false && notAa == false && notStd == true)
    {
        if ( (double)nucCount / (s.size()-nMissing) > 0.8 )
            return "rna";
        else
            return "protein";
    }
    else if ( notDna == false )
    {
        return "dna";
    }

    return "";
}


/** Try to determine if the file is likely to be in Fasta format */
bool NclReader::isFastaFile(const path& fn, std::string& dType)
{
    
    // open file
    std::ifstream fStrm( fn.string() );
    
    // read the file token-by-token looking for Fasta things
    int ch = fStrm.peek();
    std::string word = "";
    std::string seqStr = "";
    int wordNum = 0, lineNum = 0, lastCarotLine = -100;
    int numSequences = 0;
    while (ch != EOF)
    {
        safeGetline(fStrm, word);
        
        // we know that the last character is an escape character
        if ( word.size() > 0 )
        {
            word.erase(word.size()-1);
        }
        if (wordNum == 0 && word[0] == '>')
        {
            if (lineNum - lastCarotLine > 1)
            {
                lastCarotLine = lineNum;
                numSequences++;
            }
            else
            {
                return false;
            }
        }
        else if (wordNum == 0 && word[0] == ';')
        {
            // comment
        }
        else if (lineNum > 0 && word[0] != '>' && word[0] != ';')
        {
            seqStr += word;
        }
        
        wordNum++;
        ch = fStrm.peek();

        lineNum++;
        wordNum = 0;
    }
    
    // close file
    fStrm.close();
    
    if (numSequences < 1)
        return false;
    
    dType = intuitDataType(seqStr);
    
    return true;
}


/**
 * Try to determine if the file is likely to be in Nexus format. We check if the first word in the
 * file is #NEXUS. If not, then we check if the file name ending is ".nex". If neither is true, it
 * is probably not a NEXUS file.
 */
bool NclReader::isNexusFile(const path& fn)
{
    
    // open file, read first word, close file
    std::ifstream fStrm( fn.string() );

    std::string word;
    fStrm >> word;
    fStrm.close();
    
    if (word=="#NEXUS")
    {
        return true;
    }
    else
    {
        if ( fn.extension() == ".nex" )
            return true;
        else
            return false;
    }
}


/** Try to determine if the file is likely to be in Phylip format */
bool NclReader::isPhylipFile(const path& fn, std::string& dType, bool& is_interleaved)
{
    // open file
    std::ifstream fStrm( fn.string() );

    std::string seqStr = "";
    
    // read the file token-by-token looking for NEXUS things
    bool foundNumTaxa = false, foundNumChar = false;
    unsigned int num_taxa = 0;
    std::vector<std::string> taxonNames;
    int ch = fStrm.get();
    fStrm.unget();
    std::string word = "";
    int wordNum = 0, lineNum = 0;
    while (ch != EOF)
    {
        word = "";
        fStrm >> word;
        StringUtilities::toLower( word );
        
        if (lineNum == 0 && wordNum == 0 && StringUtilities::isNumber(word) == true)
        {
            std::istringstream buf(word);
            buf >> num_taxa;
            foundNumTaxa = true;
        }
        else if (lineNum == 0 && wordNum == 1 && StringUtilities::isNumber(word) == true)
            foundNumChar = true;
        else if (lineNum > 0 && wordNum == 0 && word != "")
            taxonNames.push_back( word );
        else if (lineNum > 0 && wordNum > 0)
            seqStr += word;
        
        wordNum++;
        ch = fStrm.get();
        if (ch == '\n' || ch == '\r' || ch == EOF)
        {
            lineNum++;
            wordNum = 0;
        }
    }
    
    // close file
    fStrm.close();
    
    if (foundNumTaxa == true && foundNumChar == true)
    {
        if (taxonNames.size() == 0)
            return false;
        if (taxonNames.size() % num_taxa != 0)
            return false;
        
        if (taxonNames.size() > num_taxa)
            is_interleaved = true;
        dType = intuitDataType(seqStr);
        return true;
    }
    
    return false;
}


std::vector<AbstractCharacterData*> NclReader::readMatrices(const path &fn)
{
    
    // check that the file/path name has been correctly specified
    if ( fn.empty() )
    {
        std::string errorStr = "";
        formatError(fn, errorStr);
        throw RbException()<<"Could not find file or path with name "<<fn;
    }
    
    // are we reading a single file or are we reading the contents of a directory?
    bool readingDirectory = false;
    if ( fn.filename().empty() or fn.filename_is_dot() or fn.filename_is_dot_dot())
        readingDirectory = true;
    if (readingDirectory == true)
        RBOUT( "Recursively reading the contents of a directory\n" );
    else
        RBOUT( "Attempting to read the contents of file \"" + fn.filename().string() + "\"\n" );
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<path> vectorOffile_names;
    if (readingDirectory == true)
        setStringWithNamesOfFilesInDirectory(fn, vectorOffile_names);
    else
    {
        vectorOffile_names.push_back( fn );
    }

    if (readingDirectory == true)
    {
        std::stringstream o1;
        o1 << "Found " << vectorOffile_names.size() << " files in directory";
        RBOUT( o1.str() );
    }
    
    // clear warnings from warnings buffer
    clearWarnings();
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOffile_names" because some of them may not be in a format
    // that can be read.
    std::map<path,std::string> fileMap;
    for (std::vector<path>::iterator p = vectorOffile_names.begin(); p != vectorOffile_names.end(); p++)
    {
        bool is_interleaved = false;
        std::string myFileType = "unknown", dType = "unknown";
        if (isNexusFile(*p) == true)
            myFileType = "nexus";
        else if (isPhylipFile(*p, dType, is_interleaved) == true)
            myFileType = "phylip";
        else if (isFastaFile(*p, dType) == true)
            myFileType = "fasta";
        
        if (myFileType != "unknown")
        {
            std::string suffix = "|" + dType;
            if ( myFileType == "phylip" )
            {
                if (is_interleaved == true)
                    suffix += "|interleaved";
                else
                    suffix += "|noninterleaved";
            }
            else if ( myFileType == "fasta" )
                suffix += "|noninterleaved";
            else
                suffix += "|unknown";
            myFileType += suffix;
            fileMap.insert( std::make_pair(*p,myFileType) );
        }
        else
        {
            addWarning("Unknown file type");
        }
    }
    
    // read the files in the map containing the file names with the output being a vector of pointers to
    // the character matrices that have been read
    std::vector<AbstractCharacterData*> m = readMatrices( fileMap );
    
    // print summary of results of file reading to the user
    if (readingDirectory == true)
    {
        std::stringstream o2;
        if ( m.size() == 0 )
            o2 << "Failed to read any files";
        else if ( m.size() == 1 )
            o2 << "Successfully read one file";
        else
            o2 << "Successfully read " << m.size() << " files";
        RBOUT( o2.str() );
        std::set<std::string> myWarnings = getWarnings();
        if ( vectorOffile_names.size() - m.size() > 0 && myWarnings.size() > 0 )
        {
            std::stringstream o3;
            if (vectorOffile_names.size() - m.size() == 1)
                o3 << "Did not read a file for the following ";
            else
                o3 << "Did not read " << vectorOffile_names.size() - m.size() << " files for the following ";
            if (myWarnings.size() == 1)
                o3 << "reason:";
            else
                o3 << "reasons:";
            RBOUT( o3.str() );
            for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
                RBOUT( "* "+(*it) );
        }
    }
    else
    {
        if (m.size() > 0)
            RBOUT( "Successfully read file" );
        else
        {
            std::set<std::string> myWarnings = getWarnings();
            if ( myWarnings.size() > 0 )
            {
                std::stringstream o3;
                o3 << "Did not read the file for the following ";
                if (myWarnings.size() == 1)
                    o3 << "reason:";
                else
                    o3 << "reasons:";
                RBOUT( o3.str() );
                for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
                    RBOUT( "* "+(*it) );
            }
        }
    }
    
    return m;
}

/** Read a list of file names contained in a map (with file format info too) */
std::vector<AbstractCharacterData *> NclReader::readMatrices(const path &fn, const std::string &ft)
{
    // instantiate a vector of matrices
    std::vector<AbstractCharacterData* > cmv;
    
    // The data types are as follows: Nexus, Phylip+DNA/RNA/AA/Standard+Interleaved/Noninterleaved,
    // Fasta+DNA/RNA/AA.
    
    // Check that the file exists. It is likely that this has been already checked during the formation of
    // the map that is passed into the function, but it never hurts to be safe...
    if ( not exists( fn ))
    {
        addWarning("Data file not found");
    }
    else
    {
        // Extract information on the file format from the value of the key/value pair. Note that we expect the
        // file_format string to be in the format file_type|data_type|interleave_type with pipes ('|') separating
        // the format components. It might be better to make an object value in the key/value pair that contains
        // this information.
        std::vector<std::string> file_format;
        StringUtilities::stringSplit( ft, "|", file_format );
        std::string ff = file_format[0];
        std::string dt = file_format[1];
        bool il = false;
        if ( file_format[2] == "interleaved" )
        {
            il = true;
        }
        
        // read the file
        cmv = readMatrices( fn, ff, dt, il );
        nexusReader.ClearContent();
    }
    
    return cmv;
}


/** Read a list of file names contained in a map (with file format info too) */
std::vector<AbstractCharacterData *> NclReader::readMatrices(const std::map<path,std::string>& fileMap) {
    
    // instantiate a vector of matrices
    std::vector<AbstractCharacterData* > cmv;
    
    // We loop through the map that contains as the key the path and name of the file to be read and as the
    // value the data type. The data types are as follows: Nexus, Phylip+DNA/RNA/AA/Standard+Interleaved/Noninterleaved,
    // Fasta+DNA/RNA/AA.
    for (auto p = fileMap.begin(); p != fileMap.end(); p++)
    {
        // Check that the file exists. It is likely that this has been already checked during the formation of
        // the map that is passed into the function, but it never hurts to be safe...
        if ( not exists(p->first) )
        {
            addWarning("Data file not found");
            continue;
        }
        
        // Extract information on the file format from the value of the key/value pair. Note that we expect the
        // file_format string to be in the format file_type|data_type|interleave_type with pipes ('|') separating
        // the format components. It might be better to make an object value in the key/value pair that contains
        // this information.
        std::vector<std::string> file_format;
        StringUtilities::stringSplit( p->second, "|", file_format );
        std::string ff = file_format[0];
        std::string dt = file_format[1];
        bool il = false;
        if ( file_format[2] == "interleaved" )
            il = true;
        
        // read the file
        std::vector<AbstractCharacterData* > v = readMatrices( p->first, ff, dt, il );
        if (v.size() > 0)
        {
            for (std::vector<AbstractCharacterData*>::iterator m = v.begin(); m != v.end(); m++)
                cmv.push_back( (*m) );
        }
        nexusReader.ClearContent();
    }
    
    return cmv;
}


/** Read a list of file names contained in a vector of strings */
std::vector<AbstractCharacterData*> NclReader::readMatrices(const std::vector<path> fn, const std::string file_format, const std::string data_type, const bool is_interleaved) {
    
    // instantiate a vector of matrices
    std::vector<AbstractCharacterData*> cmv;
    
    // check that the files exist
    for (auto f = fn.begin(); f != fn.end(); f++)
    {
        if ( not exists(*f) )
        {
            addWarning("Data file not found");
            return cmv;
        }
    }
    
    // read the data files
    for (auto f = fn.begin(); f != fn.end(); f++)
    {
        std::vector<AbstractCharacterData*> v = readMatrices( (*f), file_format, data_type, is_interleaved );
        if (v.size() > 0)
        {
            for (std::vector<AbstractCharacterData*>::iterator m = v.begin(); m != v.end(); m++)
                cmv.push_back( (*m) );
        }
        nexusReader.ClearContent();
    }
    
    return cmv;
}


/** Reads a single file using NCL */
std::vector<AbstractCharacterData*> NclReader::readMatrices(const path &file_name, const std::string &file_format, const std::string &data_type, bool is_interleaved)
{
    
    std::string fns = file_name.string();
    
    const char *fn = fns.c_str();
    
    // check that the file exists
    if ( not exists(file_name) )
    {
        addWarning("Data file not found");
        std::vector<AbstractCharacterData*> dummy;
        return dummy;
    }
    
    try
    {
        if (file_format == "nexus")
        {
            // NEXUS file format
            nexusReader.ReadFilepath(fn, MultiFormatReader::NEXUS_FORMAT);
        }
        else if (file_format == "fasta")
        {
            // fasta file format
            if (data_type == "dna")
            {
                nexusReader.ReadFilepath(fn, MultiFormatReader::FASTA_DNA_FORMAT);
            }
            else if (data_type == "rna")
            {
                nexusReader.ReadFilepath(fn, MultiFormatReader::FASTA_RNA_FORMAT);
            }
            else if (data_type == "protein")
            {
                nexusReader.ReadFilepath(fn, MultiFormatReader::FASTA_AA_FORMAT);
            }
            else
            {
                throw RbException() << "Unknown data type '" <<  data_type << "' for fasta formatted files.";
            }
        }
        else if (file_format == "phylip")
        {
            // phylip file format
            if (is_interleaved == false)
            {
                if (data_type == "dna")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT);
                }
                else if (data_type == "rna")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT);
                }
                else if (data_type == "protein")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT);
                }
                else if (data_type == "standard")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT);
                }
                
            }
            else
            {
                if (data_type == "dna")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT);
                }
                else if (data_type == "rna")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT);
                }
                else if (data_type == "protein")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT);
                }
                else if (data_type == "standard")
                {
                    nexusReader.ReadFilepath(fn, MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT);
                }
                
            }
            
        }
        
    }
    catch(NxsException& err)
    {
        if ( err.msg.length() == 0 )
        {
            // Basic error message if ncl fails to give something back
            addWarning("Nexus error in file \"" + file_name.filename().string() + "\"");
        }
        else
        {
            // NxsReader error message
            addWarning(err.msg);
        }
        
        // Position information
        std::stringstream errorMessage;
        errorMessage << "The error occurred while reading line ";
        errorMessage << err.line << " column " << err.col;
        errorMessage << " in file " << file_name.filename();
        addWarning(errorMessage.str());
        
        // Return empty character matrix vector
        std::vector<AbstractCharacterData*> dummy;
        return dummy;
    }
    
    std::vector<AbstractCharacterData*> cvm = convertFromNcl(file_name);
    return cvm;
}


/** Read trees */
std::vector<Tree*>* NclReader::readBranchLengthTrees(const path &fn)
{
    
    nexusReader.ClearContent();
    
    // check that the file/path name has been correctly specified
    if ( fn.empty() )
    {
        std::string errorStr = "";
        formatError(fn, errorStr);
        throw RbException()<<"Could not find file or path with name "<<fn;
    }
    
    // are we reading a single file or are we reading the contents of a directory?
    bool readingDirectory = false;
    if ( fn.filename().empty() or fn.filename_is_dot() or fn.filename_is_dot_dot())
    {
        readingDirectory = true;
    }
    if (readingDirectory == true)
    {
        RBOUT( "Recursively reading the contents of a directory" );
    }
    else
    {
        RBOUT( "Attempting to read the contents of file \"" + fn.filename().string() + "\"" );
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<path> vectorOffile_names;
    if (readingDirectory == true)
    {
        setStringWithNamesOfFilesInDirectory(fn, vectorOffile_names);
    }
    else
    {
        vectorOffile_names.push_back( fn );
    }
    if (readingDirectory == true)
    {
        std::stringstream o1;
        o1 << "Found " << vectorOffile_names.size() << " files in directory";
        RBOUT( o1.str() );
    }
    
    // clear warnings from its warnings buffer
    clearWarnings();
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOffile_names" because some of them may not be in a format
    // that can be read.
    std::vector<Tree*> *trees = NULL;
    for (auto p = vectorOffile_names.begin(); p != vectorOffile_names.end(); p++)
    {
        // we should check here the file type first and make sure it is valid
        
        // read the files in the map containing the file names with the output being a vector of pointers to
        // the character matrices that have been read
        trees = readBranchLengthTrees( *p, "nexus" );
        if ( trees == NULL || trees->size() == 0 )
        {
            delete trees;
            trees = NULL;
            try
            {
                trees = readBranchLengthTrees( *p, "phylip" );
            }
            catch (RbException& e)
            {
                ;
            }
            catch (exception& e)
            {
                ;
            }
        }
        
        if ( trees == NULL || trees->size() == 0 )
        {
            delete trees;
            trees = readBranchLengthTrees( *p, "newick" );
        }
    }
    
// We cannot reset the tip node indices in case the tree topology changed during ancestral state monitoring.
// Instead see Tree::getTaxonBitSetMap()
//    if ( trees != NULL )
//    {
//        checkTreeTaxonIndices( trees );
//    }
    
    // print summary of results of file reading to the user
    if (readingDirectory == true)
    {
        std::stringstream o2;
        std::size_t size = 0;
        if ( trees == NULL || trees->size() == 0 )
        {
            o2 << "Failed to read any tree";
        }
        else if ( trees->size() == 1 )
        {
            size = trees->size();
            o2 << "Successfully read one tree";
        }
        else
        {
            size = trees->size();
            o2 << "Successfully read " << trees->size() << " trees";
        }
        RBOUT( o2.str() );
        std::set<std::string> myWarnings = getWarnings();
        if ( vectorOffile_names.size() - size > 0 && myWarnings.size() > 0 )
        {
            std::stringstream o3;
            if (vectorOffile_names.size() - size == 1)
            {
                o3 << "Did not read a file for the following ";
            }
            else
            {
                o3 << "Did not read " << vectorOffile_names.size() - size << " files for the following ";
            }
            if (myWarnings.size() == 1)
            {
                o3 << "reason:";
            }
            else
            {
                o3 << "reasons:";
            }
            RBOUT( o3.str() );
            for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
            {
                RBOUT( "* "+(*it) );
            }
        }
    }
    else
    {
        if ( trees != NULL && trees->size() > 0)
        {
            RBOUT( "Successfully read file" );
        }
        else
        {
            std::set<std::string> myWarnings = getWarnings();
            if ( myWarnings.size() > 0 )
            {
                std::stringstream o3;
                o3 << "Did not read the file for the following ";
                if (myWarnings.size() == 1)
                {
                    o3 << "reason:";
                }
                else
                {
                    o3 << "reasons:";
                }
                RBOUT( o3.str() );
                for (std::set<std::string>::iterator it = myWarnings.begin(); it != myWarnings.end(); it++)
                {
                    RBOUT( "* " + (*it) );
                }
            }
        }
    }

    
    return trees;
}

/** Read trees */
std::vector<Tree*>* NclReader::readBranchLengthTrees(const path &fn, std::vector<std::map<int,std::string> >& translationTables)
{

    std::vector<Tree*>* trees = readBranchLengthTrees(fn);
    getTranslateTables(translationTables);
    return trees;
}

/** Read trees */
std::vector<Tree*>* NclReader::readBranchLengthTrees(const path &file_name, const std::string &file_format)
{
    std::string fns = file_name.string();

    // check that the file exists
    if ( not exists(file_name) )
    {
        addWarning("Data file not found");
        return NULL;
    }
        
    try
    {

        if (file_format == "nexus")
        {
            // NEXUS file format
            nexusReader.ReadFilepath( fns.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        else if (file_format == "phylip")
        {
            // phylip file format with long taxon names
            nexusReader.ReadFilepath( fns.c_str(), MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT);
        }
        else if (file_format == "newick")
        {
            NewickTreeReader ntr;
            std::vector<Tree*>* trees = ntr.readBranchLengthTrees(fns);
            return trees;
        }
    }
    catch(NxsException& err)
    {
        std::stringstream o;
        o << "Nexus Error: " << err.msg << " (" << err.pos << ", " << err.line << ", " << err.col << ")";
        throw RbException( o.str() );
        return NULL;
    }
        
    std::vector<Tree*>* cvm = convertTreesFromNcl();
    
    return cvm;
}


std::vector<Tree*>* NclReader::readTimeTrees( const path &treefile_name )
{
    
    std::vector<Tree*>* trees = new std::vector<Tree*>();
    std::vector<Tree*>* m = readBranchLengthTrees( treefile_name );
    
    if (m != NULL)
    {
        for (std::vector<Tree*>::iterator it = m->begin(); it != m->end(); it++)
        {
            Tree* convertedTree = TreeUtilities::convertTree( *(*it) );
            delete (*it);
            trees->push_back( convertedTree );
        }
    }
    
    return trees;
}


/** Set excluded characters and taxa */
void NclReader::setExcluded( const NxsCharactersBlock* charblock, HomologousCharacterData* cMat ) const
{
    
    // Set excluded taxa
    for (unsigned int origTaxIndex=0; origTaxIndex<charblock->GetNTax(); origTaxIndex++ )
    {
        if ( !charblock->IsActiveTaxon( origTaxIndex ) )
        {
            cMat->excludeTaxon( origTaxIndex );
        }
    }
    
    // Set excluded characters
    NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
    for ( NxsUnsignedSet::const_iterator cit = excluded.begin(); cit != excluded.end(); cit++ )
        cMat->excludeCharacter( *cit );
}


/** Translate a single NCL tree into a RevBayes tree */
Tree* NclReader::translateNclSimpleTreeToBranchLengthTree(NxsSimpleTree& nTree, const NxsTaxaBlock *tb, bool rooted) {
    
    // get the root from the ncl tree
    const NxsSimpleNode* rn = nTree.GetRootConst();
    
    // create a new tree root node
    std::string name = rn->GetName();
    if ( rn->GetTaxonIndex() < tb->GetNumTaxonLabels() )
    {
        name = tb->GetTaxonLabel( rn->GetTaxonIndex() ).BlanksToUnderscores();
    }
    TopologyNode* root = new TopologyNode(name);
    
    // create a map which holds for each node a map of name value pairs.
    //    std::map<const TopologyNode*, std::map<std::string, RevObject*> > nodeParameters;
    
    // the vector of nodes so that we remember the ordering (and matching) which node belongs to which parameter
    std::vector<TopologyNode*> nodes;
    
    
    // the vector of branch lengths that we use later to set the branch length in the tree object
    std::vector<double> brlens;
    
    
    // the vector of node indices that we use later to set the indices in the tree object
    bool reindex = false;
    
    // construct tree recursively
    constructBranchLengthTreefromNclRecursively(root, nodes, brlens, reindex, rn, tb);

    // create a new simple tree
    Tree* tau = new Tree();
    
    // initialize the topology by setting the root
    tau->setRoot(root, reindex);

    // finally set the branch lengths
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
        tau->getNode(nodes[i]->getIndex()).setBranchLength( brlens[i] );
    }
    
    // only trees with 2-degree root nodes are rooted trees.
    tau->setRooted( root->getNumberOfChildren() == 2 );

    return tau;    
}

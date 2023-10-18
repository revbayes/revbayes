/**
 * @file
 * This file contains the declaration of NclReader, which is the wrapper class for the NCL
 * for reading character matrices, sequences, and trees.
 
 * @brief Declaration of NclReader
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-16, version 1.0
 * @extends DAGNode
 *
 * $Id$
 */

#ifndef NclReader_H
#define NclReader_H

#include "ncl.h"
#include "nxsmultiformat.h"
#include "RbFileManager.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    class AminoAcidState;
    class AbstractCharacterData;
    class BranchLengthTree;
    class ContinuousCharacterData;
    class DnaState;
    class HomologousCharacterData;
    class RnaState;
    class StandardState;
    class Tree;
    class TimeTree;
    class TopologyNode;

    template<class charType>
    class HomologousDiscreteCharacterData;
    
    template<class charType>
    class NonHomologousDiscreteCharacterData;

    
    class NclReader{
        
        friend class ::NxsBlock;
        friend class ::NxsUnalignedBlock;
        
    public:
                                                                NclReader(void);                                                                //!< Default constructor
                                                                NclReader(const NclReader& r);                                                  //!< Copy constructor
        virtual                                                ~NclReader(void) { }                                                             //!< Destructor

        void                                                    addWarning(std::string s) { warningsSummary.insert(s); }                        //!< Add a warning to the warnings vector
        void                                                    clearWarnings(void) { warningsSummary.clear(); }                                //!< Clear all of the warnings from the warnings vector
        size_t                                                  getNumWarnings(void) { return warningsSummary.size(); }                         //!< Return the number of warnings
        std::set<std::string>&                                  getWarnings(void) { return warningsSummary; }                                   //!< Get a reference to the warnings vector
        void                                                    getTranslateTables(std::vector<std::map<int,std::string> >& translateTables);
        void                                                    clearContent(void) { nexusReader.ClearContent(); }                              //!< Clear the content of the NCL object
        
        // file type methods
        bool                                                    isFastaFile(const path& fn, std::string& dType);                                //!< Checks if the file is in Fasta format
        bool                                                    isNexusFile(const path& fn);                                                    //!< Checks if the file is in NEXUS format
        bool                                                    isPhylipFile(const path& fn, std::string& dType, bool& isInterleaved);          //!< Checks if the file is in Phylip format
        
        // alignment functions
        std::vector<AbstractCharacterData* >                    readMatrices(const path &fn);
        std::vector<AbstractCharacterData* >                    readMatrices(const path &fn, const std::string &ft);                            //!< Read character data from a single file (with file format info too)
        std::vector<AbstractCharacterData* >                    readMatrices(const std::map<path,std::string>& fileMap);                 //!< Read a list of file names contained in a map (with file format info too)
        std::vector<AbstractCharacterData* >                    readMatrices(const std::vector<path> fn, const std::string fileFormat, const std::string dataType, const bool isInterleaved);  //!< Read a list of file names contained in a vector of strings
        
        // functions for reading trees
        std::vector<Tree* >*                                    readBranchLengthTrees(const path &fn);                                   //!< Read trees
        std::vector<Tree* >*                                    readBranchLengthTrees(const path &fn, std::vector<std::map<int,std::string> >& translationTables);
        std::vector<Tree* >*                                    readTimeTrees(const path &treeFilename);
        
    private:
        
        void                                                    checkTreeTaxonIndices( std::vector<Tree* >* trees );                            //!< Check and change the indices of the trees to match
        HomologousDiscreteCharacterData<AminoAcidState>*        createAminoAcidMatrix(NxsCharactersBlock* charblock);                           //!< Create an object to hold amino acid data
        ContinuousCharacterData*                                createContinuousMatrix(NxsCharactersBlock* charblock);                          //!< Create an object to hold continuous data
        HomologousDiscreteCharacterData<DnaState>*              createDnaMatrix(NxsCharactersBlock* charblock);                                 //!< Create an object to hold DNA data
        HomologousDiscreteCharacterData<RnaState>*              createRnaMatrix(NxsCharactersBlock* charblock);                                 //!< Create an object to hold RNA data
        HomologousDiscreteCharacterData<StandardState>*         createStandardMatrix(NxsCharactersBlock* charblock);                            //!< Create an object to hold standard data
        NonHomologousDiscreteCharacterData<AminoAcidState>*     createUnalignedAminoAcidMatrix(NxsUnalignedBlock* charblock);                   //!< Create an object to hold amino acid data
        NonHomologousDiscreteCharacterData<DnaState>*           createUnalignedDnaMatrix(NxsUnalignedBlock* charblock);                         //!< Create an object to hold DNA data
        NonHomologousDiscreteCharacterData<RnaState>*           createUnalignedRnaMatrix(NxsUnalignedBlock* charblock);                         //!< Create an object to hold RNA data
        std::string                                             intuitDataType(std::string& s);                                                 //!< Attempt to determine the type of data
        
        // methods for reading sequence alignments
        std::vector<AbstractCharacterData* >                    convertFromNcl(const path& fileName);                                           //!< Reads the blocks stored by NCL and converts them to RevBayes character matrices
        std::vector<AbstractCharacterData* >                    readMatrices(const path &fn, const std::string &format, const std::string &type, bool isInterleaved);           //!< Reads a single file using NCL
        void                                                    setExcluded(const NxsCharactersBlock* charblock, HomologousCharacterData* cMat ) const;                         //!< Set excluded taxa and excluded characters
        
        // methods for reading trees
        void                                                    constructBranchLengthTreefromNclRecursively(TopologyNode* tn, std::vector<TopologyNode*> &nodes, std::vector<double> &brlens, bool& reindex, const NxsSimpleNode* tnNcl, const NxsTaxaBlock *tb);  //!< Constructs a tree from NCL
     // void                                                    constructTreefromNclRecursively(TopologyNode* tn, const NxsSimpleNode* tnNcl, const NxsTaxaBlock *tb);  //!< Constructs a tree from NCL
        std::vector<Tree* >*                                    readBranchLengthTrees(const path &fn, const std::string &fileFormat);            //!< Read trees
        std::vector<Tree* >*                                    convertTreesFromNcl(void);                                                       //!< Converts trees stored by NCL into RevBayes formatted trees
        Tree*                                                   translateNclSimpleTreeToBranchLengthTree(NxsSimpleTree &nTree,const NxsTaxaBlock *tb, bool rooted);  //!< Translate a single NCL tree into a RevBayes tree
        
        MultiFormatReader                                       nexusReader;                                                                     //!< The NCL object that reads the files
        std::set<std::string>                                   warningsSummary;                                                                 //!< A vector that contains the warnings that acumulate
    };
}

#endif

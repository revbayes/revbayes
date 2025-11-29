#include <fstream>
#include <ostream>
#include <sstream>
#include <random>
#include <ios>
#include <string>
#include <vector>
#include <algorithm>

#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"
#include "FastaFileToNaturalNumbersConverter.h"

// FIXME: This file shares a lot of code with CountFileToNaturalNumbersConverter.cpp

using namespace RevBayesCore;



/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
FastaFileToNaturalNumbersConverter::FastaFileToNaturalNumbersConverter( void ) 
{
    
}


/** Read Count File and Write Natural Numbers file */
void FastaFileToNaturalNumbersConverter::faconverter( const path &fi, const std::vector<std::string> &taxa, const std::vector<std::string> &alleles , const size_t& n_individuals, const path &fo )
{
  
  // getting the number of alleles and taxa
  size_t n_alleles = alleles.size(); 
  size_t n_taxa    = taxa.size();

  // open file
  std::ifstream readStream( fi.string() );
  if ( not readStream )
      throw RbException()<<"Could not open file "<<fi<<".";

  // some important quantities to parse the fasta file    
  std::vector<std::string> alignment;
  std::vector<int>      ns_taxa(n_taxa,0);
  std::vector<int>      aligment_index;
  std::string              line,seq_name;
  size_t                   index;

  std::vector<size_t> length_taxa (n_taxa);
  for (size_t i =0; i<n_taxa; ++i){
    length_taxa[i] = taxa[i].length();
  }
  //for (size_t i =0; i<n_taxa; ++i){ std::cout << "lt: " <<length_taxa[i] << "\n"; }

  // going through the alingment lines
  while (safeGetline(readStream,line)) {

    if (line.length() == 0) { break; }

    for (size_t t =0; t<n_taxa; ++t){

      seq_name = line.substr(1,length_taxa[t]);
      index    = getIndex(seq_name,taxa);
      //std::cout << "tname:" << seq_name << " index: " << index << "\n";

      // if the element is found keep the sequence (sitting in the next line) and save the taxa index
      if ( index < n_taxa )  {
        safeGetline(readStream,line);
        alignment.push_back(line);
        aligment_index.push_back(index);
        ns_taxa[index] += 1;
        break;

      // if the sequence indentifier does not include a taxa names throw an error 
      } else {
        throw RbException() << "Sequence \"" << line << "\" does not belong to any of the taxa. Make sure every sequence starts with one of the taxa names."; 
      }

    }
  }

  // getting the number of sampled individuals and sites
  size_t n_samples = alignment.size();
  size_t n_sites = alignment[0].length();

  // defining the edge vector and matrix 
  //number of edges
  size_t n_edges = (n_alleles*n_alleles - n_alleles)/2;
  
  // vector and matrix of edges 
  // vector of edges indexed as [edge_index] = "a_ia_j"
  // matrix of edges indexed as [allele_index,edge_index]
  // it indicates which edges have a certain allele: matrix of 0s and 1s
  std::vector<int> matrix_edges(n_alleles*n_edges,0);
  std::vector<std::string> vector_edges(n_edges);
  
  // generating all the possible pairwise combinations of alleles
  size_t edge = 0;
  for (size_t i=0; i<n_alleles; ++i ){
    for (size_t j=(i+1); j<n_alleles; ++j){
    
      vector_edges[edge] = std::to_string(i)+std::to_string(j);
      matrix_edges[i*n_edges+edge] = 1;
      matrix_edges[j*n_edges+edge] = 1;
      edge += 1;
      
    }
  }
 
  // cloning the taxa vector to build the final Natuaral Numbers alignment
  std::vector<std::string> ctaxa = taxa ;

  // parsing the allelic counts
  std::string allele;
  std::vector<int> counts(n_alleles,0);
  size_t state;

  for (size_t t=0; t < n_taxa; ++t){
    
    for (size_t s=0; s<n_sites; ++s){

      std::fill(counts.begin(), counts.end(), 0);

      for (size_t n=0; n<n_samples; ++n) {

        if (aligment_index[n]==t){
          allele = alignment[t].substr(s,1);
          index  = getIndex(allele,alleles);
          //std::cout << "pop: " << t << " allele: " << allele << " index: " << index << "\n";

          // counting only if the alleles that present in alleles
          // not throwing an error here because users usually have plenty of Ns and -s in their alignments
          if (index < n_alleles){ counts[index] += 1; }

        }
      }

      state = getState(counts,n_alleles,n_individuals,vector_edges,matrix_edges,s,taxa[t]);
      ctaxa[t] += " " + std::to_string(state);

    }
  }

  // close the input file connection
  readStream.close();
  //for (size_t i =0; i<n_taxa; ++i){ std::cout << ctaxa[i] << "\n"; }

  // summarizing the 
  std::cout <<     "\n  Number of taxa                  " << n_taxa <<
                   "\n  Number of sampled sequences     " << n_samples << 
                   "\n  Number of alleles               " << n_alleles << 
                   "\n  Number of sites                 " << n_sites <<
                   "\n  Number of individuals           " << n_individuals <<
                   "\n  Number of PoMo states           " << n_alleles*(1.0+(n_alleles-1.0)*(n_individuals-1.0)*0.5) << "\n\n";

  // the filestream object
  createDirectoryForFile( fo );

  // open the stream to the file a write it
  std::ofstream NaturalNumbersStream( fo.string() );

  for (size_t i=0; i<n_taxa; ++i)
  {
    NaturalNumbersStream << ctaxa[i] + "\n";
  }
  
  // close the stream
  NaturalNumbersStream.close();

}


const size_t FastaFileToNaturalNumbersConverter::getIndex(std::string& taxa_name, const std::vector<std::string>& taxa) 
{
   
  for (size_t i =0; i<taxa.size(); ++i){
    if (taxa[i] == taxa_name) { return i; }
  }

  return taxa.size();
}




const size_t FastaFileToNaturalNumbersConverter::getState(std::vector<int>& counts, size_t& n_alleles, const size_t & n_individuals, std::vector<std::string>& vector_edges, std::vector<int>& matrix_edges, size_t& s, const std::string& taxa_name )
{
  
        
        // some variables
        std::string str_index;
        size_t state=-1, int_index=-1, m=-1;


        // setting total counts, the last postive count (why last? important for state indexing), and number of non-null counts to 0
        size_t M        = 0;
        size_t n_counts = 0;
      
        // goes through the number of alleles
        // counts are comma separated
        for (size_t k=0; k<n_alleles; k++){
        
          size_t value = counts[k];
        
          // getting info from the count pattenr
          if (value > 0) {
            M        += value;
            m         = value;
            n_counts += 1;
            int_index = k;
            str_index = str_index + std::to_string(k);
          }
        
        }
      
        // pointing out some typical invalid counts: null counts (e.g., 0,0,0,0) and >2-allelic counts (e.g., 0,1,1,1)
        if (n_counts==0){
          throw RbException() << "Unexpected allelic counts at site " << s << " and taxa " << taxa_name << ". PoMos require at least one postive count."; 
        }
        if (n_counts>2){
          throw RbException() << "Unexpected allelic couts at site " << s << " and taxa " << taxa_name << ". PoMos only accept monoallelic or biallelic counts."; 
        }
      
        // sampling a 0:n_individuals frequency from the weight vector 
        size_t weight = sample_weight(M, m, n_individuals);
      
        // determining the pomo state
        // three possible situations
        // if the count is monoallelic & likely "sampled" from a fixed state
        if (weight==n_individuals){
        
          state = int_index;
          //std::cout << "  " << state << "\n";

        // if the count is monoallelic & likely "sampled" from a polymoprhic state
        } else if (n_counts==1 && weight<n_individuals ) {
        
          size_t edge = sample_edge(int_index, matrix_edges);
          state = n_alleles+edge*n_individuals-edge+weight-1;
          //std::cout << "  " << state << "\n";
      
        // if the count is biallelic, thus necessarily sampled from a polymoprhic state
        } else if (n_counts>1) {
        
          size_t edge = get_index(vector_edges,str_index);
          state = n_alleles+edge*n_individuals-edge+weight-1;
          //std::cout << "  " << state << "\n";
        
        }

  return state;
  
}


// samples a state from a pomo edge 0:(N-1) givent the weights of a 
// binomal distribution B(m|n/N,M)
// based on equation (13) of Schrempf et al. (2016) JTB
const size_t FastaFileToNaturalNumbersConverter::sample_weight(size_t& M, size_t& m, const size_t& N)
{
  
  std::vector<double> weights(N+1);
  double prob;
  
  // calculating the weight vector
  for (size_t i=0; i < N+1; ++i){
    prob = 1.0*i/N;
    weights[i] = pow(prob,m)*pow(1-prob,M-m);
  }
  
  // sampling a pomo state from the weight vector 
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  
  return d(gen);
  
}


// samples a pomo edge from a matrix of allele_indexes*edge_indexes matrix
// the matrix estipulates wich edges contain the allele allele_index
const size_t FastaFileToNaturalNumbersConverter::sample_edge(size_t& allele_index, std::vector<int>& vector){
  
  // setting the appropriate sub vector to sample from
  // the edge matrix is indexed as [allele_index,edge_index]
  // the sub vector pics the allele_index line
  std::vector<int> sub_vector = {vector.begin() + allele_index*6, vector.begin() + allele_index*6 + 6};

  // sampling a pomo edge
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d( sub_vector.begin(), sub_vector.end() );

  return d(gen);
  
}


// gets the edge of an observed polymorphic count
const size_t FastaFileToNaturalNumbersConverter::get_index(std::vector<std::string>& vector, std::string& element) {
  
  for (size_t i=0; i<vector.size(); ++i){
    if (vector[i]==element){
      return i;
    }
  }
  return -1;
}



//const size_t FastaFileToNaturalNumbersConverter::getNumberOfAlleles( void )
//{
//    return numberOfAlleles_;
//}

//const size_t FastaFileToNaturalNumbersConverter::getNumberOfIndividuals( void )
//{
//    return numberOfIndividuals_;
//}

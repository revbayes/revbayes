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
#include "CountFileToNaturalNumbersConverter.h"

// FIXME: This file shares a lot of code with FastaFileToNaturalNumbersConverter.cpp

using namespace RevBayesCore;



/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
CountFileToNaturalNumbersConverter::CountFileToNaturalNumbersConverter( void ) 
{
    
}


/** Read Count File and Write Natural Numbers file */
void CountFileToNaturalNumbersConverter::cfconverter( const path &fi, const size_t n_individuals, const path &fo )
{
  
    // open file
    std::ifstream readStream( fi.string() );
    if ( not readStream )
    {
        throw RbException()<<"Could not open file "<<fi<<".";
    }
    

    // reading the count file 
    // first ignoring all the commented lines
    std::string line = "";
    while (safeGetline(readStream,line))
    {
        if ( line.at(0) == '#' ) { continue; } 
        else { break; }
    }

    // reading the actual first line from the count file
    // saving taxa names
    size_t n_taxa;
    std::vector<std::string> taxa;

    std::string cell = "";
    std::stringstream ss0(line);

    // ignoring the first two elements CHRS and POS
    std::getline(ss0,cell,' ');
    std::getline(ss0,cell,' ');
    while ( std::getline(ss0,cell,' ') ) {
      taxa.push_back(cell);
    }

    n_taxa = taxa.size();
    
    //for (size_t i=0; i<n_taxa; ++i){ std::cout << taxa[i] << "\n"; }
    //std::cout << "n_taxa:" << n_taxa << "\n";

 
    // getting the first count patter to determine the number of alleles
    safeGetline(readStream,line);
    std::stringstream ss1(line);

    // again ignoring the first two elements
    std::getline(ss1,cell,' ');
    std::getline(ss1,cell,' ');

    // getting the first count and setting the number of alleles
    std::getline(ss1,cell,' ');
    size_t n_alleles = std::count(cell.begin(), cell.end(), ',') + 1;
    //std::cout << "N_alleles:" << n_alleles << "\n"; 
    
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

    // we are now in conditions for sampling the state for the first count 
    size_t state;
    state = getState(cell,n_alleles,n_individuals,vector_edges,matrix_edges);
    taxa[0] += " " + std::to_string(state);

    // and now the remaining counts on this first line
    size_t index = 1;
    while ( std::getline(ss1,cell,' ') ) {
      state = getState(cell,n_alleles,n_individuals,vector_edges,matrix_edges);
      taxa[index] += " " + std::to_string(state);
      index += 1;
    }
    //for (size_t i=0; i<n_taxa; ++i){ std::cout << taxa[i] << "\n"; }
       
    // and now the lines
    size_t n_sites = 1; 

    while (safeGetline(readStream,line)){

        std::stringstream ss2(line);

        // again ignoring the first two elements
        std::getline(ss2,cell,' ');
        std::getline(ss2,cell,' ');

        // and now the remaining counts on this first line
        index = 0;
        while ( std::getline(ss2,cell,' ') ) {
          state = getState(cell,n_alleles,n_individuals,vector_edges,matrix_edges);
          taxa[index] += " " + std::to_string(state);
          index += 1;
        }
        
        n_sites += 1;

    }
    //for (size_t i=0; i<n_taxa; ++i){ std::cout << taxa[i] << "\n"; }
    
    // close the input file connection
    readStream.close();

    // summarizing the 
    std::cout <<     "\n  Number of taxa                  " << n_taxa <<
                     "\n  Number of alleles               " << n_alleles << 
                     "\n  Number of sites                 " << n_sites -1 <<
                     "\n  Number of individuals           " << n_individuals <<
                     "\n  Number of PoMo states           " << n_alleles*(1.0+(n_alleles-1.0)*(n_individuals-1.0)*0.5) << "\n\n";

    std::string alignment = "";
    for (size_t i=0; i<n_taxa; ++i){
      alignment += taxa[i] + "\n";
    }

    // the filestream object
    std::fstream NaturalNumbers;
    
    createDirectoryForFile(fo);
    
    // open the stream to the file
    NaturalNumbers.open( fo.string(), std::fstream::out );
    
    NaturalNumbers << alignment;
  
    // close the stream
    NaturalNumbers.close();



}


const size_t CountFileToNaturalNumbersConverter::getState(std::string counts, size_t n_alleles, size_t n_individuals, std::vector<std::string>& vector_edges, std::vector<int>& matrix_edges )
{
  
        std::stringstream sss( counts );
        std::string count,str_index;
        
        // some variables
        size_t state = -1, int_index = -1, m = -1;

        // setting total counts, the last postive count (why last? important for state indexing), and number of non-null counts to 0
        size_t M        = 0;
        size_t n_counts = 0;
      
        // goes through the number of alleles
        // counts are comma separated
        for (size_t k=0; k<n_alleles; k++){
        
          getline( sss, count, ',' );
          size_t value = stoi(count);
        
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
          throw RbException() << "Unexpected count pattern: " << counts << ". PoMos require at least one postive count."; 
        }
        if (n_counts>2){
          throw RbException() << "Unexpected count pattern: " << counts << ". PoMos only accept monoallelic or biallelic counts."; 
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
const size_t CountFileToNaturalNumbersConverter::sample_weight(size_t M, size_t m, size_t N)
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
const size_t CountFileToNaturalNumbersConverter::sample_edge(size_t allele_index, std::vector<int>& vector){
  
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
const size_t CountFileToNaturalNumbersConverter::get_index(std::vector<std::string>& vector, std::string element) {
  
  for (size_t i=0; i<vector.size(); ++i){
    if (vector[i]==element){
      return i;
    }
  }
  return -1;
}



//const size_t CountFileToNaturalNumbersConverter::getNumberOfAlleles( void )
//{
//    return numberOfAlleles_;
//}

//const size_t CountFileToNaturalNumbersConverter::getNumberOfIndividuals( void )
//{
//    return numberOfIndividuals_;
//}

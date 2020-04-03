#ifndef InferAncestralPopSizeFunction_h
#define InferAncestralPopSizeFunction_h

#include "RlTypedFunction.h"
#include "MatrixReal.h"
// #include "TypedFunction.h"
// anything other includes?

#include <string>
#include <vector>

namespace RevBayesCore {
class DagNode;
class Tree;
class MatrixReal;
class Clade;
class Taxon;

template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
  
  /**
     * @brief Declaration of the deterministic variable for ancestral population size estimation.
     * @return The probability density of the number of hidden lineages through time, computed according to Manceau & al. 2020 algorithm (http://dx.doi.org/10.1101/755561)
     */

  class InferAncestralPopSizeFunction : public TypedFunction<MatrixReal> {

  public:
      InferAncestralPopSizeFunction(                        const TypedDagNode<double> *sa,
                                                            const TypedDagNode<double> *l,
                                                            const TypedDagNode<double> *m,
                                                            const TypedDagNode<double> *p,
                                                            const TypedDagNode<double> *o,
                                                            const TypedDagNode<double> *rho,
                                                            const TypedDagNode<double> *r,
                                                            const TypedDagNode<long> *n,

                                                            const std::string& cdt,
                                                            const TypedDagNode< RevBayesCore::RbVector<double> > *O,  
                                                            const std::vector<double> &tau,
                                                            bool uo,
                                                            TypedDagNode<Tree> *tr);
      
      virtual                                               ~InferAncestralPopSizeFunction(void);

      // public member functions
      InferAncestralPopSizeFunction*                        clone(void) const;                     //!< Create an independent clone
      void                                                  update(void);                          //!< Update the value of the function


  protected:
      // Parameter management functions
      void                                                  swapParameterInternal(const DagNode *oldP, const DagNode *newP);     //!< Swap a parameter

  private:
      // Members
      const TypedDagNode< double > *                        start_age;                             //!< Start age of the process.
      const TypedDagNode< double > *                        lambda;                                //!< The speciation rate.
      const TypedDagNode< double > *                        mu;                                    //!< The extinction rate.
      const TypedDagNode< double > *                        psi;                                   //!< The sampling probability of a just extinct species.
      const TypedDagNode< double > *                        omega;                                 //!< The occurrence sampling rate.
      const TypedDagNode< double > *                        rho;                                   //!< The sampling probability of extant taxa.
      const TypedDagNode< double > *                        removalPr;                             //!< The removal probability after sampling.
      const TypedDagNode< long > *                          maxHiddenLin;                          //!< The maximal number of hidden lineages.
      const std::string&                                    cond;                                  //!< Condition of the process ("time" or "survival")
      const std::vector<double>                             time_points;                           //!< Times at which density is computed
      const TypedDagNode< RbVector<double> > *              occurrences;                           //!< Occurrence ages of incomplete fossils
      const bool                                            useOrigin;                             //!< Start the process at the origin (otherwise root)
      const TypedDagNode< Tree > *                          timeTree;                              //!< Facultative initial tree

    };


}

#endif

#ifndef InferAncestralPopSizeFunctionPiecewise_h
#define InferAncestralPopSizeFunctionPiecewise_h

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

  class InferAncestralPopSizeFunctionPiecewise : public TypedFunction<MatrixReal> {

  public:
      InferAncestralPopSizeFunctionPiecewise(               const TypedDagNode<double> *sa,
                                                            const DagNode *inspeciation,
                                                            const DagNode *inextinction,
                                                            const DagNode *inserialsampling,
                                                            const DagNode *inoccurrence,
                                                            const DagNode *ineventsampling,
                                                            const DagNode *intreatment,
                                                            const TypedDagNode<long> *n,
                                                            const std::string& cdt,
                                                            const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                            const std::vector<double> &tau,
                                                            bool uo,
                                                            TypedDagNode<Tree> *tr,
                                                            const TypedDagNode< RbVector<double> > *ht);

      virtual                                               ~InferAncestralPopSizeFunctionPiecewise(void);

      // public member functions
      InferAncestralPopSizeFunctionPiecewise*                         clone(void) const;                     //!< Create an independent clone
      void                                                            update(void);                          //!< Update the value of the function
      void                                                            updateVectorParameters(void) const;

  protected:
      // Parameter management functions
      void                                                  swapParameterInternal(const DagNode *oldP, const DagNode *newP);     //!< Swap a parameter

  private:
      // Members
      const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous birth rates.
      const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous birth rates.
      const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous death rates.
      const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous death rates.
      const TypedDagNode<double >*                    homogeneous_r;                                         //!< The homogeneous conditional probability of death upon treatment.
      const TypedDagNode<RbVector<double> >*          heterogeneous_r;
      const TypedDagNode<double >*                    homogeneous_o;                                         //!< The homogeneous conditional probability of death upon treatment.
      const TypedDagNode<RbVector<double> >*          heterogeneous_o;                                         //!< The heterogeneous conditional probability of death upon treatment.
      const TypedDagNode<double >*                    homogeneous_phi;                                       //!< The homogeneous sampling rates.
      const TypedDagNode<RbVector<double> >*          heterogeneous_phi;                                     //!< The heterogeneous sampling rates.
      const TypedDagNode<double >*                    homogeneous_Phi;                                       //!< The probability of sampling a tip at the present.
      const TypedDagNode<RbVector<double> >*          heterogeneous_Phi;                                     //!< The probability of sampling individuals at set time intervals.
      const TypedDagNode<RbVector<double> >*          interval_times;




      const TypedDagNode< double > *                        start_age;                             //!< Start age of the process.
      const TypedDagNode< long > *                          maxHiddenLin;                          //!< The maximal number of hidden lineages.
      const std::string&                                    cond;                                  //!< Condition of the process ("time" or "survival")
      const std::vector<double>                             time_points;                           //!< Times at which density is computed
      const TypedDagNode< RbVector<double> > *              occurrences;                           //!< Occurrence ages of incomplete fossils
      const bool                                            useOrigin;                             //!< Start the process at the origin (otherwise root)
      const TypedDagNode< Tree > *                          timeTree;                              //!< Facultative initial tree

      mutable std::vector<double>                     lambda;
      mutable std::vector<double>                     mu;
      mutable std::vector<double>                     phi;
      mutable std::vector<double>                     r;
      mutable std::vector<double>                     omega;
      mutable std::vector<double>                     phi_event;
      mutable std::vector<double>                     timeline;                                              //!< The times of the instantaneous events and rate shifts.


    };


}

#endif

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
                                                            const TypedDagNode< RbVector<double> > *tau,
                                                            bool uo,
                                                            TypedDagNode<Tree> *tr);
      
      virtual                                               ~InferAncestralPopSizeFunction(void);

      // public member functions
      InferAncestralPopSizeFunction*                        clone(void) const;                                                          //!< Create an independent clone
      void                                                  update(void);                                                           //!< Update the value of the function
      void                                                  poolTimes(void) const;
      MatrixReal                                            ComputeLt(void) const;
      MatrixReal                                            ComputeMt(void) const;
      // void                                                  ComputeKt(void) const;                                //!< Compute the log-transformed probability of the current value.

  protected:
      // Parameter management functions
      void                                                  swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

      struct Event {
        Event(double d, std::string s) : time(d), type(s) {};

        double time;
        std::string type;

        double getEventTime(void){ return time; }
        std::string getEventType(void){ return type; }

      };

      // Vector of Events
      mutable std::vector<Event>         events;

      // Sorting functions
      struct AgeCompare {
        bool operator()(const Event first, const Event second) {
                  return first.time < second.time;
        }
      };

      struct AgeCompareReverse {
        bool operator()(const Event first, const Event second) {
                  return first.time > second.time;
        }
      };

  private:
      // Members
      const TypedDagNode< double > *                        tor;                                                                            //!< Time of origin.
      const TypedDagNode< double > *                        lambda;                                                                         //!< The speciation rate.
      const TypedDagNode< double > *                        mu;                                                                             //!< The extinction rate.
      const TypedDagNode< double > *                        psi;                                                                            //!< The sampling probability of a just extinct species.
      const TypedDagNode< double > *                        omega;                                                                          //!< The occurrence sampling rate.
      const TypedDagNode< double > *                        rho;                                                                            //!< The sampling probability of extant taxa.
      const TypedDagNode< double > *                        removalPr;                                                                      //!< The removal probability after sampling.
      const TypedDagNode< long > *                          maxHiddenLin;                                                                   //!< The maximal number of hidden lineages.
      const TypedDagNode< RbVector< double > > *            dn_time_points;                                                                 //!< Times at which density is computed
      const TypedDagNode< Tree > *                          timeTree;                                                                       //!< Facultative initial tree
      mutable size_t                                        extant;                                                                         //!< Number of extant taxa
    };


}

#endif

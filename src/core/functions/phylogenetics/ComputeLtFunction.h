#ifndef ComputeLtFunction_h
#define ComputeLtFunction_h

#include "TypedFunction.h"
// anything other includes?

#include <string>
#include <vector>

namespace RevBayesCore {
// anything else here?
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;

  // class Event {
  //
  // //public: Event(double time, std::<string> type, bool extant) : time(time), type(type), extant(extant) {}
  //
  // public: Event(double time, std::string type, bool extant) : time(time), type(type), extant(extant) {}
  //
  // double getEventTime(void){ return time; }
  // std::string getEventType(void){ return type; }
  // bool getIsExtant(void){ return extant; }
  // void setIsExtant(bool b){ extant = b; }
  //
  // double time;
  //
  // private:
  //   std::string type;
  //   bool extant;
  //
  // };

  class ComputeLtFunction : public TypedFunction<double> {

  public:
    ComputeLtFunction(
      const TypedDagNode< RbVector<double> > *a,
      const TypedDagNode< RbVector<double> > *b,
      const TypedDagNode< RbVector<double> > *c,
      const TypedDagNode< RbVector<double> > *d,
      const TypedDagNode< RbVector<double> > *e,
      const TypedDagNode< RbVector<double> > *f,
      const TypedDagNode< double > *t,
      const TypedDagNode< double > *l,
      const TypedDagNode< double > *m,
      const TypedDagNode< double > *p,
      const TypedDagNode< double > *o,
      const TypedDagNode< double > *rho,
      const TypedDagNode< double > *r,
      const TypedDagNode< RbVector<double> > *g
    );
    virtual                                            ~ComputeLtFunction(void);

    // public member functions
    ComputeLtFunction*                                  clone(void) const;                                                          //!< Create an independent clone
    void                                                update(void);
    void                                                poolTimes(void);

  protected:
    void                     swapParameterInternal(const DagNode *oldP, const DagNode *newP);


    struct Event {
      Event(double d, std::string s, bool b) : time(d), type(s), extant(b) {};

      double time;
      std::string type;
      bool extant;

      double getEventTime(void){ return time; }
      std::string getEventType(void){ return type; }
      bool getIsExtant(void){ return extant; }

    };

    // vector of Events
    std::vector<Event>         events;

    struct AgeCompare {
      bool operator()(const Event first, const Event second) {
                return first.time < second.time;
      }
    };

  private:

    // members
    const TypedDagNode< RbVector< double > >*           listA;
    const TypedDagNode< RbVector< double > >*           listB;
    const TypedDagNode< RbVector< double > >*           listC;
    const TypedDagNode< RbVector< double > >*           listD;
    const TypedDagNode< RbVector< double > >*           listE;
    const TypedDagNode< RbVector< double > >*           listF;
    const TypedDagNode< RbVector< double > >*           listG;

    const TypedDagNode< double > *                      tor;
    const TypedDagNode< double > *                      lambda;
    const TypedDagNode< double > *                      mu;
    const TypedDagNode< double > *                      psi;
    const TypedDagNode< double > *                      omega;
    const TypedDagNode< double > *                      rho;
    const TypedDagNode< double > *                      removalPr;

  };


}

#endif

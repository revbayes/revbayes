// #ifndef TimeCalibratedSBN_H
// #define TimeCalibratedSBN_H

// #include "MatrixReal.h"
// #include "Subsplit.h"
// #include "SBNParameters.h"
// #include "Taxon.h"
// #include "Tree.h"
// #include "TypedDagNode.h"
// #include "TypedDistribution.h"

// namespace RevBayesCore {

//     class TimeCalibratedSBN : public TypedDistribution<Tree> {

//     public:
//         TimeCalibratedSBN(const SBNParameters parameters);
// 		virtual                                            ~TimeCalibratedSBN(void);                                                                    //!< Virtual destructor

//         // public member functions
//         TimeCalibratedSBN*                        clone(void) const;                                                                                  //!< Create an independent clone
//         double                                              computeLnProbability(void);
//         double                                              computeLnProbabilityNodeTimes( void );
//         void                                                redrawValue(void);
//         virtual void                                        setValue(Tree *v, bool f=false);                                    //!< Set the current value, e.g. attach an observation (clamp)

//     protected:
//         // Parameter management functions
//         void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

//     private:

//         // helper functions
//         void                                                simulateTree(void);

//         // members
//         SBNParameters                                       parameters;  // The parameters for this SBN
//         std::vector<Taxon>                                  taxa;
//     };

// }

// #endif

#ifndef ComputeLikelihoodsLtMtFunction_H
#define ComputeLikelihoodsLtMtFunction_H

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "MatrixReal.h"

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
     * @brief Declaration of the deterministic variable for Lt/Mt likelihoods.
     */
    
    class ComputeLikelihoodsLtMtFunction : public TypedFunction<MatrixReal> {
        
    public:
        ComputeLikelihoodsLtMtFunction(                     const TypedDagNode<double> *sa,
                                                            const TypedDagNode<double> *l,
                                                            const TypedDagNode<double> *m,
                                                            const TypedDagNode<double> *p,
                                                            const TypedDagNode<double> *o,
                                                            const TypedDagNode<double> *rho,
                                                            const TypedDagNode<double> *r,
                                                            const TypedDagNode<long> *n,

                                                            const std::string& cdt,
                                                            const std::vector<double> &tau,
                                                            bool uo,
                                                            bool mt,
                                                            const TypedDagNode<Tree> *tr);                                      //!< Default constructor
        
        virtual                                             ~ComputeLikelihoodsLtMtFunction(void);                              //!< Destructor
                
        // Basic utility functions
        ComputeLikelihoodsLtMtFunction*                     clone(void) const;                                                  //!< Clone object
        void                                                update(void);                                                       //!< Clone the function

        // Likelihoods computation functions
        void                                                poolTimes(void) const;
        MatrixReal                                          ComputeLt(void) const;
        MatrixReal                                          ComputeMt(void) const;
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
        // Structure of the phylodynamic events in a tree with occurrences : type + time
        struct Event {
            Event(double d, std::string s) : time(d), type(s) {};

            std::string type;
            double time;

            std::string getEventType(void){ return type; }
            double getEventTime(void){ return time; }

        };

        // Vector of events
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
        RevBayesCore::DistanceMatrix* getDistanceMatrix(const TypedDagNode<Tree>& tree);
        

        // members
        const TypedDagNode< double > *                        start_age;                             //!< Time of origin.
        const TypedDagNode< double > *                        lambda;                                //!< The speciation rate.
        const TypedDagNode< double > *                        mu;                                    //!< The extinction rate.
        const TypedDagNode< double > *                        psi;                                   //!< The sampling probability of a just extinct species.
        const TypedDagNode< double > *                        omega;                                 //!< The occurrence sampling rate.
        const TypedDagNode< double > *                        rho;                                   //!< The sampling probability of extant taxa.
        const TypedDagNode< double > *                        removalPr;                             //!< The removal probability after sampling.
        const TypedDagNode< long > *                          maxHiddenLin;                          //!< The maximal number of hidden lineages.
        const std::string&                                    cond;                                  //!< Condition of the process ("time" or "survival")
        const std::vector<double> &                           time_points;                           //!< Times at which density is computed
        const bool                                            useOrigin;                             //!< Start the process at the origin (otherwise root)
        const bool                                            useMt;                                 //!< Forward traversal Mt algorithm (otherwise backward Lt)
        const TypedDagNode< Tree > *                          timeTree;                              //!< Facultative initial tree
        mutable size_t                                        extant;                                //!< Number of extant taxa
        MatrixReal                                            B;                                     //!< Output, likelihoods through time matrix
        
    };
    
}

#endif

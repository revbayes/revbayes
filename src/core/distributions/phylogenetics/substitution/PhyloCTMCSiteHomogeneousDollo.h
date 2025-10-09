#ifndef PhyloCTMCSiteHomogeneousDollo_H
#define PhyloCTMCSiteHomogeneousDollo_H

#include <cstddef>
#include <map>
#include <vector>

#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DiscreteTaxonData.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "PhyloCTMCSiteHomogeneousConditional.h"
#include "RbException.h"
#include "StandardState.h"
#include "TopologyNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;

    struct DolloAscertainmentBias {

        enum Coding { ALL                 = 0x00,
                      NOABSENCESITES      = 0x01,
                      VARIABLE            = 0x03,
                      INFORMATIVE         = 0x0F
                    };
    };

    class PhyloCTMCSiteHomogeneousDollo : public PhyloCTMCSiteHomogeneousConditional<StandardState> {

        public:
        PhyloCTMCSiteHomogeneousDollo(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, DolloAscertainmentBias::Coding cod = DolloAscertainmentBias::NOABSENCESITES, bool normalize = true);

        PhyloCTMCSiteHomogeneousDollo(const PhyloCTMCSiteHomogeneousDollo&);
        // public member functions
        PhyloCTMCSiteHomogeneousDollo*                          clone(void) const;

        virtual void                                            redrawValue(void);
        void                                                    setDeathRate(const TypedDagNode< double > *r);

        protected:

            void                                                computeRootLikelihood(size_t root, size_t l, size_t r);
            void                                                computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
            void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
            void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
            void                                                computeTipLikelihood(const TopologyNode &node, size_t nIdx);

            void                                                computeRootCorrection(size_t root, size_t l, size_t r);
            void                                                computeRootCorrection(size_t root, size_t l, size_t r, size_t m);
            void                                                computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
            void                                                computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
            void                                                computeTipCorrection(const TopologyNode &node, size_t nIdx);

            double                                              sumRootLikelihood( void );
            void                                                resizeLikelihoodVectors(void);
            void                                                updateTransitionProbabilities(size_t node_idx);
            void                                                getStationaryFrequencies( std::vector<std::vector<double> >& ) const;

            virtual double                                      computeIntegratedNodeCorrection(const std::vector<std::vector<std::vector<double> > >& partials, size_t nodeIndex, size_t mask, size_t mixture, const std::vector<double> &f);
            virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);

            size_t                                              dim;
            std::vector<double>                                 integrationFactors;
            std::vector< std::vector<size_t> >                  maskNodeObservationCounts;
            std::vector<double>                                 survival;
            size_t                                              activeMassOffset;
            size_t                                              massNodeOffset;

            bool                                                normalize;
            const TypedDagNode< double >*                       death_rate;

        private:
            double                                              getScaledNodeWeights(const TopologyNode &node, size_t pattern, std::vector<double>& weights);
            void                                                scale(size_t i);
            void                                                scale(size_t i, size_t l, size_t r);
            void                                                scale(size_t i, size_t l, size_t r, size_t m);
            virtual void                                        simulateDollo( const TopologyNode &node, std::vector<StandardState> &taxa, size_t rateIndex, std::map<size_t, size_t>& charCounts);
            void                                                updateTransitionProbabilityMatrices(void);
        };

}

#endif

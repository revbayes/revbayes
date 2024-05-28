#ifndef GENERATOR_TO_SITE_MODEL_H
#define GENERATOR_TO_SITE_MODEL_H

#include "SiteModel.h"
#include "RateGenerator.h"
#include "RateMatrix.h"

namespace RevBayesCore {

    class Tree;

    class GeneratorToSiteModel: public SiteModel
    {
        std::shared_ptr<const RateGenerator> generator;
        std::vector<double> frequencies;
        double _scale;

    public:
        GeneratorToSiteModel* clone() const;

	int getNumberOfStates() const;

        TransitionProbabilityMatrix calculateTransitionProbabilities(const Tree& tau, int node_index, double rate) const;
        bool simulateStochasticMapping(const Tree&, int node, int rate, std::vector<size_t>&, std::vector<double>&) const;

        std::vector<double> getRootFrequencies() const;

        std::optional<double> rate() const;
        void scale(double factor);
	// factor out rescalable class
        
        GeneratorToSiteModel& operator=(const GeneratorToSiteModel&) = default;;
        GeneratorToSiteModel& operator=(GeneratorToSiteModel&&) noexcept = default;

        GeneratorToSiteModel(const GeneratorToSiteModel&) = default;
        GeneratorToSiteModel(GeneratorToSiteModel&&) noexcept = default;
        
        GeneratorToSiteModel(const RateMatrix& m, double s = 1);
        GeneratorToSiteModel(const std::shared_ptr<const RateMatrix>& g, double s = 1);
        GeneratorToSiteModel(const RateGenerator& g, const std::vector<double>& freqs, double s = 1);
        GeneratorToSiteModel(const std::shared_ptr<const RateGenerator>& g, const std::vector<double>& freqs, double s = 1);
    };

}

#endif

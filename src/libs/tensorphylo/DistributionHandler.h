/*
 * DistributionHandler.h
 *
 *  Created on: Mar 9, 2020
 *      Author: xaviermeyer
 */

#ifndef INTERFACE_DISTRIBUTIONHANDLER_H_
#define INTERFACE_DISTRIBUTIONHANDLER_H_

#include <boost/config.hpp>
#include <string>
#include <vector>
#include <map>

namespace TensorPhylo {
namespace Interface {

typedef std::vector< double > stdVectorXd;
typedef std::vector< stdVectorXd > stdMatrixXd;
typedef std::map< std::vector<unsigned>, double > eventMap_t;

typedef size_t mapHistoriesKey_t;
typedef std::vector< std::pair<double, size_t> > mapHistoriesVal_t;
typedef std::map< mapHistoriesKey_t, mapHistoriesVal_t > mapHistories_t;
typedef std::vector< mapHistories_t > vecHistories_t;


typedef enum {
	TIME=0,
	ROOT_SURVIVAL=1,
	ROOT_MRCA=2,
	ROOT_SAMPLING=3,
	STEM_SURVIVAL=4,
	STEM_ONE_SAMPLE=5,
	STEM_TWO_EXT_SAMPLES=6,
	STEM_TWO_SAMPLES=7,
	ROOT_SAMPLING_AND_MRCA=8
} conditionalProbability_t;

typedef enum {
	AUTO_TUNING=0,
	SEQUENTIAL_OPTIMIZED=1,
	SEQUENTIAL_BRANCHWISE=2,
	PARALLEL_OPTIMIZED=3,
	PARALLEL_BRANCHWISE=4
} approximatorVersion_t;

typedef enum {
	EULER=0,
	RUNGE_KUTTA4=1,
	RUNGE_KUTTA54=2,
	RUNGE_KUTTA_DOPRI5=3
} integrationScheme_t;

typedef enum {
	DBG_NONE=0,
	DBG_PRINT=1,
	DBG_FILE=2
} debugMode_t;

typedef size_t version_t;

class BOOST_SYMBOL_VISIBLE DistributionHandler {
public:
	DistributionHandler() {};
	virtual ~DistributionHandler() {};

	virtual void setTree(const std::string &aNewickTree) = 0;
	virtual void setTree(const std::string &aNewickTree, const std::vector<bool> dirtyNodes) = 0;
	virtual void setData(const std::vector<std::string> &taxa, const std::map<std::string, std::vector<double> > &aProbabilityMap) = 0;

	virtual void forceSchedulerUpdate() = 0;
	virtual void forceApproximatorDirty() = 0;

	virtual void setApplyTreeLikCorrection(bool doApply) = 0;
	virtual void setConditionalProbCompatibilityMode(bool setActive) = 0;
	virtual void setQuasistationaryFrequencyMode(bool setActive) = 0;
	virtual void setLikelihoodApproximator(approximatorVersion_t approxVersion) = 0;
	virtual void setConditionalProbabilityType(conditionalProbability_t condProb) = 0;
	virtual void setIntegrationScheme(integrationScheme_t aIntScheme) = 0;
	virtual void setNumberOfThreads(size_t nThreads) = 0;
	virtual void setMaxNumStochMapTries(size_t aNumTries) = 0;
	virtual void setAbsoluteTolerance(double anAbsoluteTolerance) = 0;
	virtual void setRelativeTolerance(double aRelativeTolerance) = 0;

	virtual void setInitialDeltaT(double initDeltaT) = 0;

	virtual void setRootPrior(const stdVectorXd &rootPrior) = 0;

	virtual void setLambda(const stdVectorXd &times, const stdMatrixXd &lambdas) = 0;
	virtual void setMu(const stdVectorXd &times, const stdMatrixXd &mus) = 0;
	virtual void setPhi(const stdVectorXd &times, const stdMatrixXd &phis) = 0;
	virtual void setDelta(const stdVectorXd &times, const stdMatrixXd &deltas) = 0;

	virtual void setEta(const stdVectorXd &times, const std::vector< stdMatrixXd > &etas) = 0;
	virtual void setOmega(size_t aNState, const stdVectorXd &times, const std::vector< eventMap_t > &omegas) = 0;

	virtual void setMassSpeciationEvents(const stdVectorXd &massSpeciationTimes, const stdMatrixXd &massSpeciationProb) = 0;
	virtual void setMassExtinctionEvents(const stdVectorXd &massExtinctionTimes, const stdMatrixXd &massExtinctionProb) = 0;
	virtual void setMassExtinctionStateChangeProb(const std::vector< stdMatrixXd> &massExtinctionStateChangeProb) = 0;
	virtual void setMassSamplingEvents(const stdVectorXd &massSamplingTimes, const stdMatrixXd &massSamplingProb) = 0;
	virtual void setMassDestrSamplingEvents(const stdVectorXd &massDestrSamplingTimes, const stdMatrixXd &massDestrSamplingProb) = 0;

	virtual void setSyncMonitors(const std::vector< double > &synchMonitoring) = 0;

	virtual double computeLogLikelihood() = 0;
	virtual void keepSchedulerAndApproximator() = 0;
	virtual void resetSchedulerAndApproximator() = 0;

	virtual mapHistories_t drawHistory() = 0;
	virtual mapHistories_t drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges) = 0;
	virtual mapHistories_t drawAncestralStates() = 0;

	virtual vecHistories_t drawMultipleHistories(size_t nReplicas) = 0;
	virtual vecHistories_t drawMultipleAncestralStates(size_t nReplicas) = 0;

	virtual void setDebugMode(debugMode_t aDebugMode) = 0;
	virtual void setDebugMode(debugMode_t aDebugMode, const std::string &aFilePath) = 0;

	virtual void writeStateToFile(const std::string &aFilePath) = 0;
	virtual void loadStateFromFile(const std::string &aFilePath) = 0;

	virtual size_t getVersion() const = 0;

	virtual void setSeed(size_t aSeed) const = 0;

};

} /* namespace Interface */
} /* namespace TensorPhylo */

#endif /* INTERFACE_DISTRIBUTIONHANDLER_H_ */

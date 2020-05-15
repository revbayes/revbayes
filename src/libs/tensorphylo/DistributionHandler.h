/*
 * DistributionHandlerImpl.h
 *
 *  Created on: Mar 9, 2020
 *      Author: xaviermeyer
 */

#ifndef INTERFACE_DISTRIBUTIONHANDLERIMPL_H_
#define INTERFACE_DISTRIBUTIONHANDLERIMPL_H_


#include "DistributionHandler.h"

#include "Tensor/IncFwdTensor.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Parameters/IncFwdParameterContainer.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"
#include "Test/Utils.h"

#include <fstream>
#include <Eigen/Core>

#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION
#include <boost/dll/alias.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace TensorPhylo {
namespace Interface {

class DistributionHandlerImpl : public DistributionHandler {
public: // CTOR-DTOR

	DistributionHandlerImpl();
	~DistributionHandlerImpl();

public: // CREATOR

	static boost::shared_ptr<DistributionHandlerImpl> create() {
		return boost::shared_ptr<DistributionHandlerImpl>(
			new DistributionHandlerImpl());
	}

public: // Implementation

	void setTree(const std::string &aNewickTreeStr);
	void setData(const std::vector<std::string> &aTaxa, const std::map<std::string, std::vector<double> > &aProbabilityMap);

	void setApplyTreeLikCorrection(bool aDoApply);
	void setConditionalProbCompatibilityMode(bool setActive);
	void setLikelihoodApproximator(approximatorVersion_t aApproxVersion);
	void setConditionalProbabilityType(conditionalProbability_t aCondProb);
	void setNumberOfThreads(size_t aNThreads);

	void setInitialDeltaT(double aInitDeltaT);

	void setRootPrior(const stdVectorXd &aRootPrior);

	void setLambda(const stdVectorXd &aTimes, const stdMatrixXd &aLambda);
	void setMu(const stdVectorXd &aTimes, const stdMatrixXd &aMu);
	void setPhi(const stdVectorXd &aTimes, const stdMatrixXd &aPhi);
	void setDelta(const stdVectorXd &aTimes, const stdMatrixXd &aDelta);

	void setEta(const stdVectorXd &aTimes, const std::vector< stdMatrixXd > &aEta);
	void setOmega(size_t aNState, const stdVectorXd &aTimes, const std::vector< eventMap_t > &aOmega);

	void setMassSpeciationEvents(const stdVectorXd &massSpeciationTimes, const stdMatrixXd &massSpeciationProb);
	void setMassExtinctionEvents(const stdVectorXd &massExtinctionTimes, const stdMatrixXd &massExtinctionProb);
	void setMassExtinctionStateChangeProb(const std::vector< stdMatrixXd> &massExtinctionStateChangeProb);
	void setMassSamplingEvents(const stdVectorXd &massSamplingTimes, const stdMatrixXd &massSamplingProb);
	void setMassDestrSamplingEvents(const stdVectorXd &massDestrSamplingTimes, const stdMatrixXd &massDestrSamplingProb);

	void setSyncMonitors(const std::vector< double > &aSynchMonitoring);

	double computeLogLikelihood();
	mapHistories_t drawHistory();
	mapHistories_t drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges);
	mapHistories_t drawAncestralStates();
	vecHistories_t drawMultipleHistories(size_t nReplicas);
	vecHistories_t drawMultipleAncestralStates(size_t nReplicas);

	void setDebugMode(debugMode_t aDebugMode);
	void setDebugMode(debugMode_t aDebugMode, const std::string &aFilePath);

	void writeStateToFile(const std::string &aFilePath);
	void loadStateFromFile(const std::string &aFilePath);

	size_t getVersion() const;

	void setSeed(size_t aSeed) const;

private:

	typedef enum {NONE=0, UPDATE=1, RESET=2} schedulerOperation_t;
	bool dirtyAsyncRateShifts, dirtySyncEventTimes, dirtyApproximator;
	bool dirtyLambda, dirtyMu, dirtyDelta, dirtyPhi, dirtyEta, dirtyOmega;

	debugMode_t debugMode;
	std::string debugFilePath;
	std::ofstream debugFile;

	schedulerOperation_t schedulerOperation;

	size_t nThreads;
	double initDeltaT;

	bool applyTreeLikCorrection, compatibilityMode;
	approximatorVersion_t approxVersion;
	conditionalProbability_t condProbType;

	Eigen::VectorXd rootPrior;
	std::vector<double> synchMonitoring;

	Phylogeny::Data::ContainerSharedPtr ptrData;
	Parameters::AsyncContainerSharedPtr ptrAsyncParams;
	Parameters::SyncContainerSharedPtr ptrSyncParams;

	Tensor::ContainerSharedPtr ptrTensors;
	SynchronousEvents::ContainerSharedPtr ptrSyncEvents;

	Phylogeny::Structure::TreeSharedPtr ptrTree;
	Likelihood::Approximator::ApproximatorSharedPtr approximator;
	boost::shared_ptr<Likelihood::Scheduler::BaseScheduler> scheduler;

	void updateSchedulerOperationRateChange();
	void updateSchedulerOperationSyncTimeChange();
	void updateSchedulerOperationTimeChange();

	void updateParameters();
	void updateSyncEvents();

	bool copyVectorWithCheck(const stdVectorXd &inVector, stdVectorXd &outVector);
	bool copyMatrixWithCheck(const stdMatrixXd &inMatrix, std::vector< Eigen::VectorXd > &outMatrix);
	bool copyVectorOfRateMatrixWithCheck(const std::vector< stdMatrixXd>  &inVecMatrix, std::vector< Eigen::MatrixXd > &outVecMatrix);
	bool copyVectorOfProbMatrixWithCheck(const std::vector< stdMatrixXd>  &inVecMatrix, std::vector< Eigen::MatrixXd > &outVecMatrix);
	bool copyVectorOfSparseTensorWithCheck(size_t aN, const std::vector< eventMap_t> &inVecSpTensors, std::vector< eventMap_t> &outVecSpTensors);

	void printParamVectorDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< Eigen::VectorXd > &params);
	void printParamMatrixDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< Eigen::MatrixXd > &params);
	void printParamTensorDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< eventMap_t > &params);

	double debugLikelihoodEvaluation();
	void debugChoseOutputStream(const std::string &string);

	Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr createStochasticMappingApprox();

	// For testing purpose
	friend double Test::SequentialCPU::computeLikForInterfaceTests(int idApproximator,
			   	  Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
			      const std::string &parametersFile);

	friend void Test::SequentialCPU::runBenchmarkNew(size_t nLikApproximation,
				Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
				const std::string &parametersFile,
				const std::string &logFile);

	friend double Test::SequentialCPU::computeLikForInterface(
			Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
			const std::string &parametersFile);

	friend std::vector< TensorPhylo::Interface::mapHistories_t >
		Test::SequentialCPU::computeStochasticMappingForInterface(
				size_t nHistories,
				Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
				const std::string &parametersFile);
};

BOOST_DLL_ALIAS(
	TensorPhylo::Interface::DistributionHandlerImpl::create,
    createTensorPhyloDistributionHandler
);

} /* namespace Interface */
} /* namespace TensorPhylo */

#endif /* INTERFACE_DISTRIBUTIONHANDLERIMPL_H_ */

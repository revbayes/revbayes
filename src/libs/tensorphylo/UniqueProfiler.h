/*
 * UniqueProfiler.h
 *
 *  Created on: Nov 30, 2011
 *      Author: meyerx
 */
#ifndef UNIQUEPROFILER_H_
#define UNIQUEPROFILER_H_


#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>

#include "Singleton.h"
#include <iosfwd>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <map>

namespace Utils {
namespace Profiling {

#define PROFILING		0
#define NB_TIMER		40

#if PROFILING
#define _START_EVENT(profiler, event) 	profiler->startEvent(event);
#define _END_EVENT(profiler, event)		profiler->endEvent(event);
#define _REPORT_TO_FILE(profiler)		profiler->writeReportToFile();
#else
#define _START_EVENT(profiler, event) 	void();
#define _END_EVENT(profiler, event)		void();
#define _REPORT_TO_FILE(profiler)		void();
#endif

class UniqueProfiler : public Singleton<UniqueProfiler>{

	friend class Singleton<UniqueProfiler>;

private:
	// Constructor destructor
	UniqueProfiler();
	~UniqueProfiler();

public:

	void startEvent(std::string event);
	void endEvent(std::string event);
	double reportEventAvgPerCallDuration(std::string event) const;
	double reportEventTotalDuration(std::string event) const;

	std::string reportToString() const;
	void writeReportToFile() const;
	void writeReportToFile(std::string filename) const;


private:

	typedef struct {
		timeval startTime;
		timeval endTime;
		std::string event;
		bool running;
	} timer_t;

	typedef struct {
		std::string event;
		std::string parentEvent;
		boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::sum, boost::accumulators::tag::mean> > accTime;
	} event_tracker_t;

	typedef std::map<std::string, event_tracker_t> trackerMap_t;
	trackerMap_t trackerMap;

	unsigned short int level;
	timer_t timers[NB_TIMER];

	double duration();

};

} // Profiling
} // Utils

#endif /* UNIQUEPROFILER_H_ */

/**
 * \file UniqueProfiler.cpp
 * \brief
 * \author Xavier Meyer
 * \date Nov 30, 2011
 *
 */

#include "UniqueProfiler.h"

#include <cstdio>
#include <utility>

namespace Utils {
namespace Profiling {

UniqueProfiler::UniqueProfiler() {

	level = 0;
	for(int i=0; i<NB_TIMER; i++)
		timers[i].running = false;

}

UniqueProfiler::~UniqueProfiler(){
}

/*
void UniqueProfiler::startTime(){
	 startTimer(0);
}

void UniqueProfiler::endTime(){
	endTimer(0);
}

double UniqueProfiler::duration(){
	return duration(0);
}
*/


void UniqueProfiler::startEvent(std::string event){

	//cout << "[" << level << "] " <<  event << endl;

	// if this current timer is already used, start a new one
	if(timers[level].running)
		level++;

	// Save the event name
	timers[level].event = event;

	// Set timer to running state
	timers[level].running = true;

	// Get the start time
	gettimeofday(&timers[level].startTime, NULL);
}

void UniqueProfiler::endEvent(std::string event){

	// Only capture event if the current timer is running
	if(timers[level].running){

		if(timers[level].event.compare(event) != 0){
			std::cout << "Profiler error (not same event)!" << std::endl;
			std::cout << std::scientific <<  duration() << " ";
			for(int i=0; i<level; i++) {
				std::cout << timers[i].event.c_str() << ">";
			}
			std::cout << timers[level].event.c_str() << std::endl;
		}

		// Capture end time
		gettimeofday(&timers[level].endTime, NULL);

		std::string &strEvent = event;
		/*std::stringstream ss;
		for(size_t i=0; i<level; ++i) {
			ss << timers[i].event << "==>";
		}
		ss << event;
		std::string strEvent = ss.str();*/

		// Write to file
		trackerMap_t::iterator it = trackerMap.find(strEvent);
		if(it == trackerMap.end()) {
			event_tracker_t eventTracker;

			eventTracker.event = strEvent; //event;
			if(level > 0) {
				eventTracker.parentEvent = timers[level-1].event;
			}
			eventTracker.accTime(duration());
			std::pair< std::string, event_tracker_t > newTracker = std::make_pair(strEvent, eventTracker);
			trackerMap.insert(newTracker);
		} else {
			it->second.accTime(duration());
		}

		/*myFile << std::scientific <<  duration() << " ";
		for(int i=0; i<level; i++) {
			myFile << timers[i].event.c_str() << ">";
		}
		myFile << timers[level].event.c_str() << std::endl;*/

		timers[level].running = false;
		if(level > 0)
			level--;
	} else {
		std::cout << "Profiler error (no running event)!" << std::endl;
		std::cout << std::scientific <<  duration() << " ";
		for(int i=0; i<level; i++) {
			std::cout << timers[i].event.c_str() << ">";
		}
		std::cout << timers[level].event.c_str() << std::endl;
	}
}

std::string UniqueProfiler::reportToString() const {
	std::stringstream ss;
	for(trackerMap_t::const_iterator it=trackerMap.begin(); it != trackerMap.end(); ++it) {
		ss << "Event : '" << it->second.event;
		ss << "', total time = " << std::scientific << boost::accumulators::sum(it->second.accTime);
		ss << ", time per call = " << std::scientific << boost::accumulators::mean(it->second.accTime);
		ss << ", nb calls = " << std::scientific << boost::accumulators::count(it->second.accTime);
		if(!it->second.parentEvent.empty()) {
			ss << ", parent = '" << it->second.parentEvent << "'";
		}
		ss << std::endl;
	}

	return ss.str();
}

double UniqueProfiler::reportEventAvgPerCallDuration(std::string event) const {
	trackerMap_t::const_iterator it = trackerMap.find(event);
	if(it != trackerMap.end()) {
		return boost::accumulators::mean(it->second.accTime);
	} else {
		assert(false && "Unexisting event");
	}
}

double UniqueProfiler::reportEventTotalDuration(std::string event) const {
	trackerMap_t::const_iterator it = trackerMap.find(event);
	if(it != trackerMap.end()) {
		return boost::accumulators::sum(it->second.accTime);
	} else {
		assert(false && "Unexisting event");
	}
}

void UniqueProfiler::writeReportToFile() const {
	writeReportToFile("custProf.txt");
}

void UniqueProfiler::writeReportToFile(std::string filename) const {
	std::ofstream myFile(filename.c_str(), std::ios::out);
	myFile << reportToString();
	myFile.close();
}


double UniqueProfiler::duration(){
	std::int64_t s = timers[level].endTime.tv_sec - timers[level].startTime.tv_sec;
	std::int64_t us = timers[level].endTime.tv_usec - timers[level].startTime.tv_usec;
	return ((double)s + us*1e-6);
}

} // Profiling
} // Utils

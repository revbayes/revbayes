/*
 * Loader.cpp
 *
 *  Created on: Mar 12, 2020
 *      Author: xaviermeyer
 */

#include "Loader.h"
#include "RbException.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <filesystem>

#include <boost/dll.hpp>
#include <boost/function.hpp>
#include <boost/dll/import.hpp> // for import_alias
#include <boost/system/error_code.hpp>

using boost::system::error_code;

namespace Plugin {

const std::string Loader::DEFAULT_PLUGIN_PATH(".plugins");
const std::string Loader::PLUGIN_TENSORPHYLO_NAME("libTensorPhylo");

Loader::Loader() {
}

Loader::~Loader() {
}

bool Loader::isTensorPhyloLoaded() {
	return pluginTensorPhylo.is_loaded();
}

bool Loader::loadTensorPhylo() {

	const char * home = std::getenv ("HOME");
	if (home == NULL) {
	  throw RbException("Path to your home folder not found.\nProvide the full path to the folder containing the TensorPhylo library.");
	}
	std::string pluginPath(home);
	pluginPath += "/";
	pluginPath += DEFAULT_PLUGIN_PATH;

	return loadTensorPhylo(pluginPath);
}


bool Loader::loadTensorPhylo(const std::string &aPluginFolder) {

	// Checking for the plugin folder
	std::filesystem::path pluginPath(aPluginFolder);
	if(!std::filesystem::is_directory(pluginPath)) {
		throw RbException("The folder doesn't exist.");
		return false;
	}

	// Listing all files
	std::vector<std::filesystem::directory_entry> vecPlugins;
    copy(std::filesystem::directory_iterator(pluginPath), std::filesystem::directory_iterator(), std::back_inserter(vecPlugins));

    // Look for TensorPhylo
    bool found = false;
	std::filesystem::path tensorPhyloPath;
	for(std::vector<std::filesystem::directory_entry>::const_iterator it = vecPlugins.begin(); it != vecPlugins.end();  ++ it ) {
		std::string filePath(it->path().string());
		if(filePath.find(PLUGIN_TENSORPHYLO_NAME) != std::string::npos &&
			 (filePath.find(".so") != std::string::npos ||
		 	 filePath.find(".dylib") != std::string::npos ||
		 	 filePath.find(".dll") != std::string::npos)) {
			found =  true;
			tensorPhyloPath = filePath;
		}
	}

	if(!found) {
		throw RbException("Library libTensorPhylo not found in the given folder.");
		return false;
	}

	// try to load TensorPhylo
	try {
                pluginTensorPhylo.load(tensorPhyloPath, boost::dll::load_mode::append_decorations);
	} catch(const std::exception& e) {
		throw RbException("TensorPhylo failed to load.");
		return false;
	}

	return true;
}

TensorPhylo::DistributionHandlerSharedPtr Loader::createTensorPhyloLik() const {
	assert(pluginTensorPhylo.is_loaded() && "TensorPhylo must be loaded properly before trying to create a DistributionHandler (lik approximator).");

	typedef boost::shared_ptr<TensorPhylo::Interface::DistributionHandler> (TensorDistributionHandler_create_t)();
	boost::function<TensorDistributionHandler_create_t> creator;

	try {
		creator = boost::dll::import_alias<TensorDistributionHandler_create_t>(pluginTensorPhylo, "createTensorPhyloDistributionHandler");
	} catch(const std::exception& e) {
		std::cerr << "Something unexpected happened while trying to create a TensorPhylo likelihood: " << e.what() << std::endl;
	}

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandler> tpHandler = creator();

	return tpHandler;

}


} /* namespace Plugin */

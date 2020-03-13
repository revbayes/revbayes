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
#include <boost/filesystem.hpp>

#include <boost/dll.hpp>
#include <boost/function.hpp>
#include <boost/dll/import.hpp> // for import_alias
#include <boost/system/error_code.hpp>

using boost::system::error_code;

namespace Plugin {

const std::string Loader::DEFAULT_PLUGIN_PATH("./plugins");
const std::string Loader::PLUGIN_TENSORPHYLO_NAME("libTensorPhylo");

Loader::Loader() {
}

Loader::~Loader() {
}

bool Loader::isTensorPhyloLoaded() {
	return pluginTensorPhylo.is_loaded();
}

bool Loader::loadTensorPhylo() {
	return loadTensorPhylo(DEFAULT_PLUGIN_PATH);
}


bool Loader::loadTensorPhylo(const std::string &aPluginFolder) {

	// Checking for the plugin folder
	boost::filesystem::path pluginPath(aPluginFolder);
	if(!boost::filesystem::is_directory(pluginPath)) {
		std::cerr << "Path not found : " << pluginPath << std::endl;
		return false;
	}

	// Listing all files
	std::vector<boost::filesystem::directory_entry> vecPlugins;
    copy(boost::filesystem::directory_iterator(pluginPath), boost::filesystem::directory_iterator(), std::back_inserter(vecPlugins));

    // Look for TensorPhylo
    bool found = false;
	boost::filesystem::path tensorPhyloPath;
	for(std::vector<boost::filesystem::directory_entry>::const_iterator it = vecPlugins.begin(); it != vecPlugins.end();  ++ it ) {
		std::string filePath(it->path().string());
		if(filePath.find(PLUGIN_TENSORPHYLO_NAME) != std::string::npos) {
			found =  true;
			tensorPhyloPath = filePath;
		}
	}

	if(!found) {
		std::cerr << "Library libTensorPhylo not found in : " << pluginPath << std::endl;
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

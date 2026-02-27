/*
 * Loader.h
 *
 *  Created on: Mar 12, 2020
 *      Author: xaviermeyer
 */

#ifndef REVBAYES_LOADER_H_
#define REVBAYES_LOADER_H_

#define BOOST_DLL_USE_STD_FS

#include "DistributionHandler.h"

#include <string>
#include <memory>
#include <boost/dll/shared_library.hpp>
#include <filesystem>

namespace TensorPhylo {
	typedef std::shared_ptr<TensorPhylo::Interface::DistributionHandler> DistributionHandlerSharedPtr ;
}

namespace Plugin {

class Loader {
public:

    bool isTensorPhyloLoaded();

    bool loadTensorPhylo();
    bool loadTensorPhylo(const std::filesystem::path &aPluginFolder);

    TensorPhylo::DistributionHandlerSharedPtr createTensorPhyloLik() const;

private:

    friend Loader& loader() ;

    Loader();
    ~Loader();

    Loader(const Loader&) = delete;
    Loader& operator=(const Loader&) = delete;

    static const std::string DEFAULT_PLUGIN_PATH;
    static const std::string PLUGIN_TENSORPHYLO_NAME;


    boost::dll::shared_library pluginTensorPhylo;

};

inline Loader& loader() {
    static Loader instance;
    return instance;
}

} /* namespace Plugin */

#endif /* REVBAYES_LOADER_H_ */

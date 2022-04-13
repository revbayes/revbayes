#!/bin/bash

# keep track of base folder
root_folder=`pwd`
logfile=${root_folder}/log.txt

uname_out="$(uname -s)"
case "${uname_out}" in
    Linux*)     machine=Linux; compiler=g++;distr="";;
    Darwin*)    machine=Mac; compiler=clang++;distr="";;
    CYGWIN*)    machine=Windows-cygwin; toolchain="x86_64-w64-mingw32"; compiler="${toolchain}-g++";distr="";;
    *)          machine="UNKNOWN"
esac

# Checking 1) the OS, 2) the number of processors available for compilation
# and 3) if openmp is installed (macos only)
if [ "${machine}" == "UNKNOWN" ]; then
  echo "This type of OS is not supported. Follow the manual installation instructions."
  echo "[ERROR] ${machine} - ${compiler}" > ${logfile}
  exit
else
  echo "The installation will proceed for a '${machine}' system and will take several minutes."

  if [ "${machine}" == "Linux" ]; then
    foundGentooString=`uname -a | grep -i gentoo`
    if [ "${foundGentooString}" != "" ]; then
      distr="Gentoo"
    fi
  fi
  echo "${machine} - ${compiler} - ${distr}" > ${logfile}
fi

echo "########################################"
echo "Checking your installation"

# Prepare folders
tmp_folder=${root_folder}/tmp
mkdir -p ${tmp_folder}

# Check for mandatory requirements
# Check connection
con_check=`ping -q -c1 google.com &>/dev/null && echo 1 || echo 0`
if [ "${con_check}" == "0" ]; then
  echo "No connection available. Retry whith an internet connection."
  echo "No connection available. Retry whith an internet connection." >> ${logfile} 2>&1
  exit
else
  echo "> connection established successfully"
  echo "> connection established successfully" >> ${logfile} 2>&1
fi

# Check curl
curl_bin=`command -v curl`
if [ "${curl_bin}" == "" ]; then
  echo "'curl' is unavailable on your system and is required to download dependencies."
  echo "Please install curl and retry (e.g., install with homebrew for MacOS or aptitude for Linux)."
  echo "> no curl" >> ${logfile} 2>&1
  exit
else
  echo "> curl found: ${curl_bin}"
  echo "> curl found: ${curl_bin}" >> ${logfile} 2>&1
fi

# Check cmake
cmake_bin=`command -v cmake`
if [ "${cmake_bin}" == "" ]; then
  echo "'cmake' is unavailable on your system and is required to download dependencies."
  echo "Please install autoconf/automake and retry (e.g., install with homebrew for MacOS or aptitude for Linux)."
  echo "> no cmake" >> ${logfile} 2>&1
  exit
else
  echo "> cmake found: ${cmake_bin}"
  echo "> cmake found: ${cmake_bin}" >> ${logfile} 2>&1
fi

# Check make
make_bin=`command -v make`
if [ "${make_bin}" == "" ]; then
  echo "'make' is unavailable on your system and is required to download dependencies."
  echo "Please install make and retry (e.g., install with homebrew for MacOS or aptitude for Linux)."
  echo "> no make" >> ${logfile} 2>&1
  exit
else
  echo "> make found: ${make_bin}"
  echo "> make found: ${make_bin}" >> ${logfile} 2>&1
fi

# Check compiler
compiler_bin=`command -v ${compiler}`
if [ "${compiler_bin}" == "" ]; then
  echo "'${compiler}' is unavailable on your system and is required to download dependencies."
  echo "Please install ${compiler} and retry (e.g., install with homebrew for MacOS or aptitude for Linux)."
  echo "> no ${compiler}" >> ${logfile} 2>&1
  exit
else
  echo "> compiler found: ${compiler_bin}"
  echo "> compiler found: ${compiler_bin}" >> ${logfile} 2>&1
  # Check if c++11 ready
  isCompatible11=`${compiler_bin} checkCPPVersion.cpp -o checkCPPVersion`
  if [ "${isCompatible11}" == "0" ]; then
    echo "You need a compiler with c++11 enabled (e.g., g++ >=4.8 or Clang/llvm >= 3.3)."
    echo "You need a compiler with c++11 enabled (e.g., g++ >=4.8 or Clang/llvm >= 3.3)." >> ${logfile} 2>&1
    exit
  else
    echo "> compiler is c++11 compatible"
    echo "> compiler is c++11 compatible" >> ${logfile} 2>&1
  fi
fi


echo "########################################"
echo "Dowloading and installing boost locally"
echo "... will take a few minutes..."
if [ "${machine}" == "Windows-cygwin" ]; then
  boost_folder=/usr/${toolchain}/sys-root/mingw
  if [ ! -d ${boost_folder}/include/boost ]; then
    echo "'boost' is unavailable on your system and is required."
    echo "Please install boost and retry (i.e., install with cygwin - see instructions)."
    echo "> no boost in ${boost_folder}" >> ${logfile} 2>&1
    exit
  fi
else
  boost_build_folder=${tmp_folder}/boost_1_73_0
  boost_folder=${tmp_folder}/boost
  curl -fsSL https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.gz -o ${tmp_folder}/boost.tar.gz
  cd ${tmp_folder}
  tar -xzf boost.tar.gz
  echo "Compiling boost"
  echo "... will take a few minutes..."
  cd ${boost_build_folder}
  mkdir -p ${boost_folder}
  ./bootstrap.sh --prefix=${boost_folder}  >> ${logfile} 2>&1
  if [ "${distr}" == "Gentoo" ]; then
    ./b2 install --ignore-site-config  >> ${logfile} 2>&1
  else
    ./b2 install  >> ${logfile} 2>&1
  fi
  cd ${root_folder}
fi
echo "> boost available in: ${boost_folder}"
echo "> boost available in: ${boost_folder}" >> ${logfile} 2>&1
echo "> done"
echo "########################################"
echo ""

# Compiling revbayes
echo "############################"
echo "Compiling RevBayes"
if [ "${machine}" == "Windows-cygwin" ]; then
  echo "... will take a several minutes (most probably a few hours - hang tight) ..."
else
  echo "... will take a several minutes ..."
fi
echo "I am here!!!!"
export MY_OS=${machine}
export MY_BOOST_ROOT=${boost_folder}
echo "I am building!!!!"
bash build.sh ${build_args}  >> ${logfile} 2>&1
if [ -f "rb" ]; then
  echo "> RevBayes was sucessfully installed"
  echo "> RevBayes was sucessfully installed" >> ${logfile} 2>&1
else
  echo "> [WARNING] RevBayes was not sucessfully installed. [WARNING] <"
  echo "> RevBayes was not sucessfully installed." >> ${logfile} 2>&1
  exit
fi
echo "> done"
echo "############################"
echo ""

echo "############################"
echo "Cleanup..."
cd ${root_folder}
rm -r ${tmp_folder}
rm checkCPPVersion
echo "> done"
echo "############################"
echo ""
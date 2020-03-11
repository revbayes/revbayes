#!/bin/bash
set -ev

# Build
cd projects/cmake
mkdir build
./regenerate.sh -travis true -mpi ${USE_MPI} -help true -exec_name rb
cd build
cmake .
make -j4
cd ..

# Run tests
if [ ${TRAVIS_BUILD_STAGE_NAME} = "Test" ]
then
  cd ../..
  echo "\"Hello World\"" | projects/cmake/rb
  cd tests
  ./run_integration_tests.sh -mpi ${USE_MPI} ${TRAVIS_BUILD_DIR}/projects/cmake/rb
  # Run testiphy
  export PATH=${TRAVIS_BUILD_DIR}/projects/cmake:$PATH
  cd
  git clone https://gitlab.com/testiphy/testiphy.git
  cd testiphy
  ./testiphy rb
fi
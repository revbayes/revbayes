#!/bin/sh
set -e

all_args="$@"

#################
# command line options
# set default values
debug="false"
mpi="false"
cmd="false"
help2yml="false"
boost_root=""
boost_lib=""
boost_include=""
install_dir="$HOME/Applications/RevBayes"

meson_args=""
# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo 'Command line options are:
-h                              : print this help and exit.
-debug          <true|false>    : set to true to build in debug mode. Defaults to false.
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
-cmd            <true|false>    : set to true if you want to build RevStudio with GTK2+. Defaults to false.
-help2yml       <true|false>    : update the help database and build the YAML help generator. Defaults to false.
-boost_root     string          : specify directory containing Boost headers (e.g. `/usr/include`). Defaults to unset.
-boost_lib      string          : specify directory containing Boost libraries. (e.g. `/usr/lib`). Defaults to unset.
-boost_include  string          : specify directory containing Boost libraries. (e.g. `/usr/include`). Defaults to unset.
-install_dir    string          : defaults to ~/Applications/RevBayes
-j              integer         : the number of threads to use when compiling RevBayes. Defaults to automatic.

You can also specify meson options as -Dvar1=value1 -Dvar2=value2

Examples:
  ./build.sh -mpi true -help2yml true
  ./build.sh -boost_root /home/santa/installed-boost_1.72 -install_dir $HOME/local/RevBayes
  ./build.sh -boost_include /home/santa/boost_1_72_0/ -boost_lib /home/santa/boost_1_72_0/stage/lib
  ./build.sh -mpi true -Dhelp2yml=true
  export BOOST_ROOT=/home/santa/installed-boost_1.72; ./build.sh'
        exit
    fi

    case "$1" in -D*)
                     meson_args="$meson_args $1"
                     shift
                     continue
                     ;;
                 -*=*)
                     echo "$0: I don't understand '$1' - did you mean '$(echo $1 | sed 's/=/ /')'?"
                     exit 1
                     ;;
    esac

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done

have()
{
    # This works for many, but not all programs.
    ok=false;
    prog=$1;
    if $prog --version </dev/null >/dev/null 2>&1 ; then
        ok=true;
    elif $prog -v </dev/null >/dev/null 2>&1 ; then
        ok=true;
    elif $prog -V </dev/null >/dev/null 2>&1 ; then
        ok=true;
    fi
    echo "$ok"
}

create_venv()
{
    if [ -f meson_env/bin/activate ] ; then
        echo "Using python virtual environment 'meson_venv'"
    else
        echo "Creating python virtual environment 'meson_venv'"
        python3 -m venv meson_env
    fi
    . meson_env/bin/activate
}


if [ `have python3` = "false" ] ; then
    echo "You must have python3 installed to run meson."
    exit 1
fi

if [ `have ccache` = "false" ] ; then
    echo "Warning: program 'ccache' not found.  You probably want to install this to speed up compilation.";
    echo "You might try one of these methods:"
    if [ `have apt-get` = true ] ; then
        echo "  sudo apt-get install ccache"
    fi
    if [ `have brew` = true ] ; then
        echo "  brew install ccache"
    fi
fi

# Install meson + ninja if not already around
if [ `have meson` = "false" ] || [ `have ninja` = "false" ] ; then
    create_venv
fi

if [ `have meson` = "false" ] ; then
    echo "Can't find meson -- installing into 'meson_venv'"
    pip3 install meson
    echo
fi

if [ `have ninja` = "false" ] ; then
    echo "Can't find ninja -- installing into 'meson_venv'"
    pip3 install ninja
    echo
fi

if [ `have meson` = "false" ] ; then
    echo "Program 'meson' not found.  Please install it!";
    echo "You might try one of these methods:"
    if [ `have apt_get` = true ] ; then
        echo "  sudo apt-get install meson"
    fi
    if [ `have brew` = true ] ; then
        echo "  brew install meson"
    fi
    echo "  See  https://mesonbuild.com/Getting-meson.html"
    echo "You need python3 version 3.7 or higher to run meson"
    echo "  python version = $(python3 --version)"
    exit 1
fi

if [ `have ninja` = "false" ] ; then
    echo "Program 'ninja' not found.  Please install it!";
    echo "You might try one of these methods:"
    if [ `have apt_get` = true ] ; then
        echo "  sudo apt-get install ninja"
    fi
    if [ `have brew` = true ] ; then
        echo "  brew install ninja"
    fi
    exit 1
fi
if [ "$mpi" = "true" ] ; then
    BUILD_DIR="build-mpi"
else
    BUILD_DIR="build"
fi

if [ -z "${exec_name}" ] ; then
    if [ "$mpi" = "true" ] ; then
        exec_name="rb-mpi"
    else
        exec_name=rb
    fi
fi

if [ "$travis" = "true" ]; then
    BUILD_DIR="build"
    CC=${C_COMPILER}
    CXX=${CXX_COMPILER}
    exec_name=rb
    help2yml=true
fi

if [ "$debug" = "true" ] ; then
    meson_args="--buildtype=debug $meson_args"
fi

if [ "$mpi" = "true" ] ; then
    meson_args="-Dmpi=true $meson_args"
fi

if [ "$cmd" = "true" ] ; then
    meson_args="-Dstudio=true $meson_args"
fi

if [ -n "$boost_root" ] ; then
    export BOOST_ROOT="${boost_root}"
fi

if [ "$help2yml" = "true" ] ; then
    meson_args="-Dhelp2yml=true $meson_args"
fi

if [ -n "${install_dir}" ] ; then
    meson_args="-Dprefix=${install_dir} $meson_args"
fi

meson_args="-Drb-exe-name=${exec_name} $meson_args"

if [ "$1" = "clean" ]
then
    rm -rf "${BUILD_DIR}"
    exit 1
fi

######### Generate help database
if [ "$help2yml" = "true" ]
then
    (
        cd ../../src
        echo "Generating help database"
        perl ../help/md2help.pl ../help/md/*.md > core/help/RbHelpDatabase.cpp
    )
fi

######## Generate some files for meson
echo "Running './generate.sh' in $(pwd)"
./generate.sh
echo
echo


######### Print some environment variables
# * This can alert the user if some weird values have been set.
# * This also helps us replicate the call to cmake.
echo "Note these environment variables:"
for var in CC CXX CFLAGS CPPFLAGS CXXFLAGS LDFLAGS BOOST_ROOT BOOST_INCLUDEDIR BOOST_LIBRARYDIR ; do
    cmd="if [ -n \"\${$var}\" ] ; then echo \"  ${var}=\${$var}\"; fi"
    eval $cmd
done
echo

######### Do the build
ninja_args=install
if [ -n "$j" ] ; then ninja_args="-j $j ${ninja_args}" ; fi
echo "Trying to run 'ninja -C ${BUILD_DIR} ${ninja_args}' in $(pwd)"

if [ -e "${BUILD_DIR}/compile_commands.json" ] && (cd "${BUILD_DIR}" ; meson configure $meson_args) && ninja -C "${BUILD_DIR}" ${ninja_args} ; then
    true
else
    rm -rf "${BUILD_DIR}"
    mkdir "${BUILD_DIR}"
    ######### Actually run meson
    echo "Running 'meson ${BUILDDIR} ../.. $meson_args' in $(pwd)"
    meson ${BUILD_DIR} ../.. $meson_args
    echo
    echo "Running 'ninja -C ${BUILD_DIR} ${ninja_args}' in $(pwd)"
    ninja -C "${BUILD_DIR}" ${ninja_args}
fi




####### Restore GitVersion.cpp from backup
if [ -e  GitVersion_backup.cpp ] ; then
    cp GitVersion_backup.cpp ../../src/revlanguage/utils/GitVersion.cpp
    rm GitVersion_backup.cpp
fi

echo "RevBayes executable is at '${prefix}/bin/${exec_name}'"
echo "Please at the directory ${prefix}/bin to your PATH."

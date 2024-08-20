project('RevBayes', ['cpp', 'c'],
        version: '1.2.5-preview',
        default_options: [
          'buildtype=release',
          'cpp_std=c++17',
          'b_ndebug=if-release'
        ],
        meson_version: '>= 0.57',
        license: 'GPLv2')

cpp = meson.get_compiler('cpp')
cc = meson.get_compiler('c')

# For cross-compilation, `host_machine`  is where the binary will RUN
#                        `build_machine` is where the binary is   COMPILED.

if host_machine.cpu_family().startswith('x86_64')
  add_project_arguments(['-msse','-msse2','-msse3'], language : 'cpp')
  add_project_arguments(['-DRB_x86_64'], language: 'cpp')
elif host_machine.cpu_family().startswith('x86')
  add_project_arguments(['-msse','-msse2','-msse3'], language : 'cpp')
  add_project_arguments(['-DRB_x86'], language: 'cpp')
elif host_machine.cpu_family().startswith('arm') or host_machine.cpu_family().startswith('aarch')
  # Only add sse flags for Intel architecture -- not on (say) ARM
  add_project_arguments(['-DRB_ARM'], language: 'cpp')
endif

if get_option('buildtype') == 'debug' or get_option('buildtype') == 'debugoptimized'
  # FIXME: -Woverloaded-virtual is useful, but has too many false-positives with GCC?
  add_project_arguments(['-Wall',
                         # These bugs are very hard to track down:
                         '-Wstrict-aliasing'],
                        language : 'cpp')
endif

# Disable some warnings
add_project_arguments([# These just add noise:
                       '-Wno-sign-compare','-Wno-unknown-pragmas',
                       # These could be bugs, but there are so many we can't see the really bad ones. So, disable them for now:
                       '-Wno-reorder',
                       # Not necessarily a bug:
                       '-Wno-unused-variable','-Wno-unused-but-set-variable'],
                      language : 'cpp')

# Handle file too large error when cross-compiling to windows
if host_machine.system() == 'windows'
  add_project_arguments(['-Wa,-mbig-obj'], language : 'cpp')
endif

boost_modules = [ 'regex',
                  'program_options',
                  'thread',
                  'system',
                  'filesystem',
                  'date_time',
                  'serialization']

# This doesn't build an internal boost from the internal copy, yet.
if get_option('static_boost')
  boost = dependency('boost', modules : boost_modules, version: '>=1.71', static: true)
else
  boost = dependency('boost', modules : boost_modules, version: '>=1.71')
endif

if get_option('openlibm')
  openlibm = dependency('openlibm')
else
  openlibm = dependency('', required: false)
endif

rb_name = 'rb'
if get_option('mpi')
  add_project_arguments(['-DRB_MPI'], language: 'cpp')
  mpi = dependency('mpi', language: 'cpp')
else
  mpi = dependency('', required: false)
endif

if get_option('rb-exe-name') != 'default'
  rb_name = get_option('rb-exe-name')
endif

####
subdir('src')

##### rpath: where to find shared libraries #####
extra_rpath = ''
if target_machine.system() == 'linux'
  # This allows an installed structure of
  #
  # RevBayes/
  #    bin/
  #       rb
  #    lib/
  #       RevBayes/
  #          *.so
  # To find the shared libs from inside bin/ we use $ORIGIN/../lib/RevBayes/
  #
  extra_rpath = '$ORIGIN/../lib/RevBayes/'
elif target_machine.system() == 'darwin'
  # If ever we make an App bundle, the structure might be
  #
  # RevBayes.app/
  #    Contents/
  #       MacOS/
  #          rb
  #       PlugIns/
  #          *.dylib
  #
  # So, at that point, we might add "@loaderpath/../Plugins/"
  #
  extra_rpath = '@loaderpath/../lib/RevBayes/;@loaderpath/../Plugins/'
endif


############# libraries #################
libdl = cpp.find_library('dl')

core = static_library('rb-core',
                      core_sources,
                      include_directories: [src_inc],
                      dependencies: [boost,mpi])

revlanguage = static_library('rb-revlanguage',
                             revlanguage_sources,
                             include_directories: [src_inc],
                             dependencies: [boost,mpi])

libs = static_library('rb-lib',
                      libs_sources,
                      include_directories: [src_inc],
                      dependencies: [boost, libdl])

help2yml = static_library('rb-help2yml',
                          help2yml_sources,
                          include_directories: [src_inc],
                          dependencies: [boost])

############# executables #################

rb = executable(rb_name,
                ['src/revlanguage/main.cpp'],
                link_with: [core, revlanguage, libs],
                include_directories: [src_inc],
                dependencies: [boost, mpi, openlibm],
                install_rpath: extra_rpath,
                install: true)

if get_option('studio')
  gtk2 = dependency('gtk+-2.0')
  cmd = static_library('rb-cmd',
                       cmd_sources,
                       include_directories: [src_inc],
                       dependencies: [gtk2])

  #  Add cmd/main.cpp cmd/RbGTKGui.cpp
  #  includes: RbGTRGui.ui
  revstudio = executable('RevStudio',
                         ['src/cmd/main.cpp'],
                         link_with: [core, revlanguage, libs, cmd],
                         include_directories: [src_inc],
                         dependencies: [boost, mpi, gtk2, openlibm],
                         install_rpath: extra_rpath,
                         install: true)

endif

subdir('tests')

if get_option('help2yml')

  executable(rb_name+'help2yml',
             ['src/help2yml/main.cpp'],
             link_with: [core, revlanguage, libs, help2yml],
             include_directories: [src_inc],
             dependencies: [boost, mpi],
             install_rpath: extra_rpath,
             install: true)
endif

####

summary({'host': host_machine.system()
        }, section: 'Architecture')

summary({'exe name': rb_name,
         'mpi': get_option('mpi'),
         'studio': get_option('studio'),
        },section: 'Options')

summary({'prefix': get_option('prefix'),
        },section: 'Directories')

# 3.3 Enable/Disable extra debugging for dynamic programming
assertions_off = get_option('b_ndebug') == 'true' or (get_option('b_ndebug') == 'if-release' and
                                                      (get_option('buildtype') == 'release' or
                                                       get_option('buildtype') == 'plain')
                                                     )
assertions_enabled = 'enabled'
if assertions_off
  assertions_enabled = 'disabled'
endif

summary({'optimization': get_option('optimization'),
         'debug': get_option('debug'),
         'assertions': assertions_enabled,
         'openlibm': get_option('openlibm')
        },section: 'Configuration')

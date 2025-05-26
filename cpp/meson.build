project('bvg', 'cpp',
  default_options: [
    'buildtype=release',
    'cpp_std=c++20',
    'optimization=3',
    'cpp_args=-march=native -flto',
    'cpp_link_args=-flto',
    'warning_level=3',
  ]
)

# Add strict warnings for your code
add_project_arguments([
  '-Wpedantic',
  '-Wconversion',
  '-Wshadow',
  '-Wno-sign-conversion',  # Suppress size_t vs int warnings
  '-Wno-shorten-64-to-32', # Suppress 64-bit to 32-bit conversion warnings
], language: 'cpp')

# Suppress warnings from all Homebrew-installed libraries
add_project_arguments('-isystem', '/opt/homebrew/include', language: 'cpp')
add_project_arguments('-isystem', '/opt/homebrew/opt/tbb/include', language: 'cpp')

if host_machine.system() == 'darwin'
  add_project_arguments('-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST', language: 'cpp')
  add_project_arguments('-Wno-#warnings', language: 'cpp')
endif

boost_dep = dependency('Boost', method: 'cmake')
tbb_dep = dependency('tbb')

executable('bvg-cpp',
  'bvg.cpp',
  dependencies: [boost_dep, tbb_dep]
)

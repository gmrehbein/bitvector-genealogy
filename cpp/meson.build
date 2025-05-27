project('bvg', 'cpp',
  version: '1.0.0',
  default_options: [
    'buildtype=release',
    'cpp_std=c++20',
    'optimization=3',
    'warning_level=3',
    'b_ndebug=if-release'
  ]
)

# Platform detection
is_m1 = host_machine.cpu_family() == 'aarch64' and host_machine.system() == 'darwin'
is_darwin = host_machine.system() == 'darwin'

# Architecture-specific optimizations
arch_flags = []
if is_m1
  arch_flags += [
    '-march=armv8.4-a+crypto+dotprod+fp16+rcpc',
    '-mtune=apple-m1',
    '-fvectorize',
    '-fslp-vectorize',
    '-fno-math-errno',     # Safe math optimization
    '-fno-trapping-math',  # Safe math optimization
    '-DUSE_NEON=1'
  ]
  
  message('Building for Apple M1 with safe optimizations')
else
  arch_flags += ['-march=native']
  message('Building with -march=native')
endif

add_project_arguments(arch_flags, language: 'cpp')

# Strict warnings for your code
warning_flags = [
  '-Wpedantic',
  '-Wconversion', 
  '-Wshadow',
  '-Wno-sign-conversion',
  '-Wno-shorten-64-to-32'
]

# macOS-specific flags
if is_darwin
  warning_flags += [
    '-Wno-#warnings',
    '-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST'
  ]
endif

add_project_arguments(warning_flags, language: 'cpp')

# Suppress warnings from Homebrew libraries
homebrew_includes = [
  '/opt/homebrew/include',
  '/opt/homebrew/opt/tbb/include',
  '/opt/homebrew/opt/boost/include'
]

foreach inc : homebrew_includes
  add_project_arguments('-isystem', inc, language: 'cpp')
endforeach

# Dependencies
boost_dep = dependency('Boost', 
  method: 'cmake',
  modules: ['system', 'graph'],
  required: true
)

tbb_dep = dependency('tbb', required: true)
thread_dep = dependency('threads')

# Main executable
bvg_exe = executable('bvg-cpp',
  'bvg.cpp',
  dependencies: [boost_dep, tbb_dep, thread_dep],
  install: true
)

# Print build configuration
summary({
  'Architecture': host_machine.cpu_family(),
  'System': host_machine.system(),
  'M1 Optimizations': is_m1,
  'Build Type': get_option('buildtype'),
  'Optimization Level': get_option('optimization'),
}, section: 'Build Configuration')

summary({
  'Boost': boost_dep.version(),
  'TBB': tbb_dep.found() ? 'found' : 'not found',
  'Threads': thread_dep.found() ? 'found' : 'not found',
}, section: 'Dependencies')
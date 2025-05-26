project('bvg', 'cpp',
  default_options: [
    'buildtype=release',
    'cpp_std=c++20',
    'optimization=3',
    'buildtype=release',
    'cpp_args=-march=native -flto',
    'cpp_link_args=-flto'
  ]
)


boost_dep = dependency('Boost', method: 'cmake')

tbb_dep = dependency('tbb')

executable('bvg-cpp',
  'bvg.cpp',
  dependencies: [boost_dep, tbb_dep]
)

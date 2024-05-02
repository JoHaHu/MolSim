MolSim2024 - Group D
===

## Team members

- Johannes Hupe
- Julius Kramer
- Tim Scholl


## Prerequisites

Developed and tested with these versions, other versions might work

- clang 17.0.6
- lldb 17.0.6
- cmake 3.27.7
- ninja 1.11.1

These environment variables are assumed to be always set. They are used to select which environment should be build. 
There purpose is simplifying the documentation only.

```shell
RELEASE_BUILD_DIR=build-release
DEBUG_BUILD_DIR=build-debug

BUILD_DIR=RELEASE_BUILD_DIR # or DEBUG_BUILD_DIR if debugging
```

Run this in the root folder of this repository

```shell
/usr/bin/cmake  -DCMAKE_BUILD_TYPE=Release \
                -DCMAKE_MAKE_PROGRAM=/usr/bin/ninja \
                -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
                -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
                -G Ninja \
                -S $(pwd) \
                -B $(pwd)/${BUILD_DIR}
```

`-G 'Unix Makefile'` would generate make files if required

## Building documentation

```shell
cd ${BUILD_DIR}
ninja doc_doxygen
```

## Build and run

```shell
cd ${BUILD_DIR}
ninja
./MolSim input t_end delta_t
```

## Clean

```shell
cd ${BUILD_DIR}
ninja clean
```

## Static analysis

```shell 
clang-tidy $files...$ -checks="cppcoreguidelines-*,modernize-*,performance-*,readability-*"  -p ./${BUILD_DIR}/compile_commands.json  
```

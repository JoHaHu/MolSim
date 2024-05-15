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

## Fetching git submodules

This project uses boost and spdlog as submodules. Make sure to fetch them before building the project.

Run this in the root folder of this repository:

```shell
git submodule update --init --recursive
```

Run this in the root folder of this repository

```shell
/usr/bin/cmake  -DCMAKE_BUILD_TYPE=Release \
                -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
                -G 'Unix Makefiles' \
                -S $(pwd) \
                -B $(pwd)/build
```

## Building documentation

```shell
cd ${BUILD_DIR}
make doc_doxygen
```

## Build and run

```shell
cd ${BUILD_DIR}
make
./MolSim --help
SPDLOG_LEVEL=info ./MolSim --input <input> --end_time <t_end> --delta_t <delta_t>
```

## Clean

```shell
cd ${BUILD_DIR}
make clean
```

## Static analysis

```shell 
clang-tidy $files...$ -checks="cppcoreguidelines-*,modernize-*,performance-*,readability-*"  -p ./${BUILD_DIR}/compile_commands.json  
```

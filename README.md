MolSim2024 - Group D
===

## Team members

- Johannes Hupe
- Julius Kramer
- Tim Scholl

## Prerequisites

Developed and tested with these versions, other versions might work

- clang 18.1.6
- lld 18.1.6
- lldb 18.1.6
- cmake 3.27.7

this repository includes `.idea/cmake.xml`.
This is a shareable IDE config file that allows all team members to use the same config.
The way Clion implements it, still allows to use your own config.
The only necessary prerequisite is that there is one Toolchain named "System" in your CLion config otherwise you might
get error messages.

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
SPDLOG_LEVEL=info ./MolSim --input <input_xml_file>
```

## Tests
```shell
cd test
./Test
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

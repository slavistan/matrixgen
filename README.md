# Matrixgen

## Install `Matrixgen`

External dependencies for the `Matrixgen` are:

1. [CMakeshift](https://github.com/mbeutel/CMakeshift)
2. [Microsoft's GSL implementation](https://github.com/microsoft/GSL)
3. [Intel's Threading Building Blocks (TBB)](https://github.com/intel/tbb)
4. [Eigen](https://github.com/libigl/eigen)

In order to build the unittests a recent version of [doctest](https://github.com/onqtam/doctest) is required. To enable unittests set the `BUILD_TESTS` CMake-variable accordingly.

`Makeshift` makes use of C++20 features. Tested compilers include `GCC 9.2.1`

Using `Matrixgen` in another project can be done by way of a regular `sudo make install` which will copy the pertinent files into your system's *include*, *lib* and *share* directories. The base path defaults to */usr/local/* and may be changed via CMake's `CMAKE_INSTALL_PREFIX` variable in the *CMakeCache.txt* file inside `Matrixgen`'s build directory.

```shell
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
# Optionally set `CMAKE_BUILD_TYPE` to match your project's needs.
sudo make install
```

Alternatively, you may also export `Matrixgen`'s build directory without copying any files by settings CMake's `EXPORT_BUILD_DIR` variable in the *CMakeCache.txt* file.

```shell
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
sed -i 's/EXPORT_BUILD_DIR:BOOL=.*/EXPORT_BUILD_DIR:BOOL=ON/g' CMakeCache.txt
make
```

# Make `Matrixgen` available to your project

After successfully performing one of the above methods of installation you may use `Matrixgen` like any other library via CMake's `find_package(Matrixgen 1.0.0 REQUIRED)`. The CMake target is called `Matrixgen::matrixgen`.

For a minimum working example create a file `main.cpp` ...

```cpp
// main.cpp
#include <matrixgen/core>
#include <iostream>

int main() {
  using Matrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  const auto mat = matrixgen::create<Matrix_t>(2, 2,
      { 1, 3,
        3, 7 });
  std::cout << mat << std::endl;
}
```

and alongside it define your CMake configuration in a `CMakeLists.txt` in the same directory ...

```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
find_package(Matrixgen 1.0.0 REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE Matrixgen::matrixgen)
```

Then build and run the binary.

```sh
mkdir build
cd build
cmake ..
make && ./main
```

# How to get started

Set the CMake-variable `BUILD_EXAMPLES` to build the examples and check out their verbosely commented source code in `examples/`.


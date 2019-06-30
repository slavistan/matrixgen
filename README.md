# Matrixgen - Header-only library

## Make `Matrixgen` available to your project

Using `Matrixgen` in another project can be done by way of a regular `sudo make install` which will copy the pertinent files into your system's *include*, *lib* and *share* directories. The base path defaults to */usr/local/* and may be changed via CMake's `CMAKE_INSTALL_PREFIX` variable in the *CMakeCache.txt* file inside `Matrixgen`'s build directory.

```shell
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
sudo make install
```

Alternatively, you may also export `Matrixgen`'s build directory without copying any files by settings CMake's `EXPORT_BUILD_DIR` variable in the *CMakeCache.txt* file. This option is a CMake-exclusive feature.

```shell
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
sed -i 's/EXPORT_BUILD_DIR:BOOL=.*/EXPORT_BUILD_DIR:BOOL=ON/g' CMakeCache.txt
make
```

After successfully performing one of the above methods of installation you may use `Matrixgen` like any other library via CMake's `find_package(Matrixgen REQUIRED)`. Note that, in addition to it's regular library dependencies, `Matrixgen` requires [`CMakeshift`][cmakeshift-url].

[cmakeshift-url]: https://github.com/mbeutel/CMakeshift


# matrixgen - Header-only library

## Make `matrixgen` available to your project

Using `matrixgen` in another project can be done by way of a regular `sudo make install` which will copy the pertinent files into your system's *include*, *lib* and *share* directories. The base path defaults to */usr/local/* and may be changed via CMake's `CMAKE_INSTALL_PREFIX` variable in the *CMakeCache.txt* file inside `matrixgen`'s build directory.

```sh
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
sudo make install
```

Alternatively, you may also export `matrixgen`'s build directory without copying any files by settings CMake's `EXPORT_BUILD_DIR` variable in the *CMakeCache.txt* file. This option is a CMake-exclusive feature.

```sh
git clone https://github.com/slavistan/matrixgen.git
cd matrixgen
mkdir build
cd build
cmake ..
sed -i 's/EXPORT_BUILD_DIR:BOOL=.*/EXPORT_BUILD_DIR:BOOL=ON/g' CMakeCache.txt
make
```



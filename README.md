# ROBUSfT

ROBUSFT is a C++ library for monocular real-time 3D shape tracking of isometrically deforming objects. 



## Build

In order to build the library you can run the following commands:

```cmake
mkdir build && cd build
cmake ..
make
make install
```
## Usage

Two examples of the usage of this library are presented in the Folder "examples". 

* `Example1_RegularTemplate.cpp` this code shows how to set up a template with rectangular shape and regular mesh.

* `Example2_RegularTemplate.cpp` this code shows how to set up a template with nonrectangular shape and irregular mesh.

After building the files, the examples can run by executing the following the commands.

```cmake
# for the example with regular template
cd build
./Example1_RegularTemplate

# for the example with regular template
cd build
./Example1_IrregularTemplate
```



### Using ROBUSfT as third party

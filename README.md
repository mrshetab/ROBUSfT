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

# for the example with irregular template
cd build
./Example1_IrregularTemplate
```

The texture-maps for both examples can be found in the folder `Texturemaps`. The texturemap for the template with regular mesh is a Spiderman poster. As for the template with a nonrectangular shape and an irregular shape, a shoe sole is chosen. As done in this example, in these cases, the texturemap of the object should be placed in the center of a totally white image.


### Using ROBUSfT as third party

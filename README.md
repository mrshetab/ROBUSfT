# ROBUSfT

ROBUSFT is a ready-to-use C++ library for monocular real-time 3D shape tracking of isometrically deforming objects. 
ROBUSFT forms a template from a deforming object and using that, it infers the deformation of the object in an input image (or a series of input images in a video). Each frame is processed indivisually. This makes ROBUSfT wide-baseline and resistant against losing the object in the video.  
This library is prsented in two versions: CPU-GPU (the current library) and Fully CPU (can be found in this [link](https://popsift.readthedocs.io/)).
The current version of ROBUSfT is capable of tracking surfaces up to 30fps. This speed is 20fps for the fully CPU version. 

## Dependencies

ROBUSfT depends on:

* OpenCV 

* Eigen3 

* PopSift (an open-source implementation of the SIFT algorithm in CUDA). To check the installation procedure and dependencies please visit this link. 

If you don't have a suitable cuda-adapted hardware and wish to use a pure CPU version please visit our full CPU ROBUSfT in this [link](https://popsift.readthedocs.io/).  

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
cd build/examples
./Example1_RegularTemplate

# for the example with irregular template
cd build/examples
./Example2_IrregularTemplate
```

The texture-maps for both examples can be found in the folder `Texturemaps`. 
Regarding the template with regular mesh, a Spiderman poster is selected as the texturemap. In the code, the number of mesh points in X and Y directions are chosen as 6 and 10, respectively. Foretheremore, the actual width of the object in X direction is set as 0.192m which is the width of the poster when is printed on an A4 paper.

As for the template with a nonrectangular shape and an irregular mesh, a shoe sole is chosen. Generally in these cases with an object with nonrectangular texturemaps, the texturemap of the object should be placed in the center of a totally white image as it is done for the shoe sole. The member function build_template will detach the textured region and generate a triangular mesh inside of that. Inside the code, the number of nodes on the boundary and inside of the texturemap will be chosen by the user.

# ROBUSfT

ROBUSFT is a C++ library for monocular real-time 3D shape tracking of isometrically deforming objects. 



## Build

adfsdfdsfadsfsdaf f sdaf dsaf dsf df df

```cmake
# Find the package from the PopSiftConfig.cmake
# in <prefix>/lib/cmake/PopSift/. Under the namespace PopSift::
# it exposes the target popsift that allows you to compile
# and link with the library
find_package(PopSift CONFIG REQUIRED)
...
# suppose you want to try it out in a executable
add_executable(poptest yourfile.cpp)
# add link to the library
target_link_libraries(poptest PUBLIC PopSift::popsift)
```



## Usage





### Using ROBUSfT as third party

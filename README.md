# EM Planner

This branch using foxglove-studio as visualization tool.

## 1. compile method

```shell
$ mkdir build
$ cd build
$ conan install .. --build=missing
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=/home/linkx/linkx_project/foxglove-cpp/examples/build/Debug/generators/conan_toolchain.cmake -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_BUILD_TYPE=Debug
$ make -j 4
$ ./main
```

## 2. Demo

<center class="half">
   <img src = "./figures/foxglove-1.png" width="400" height="300"/>
   <img src = "./figures/foxglove-2.png" width="400" height="300"/>
</center>

<center class="half">
   <img src = "./figures/foxglove-3.png" width="400" height="300"/>
   <img src = "./figures/foxglove-4.png" width="400" height="300"/>
</center>

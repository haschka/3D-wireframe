# 3D-wireframe
A tool to render and view basic 3D wireframes.
(c) Thomas Haschka 2024 - MIT license available. 

### Installation:
You need sdl2 libraries and headers installed, as well as a c compilier.
On a debian based system you can accomplish this by:
```
sudo apt-get install build-essential libsdl2-dev
```
Compilation of the tool works as follows:
```
./build.sh
```
You might edit `build.sh` in order to compile the tool on your choice of operating system.

### Usage:

In order to render an object you need two basic files.
One contains the vertices or nodes, the other the edges or links.
For examples look at `cube-nodes` and `cube-links`.

You can run the program by execting for instance:
```
./simple-3d cube-nodes cube-links
```

Rotate around:
- x axes with `q` and `e` keys
- y axes with `a` and `d` keys
- z axes with `z` and `c` keys

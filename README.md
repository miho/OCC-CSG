# OCC-CSG

Simple but powerful CSG tool based on the OpenCascade CAD Kernel ([OCE](https://github.com/tpaviot/oce) edition)

This tool provides a simple command-line interface for applying boolean operations and transformations to 3D objects specified in either BREP/STEP or STL format. It contains a convenient STL importer for converting mesh based geometries to an eqivalent BREP representation (might cause performance issues for large STL files). Geometries can be exported as STL files as well.

## Sample

To get started quickly download a [binary release](https://github.com/miho/OCC-CSG/releases) (available for Linux, macOS and Windows on x64).

These three lines

    ./occ-csg --create box -5,-5,-5,5,5,5 box.stp
    ./occ-csg --create sphere 0,0,0,6.5 sphere.stp
    ./occ-csg --csg difference box.stp sphere.stp cut.stp
    
produce this hollow cube:

<img src="resources/img/sample.jpg" width="400px">

To convert this resolution independent geometry to a triangulated STL file use:

    ./occ-csg --convert cut.stp cut.stl 0.1

The number at the end specifies the tolerance for the STL triangulation. Smaller values lead to better approximations. If the value is ommitted a default value (currently 0.5) will be used. This is how the resulting STL might look like:

<img src="resources/img/sample-stl.jpg" width="400px">

## CLI

To get an overview over the CLI type `./occ-csg --help`:

```
------------------------------------------------------------
------        CSG based on the OCE CAD Kernel         ------
------                  Version 0.3                   ------
------ 2018 by Michael Hoffer (info@michaelhoffer.de) ------
------                www.mihosoft.eu                 ------
------------------------------------------------------------
> USAGE:

Help & Info:

 occ-csg --help
 occ-csg --version

Creating Primitives:

 occ-csg --create box x1,y1,z1,x2,y2,z2                    box.stp
 occ-csg --create sphere x1,y1,z1,r                        sphere.stp
 occ-csg --create cyl x1,y1,z1,r1,h                        cyl.stp
 occ-csg --create extrusion ex,ey,ez,x1,y1,z1,x2,y2,z2,... extrude.stp

Format Conversion:

 occ-csg --convert file1.stl file1.stp
 occ-csg --convert file1.stp file1.stl 0.1

Geometric Transformation:

 occ-csg --transform matrix    t1,t2,t3,...,t12 file1.stp file1-transformed.stp
 occ-csg --transform translate x,y,z            file1.stp file1-translated.stp
 occ-csg --transform scale     sx,sy,sz         file1.stp file1-scaled.stp

Boolean Operators, Constructive Solid Geometry (CSG):

 occ-csg --csg union file1.stp file2.stp file-out.stp
 occ-csg --csg difference file1.stp file2.stp file-out.stp
 occ-csg --csg intersection file1.stp file2.stp file-out.stp

Bounds:

 occ-csg --bounds file.stp
 occ-csg --bounds file.stp bounds.stp
```

## How to build OCC-CSG

### Requirements

- C+\+ Compiler with C+\+11 support (tested with Clang, GCC, MSVC+\+ 2017 on x64)
- CMake >= 3.1
- [OCE](https://github.com/tpaviot/oce) >= 0.18.3*

For using `occ-csg` as independent command-line tool it is recommended to compile OCE as static library. This will increase the `occ-csg` file size but ensures the tool can be used without carrying too many additional libraries around.

#### Bash (Linux/macOS/Cygwin/other Unix-like Shell)

    cd /path/to/project
    mkdir build && cd build
    cmake .. -DOCE_DIR=/path/to/oce
    make -j4
    
#### Windows (CMD)

	"C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\vsdevcmd" -arch=x64
    cd \path\to\project
    mkdir build
    cd build
    cmake .. -DOCE_DIR=\path\to\oce
    MSBuild .\occ-csg-prj.sln  /property:Configuration=Release /property:Platform=x64
    

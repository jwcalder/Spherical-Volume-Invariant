This code computes the spherical volume invariant and PCA on local neighborhoods using the boundary integral method described in the paper

"Computation of circular area and spherical volume invariants via boundary integrals",  Riley O'Neill, Pedro Angulo-Umana, Jeff Calder, Bo Hessburg, Peter Olver, Cheri Shakiban, and Katrina Yezzi-Woodley, preprint, 2019.

COMPILATION

The core subroutines are coded in C, and need to be compiled. We provide Matlab and Python interfaces for the C code. For Matlab, we use the MEX interface. If MEX is not configured, enter "mex -setup" in the Matlab command window. Then navigate within Matlab to the c_code/ folder and execute the mexmake.m script to compile the c code. For Python, we use an extension module. Run svi_setup.py with options "build_ext --inplace" to compile the code for python. From the command line in a UNIX environment the compilation command is

>> python3 svi_setup.py build_ext --inplace

From within an interactive Python shell it is

>> run svi_setup.py build_ext --inplace

IMPORTANT: When compiling the code for Python, you must comment out the lines

#include "mex.h"
#define MEX

at the start of svi_computations.c.

USAGE

Matlab: The code takes as input a triangulated mesh in the Matlab triangulation format

>> T = triangulation(F,V);

where F is the triangle list and V is vertex list. To compute and plot the spherical volume invariant, run the code below from the main directory

>> addpath('c_code');
>> load('meshes/dragon.mat');
>> S = svi(T,1); 
>> color_surf(T,S,1,[-15 15]);

This code is contained in svi_demo.m.

For PCA on local neighborhoods run from the main directory

>> addpath('c_code');
>> load('meshes/dragon.mat');
>> [S,K1,K2] = svipca(T,1);
>> color_surf(T,K1,0.6,[-15 15]); title('K1');
>> color_surf(T,K2,0.6,[-15 15]); title('K2');
>> color_surf(T,K1.*K2,0.6,[-15 15]); title('Gauss Curvature');

This code is contained in svipca_demo.m

Python: The script svi_demo.py gives an example in Python.

Important Note: It is important to provide a manifold mesh to the algorithm, otherwise the results can be unpredictable. The code does some basic checks, and can handle stray vertices not connected to a triangle, and meshes with boundaries (e.g., edges that belong to only one triangle). The code does not check for non-manifold geometry, such as edges that are shared by 3 or more triangles, and in such cases the method is not reliable, since the surface integral formulation described in the paper is no longer valid. To avoid issues like this, one can use a mesh produced by a reliable algorithm, like marching cubes, or run some standard mesh cleaning software to remove non-manifold geometry, which is available in many free software packages, such as MeshLab.


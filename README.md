# interactive_diffraction
This project is an open-source implementation of the paper : Interactive Diffraction from Biological Nanostructures by Dhillon et al. [2014].
It is distributed under the LGPL license.

# Principle

Given a measured microstructure (a height map) the program computes lookup tables that can be used in a fragment shader to render realistic diffraction effects at real-time framerates.

# Libraries

This program has been compiled with : 

* CUDA 7.0
* OpenCV 2.4.11
* OpenMP

# Implementation

A naive implementation of the algorithm described in the paper takes days to compute the lookup tables even at a low resolution.
This implementation can handle high resolution height maps (tested with 1024x1024) and high resolution look up tables (tested with 1001x1001).
It has several optimisations that make the computation faster :

* Computation of Fast Fourier Transforms on the GPU with CUDA
* Computation of the convolution by using the separability of the Gaussian kernel
* Computation of half of the sum by using the commutativity of the $D_{n}^{m}$ functions
* Parallel computation of the convolution and the integration with OpenMP

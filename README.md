# Interactive diffraction
This project is an open-source implementation of the paper : Interactive Diffraction from Biological Nanostructures by Dhillon et al. [2014].
It is distributed under the LGPL license.

## Version

Version 1.0

## Principle
 
Given a measured microstructure (a height map) the program computes lookup tables that can be used in a fragment shader to render realistic diffraction effects at real-time framerates. Please refer to the paper for more information.

## Libraries

This software has been compiled with : 

* CUDA 7.0
* OpenCV 2.4.11
* OpenMP


A "interactive_diffraction.pro" file is provided to compile it with QtCreator. Please update the libraries' paths to match your installation.

## Implementation

A naive implementation of the algorithm described in the paper takes days to compute the lookup tables even at a low resolution.
This implementation can handle high resolution height maps (tested with 1024x1024) and high resolution look up tables (tested with 1001x1001).
It has several optimisations that make the computation faster :

* Computation of Fast Fourier Transforms on the GPU with CUDA
* Computation of the convolution by using the separability of the Gaussian kernel
* Computation of half of the sum by using the commutativity of the $D_{n}^{m}$ functions
* Parallel computation of the convolution and the integration with OpenMP

## License

Interactive Diffraction. Author :  Antoine TOISOUL LE CANN. Copyright © 2016 Antoine TOISOUL LE CANN. All rights reserved.

Interactive Diffraction is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. PFM_ReadWrite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU Lesser General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.



/*
*     Interactive Diffraction
*
*     Author:  Antoine TOISOUL LE CANN
*
*     Copyright © 2016 Antoine TOISOUL LE CANN
*              All rights reserved
*
*
* Interactive Diffraction is free software: you can redistribute it and/or modify
*
* it under the terms of the GNU Lesser General Public License as published by
*
* the Free Software Foundation, either version 3 of the License, or
*
* (at your option) any later version.
*
* Interactive Diffraction is distributed in the hope that it will be useful,
*
* but WITHOUT ANY WARRANTY; without even the implied warranty of
*
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
*
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

 /**
  * This program is an implementation of the
  * Interactive Diffraction from Biological Nanostructures
  * by Dhillon et al. [2014].
  *
  * \file fourier.h
  * \brief Implementation of Fourier transforms with CUDA.
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#ifndef FOURIER_H
#define FOURIER_H

#define M_PI 3.14159265358979323846

#include <iostream>
#include <complex>
#include <cmath>
#include <time.h>
#include <limits>
#include <vector>

//OpenCV
#include <opencv2/core/core.hpp>
#include <opencv/highgui.h>

/*---- CUDA ----*/
#include <cuda_runtime.h>
#include <cuda.h>
#include <cufft.h>

/**
 * Compute the the 2D FFT using CUDA.
 * Use Square images to calculate Fourier transforms otherwise it does not work.
 * @brief cudaFFT2D_C2C
 * @param input
 * @param width
 * @param height
 * @param output
 */
void cudaFFT2D_C2C(cufftComplex *input, int width, int height, cufftComplex *output);


/**
 * Compute the inverse 2D FFT using CUDA.
 * Use Square images to calculate Fourier transforms otherwise it does not work.
 * @brief cudaInvFFT2D_C2C
 * @param input
 * @param width
 * @param height
 * @param output
 */
void cudaInvFFT2D_C2C(cufftComplex *input, int width, int height, cufftComplex *output);

/**
 * Display the module of the Fourier transform using a log scale.
 * fourierImage has to be a CV_32C1 image that contains the module of the Fourier transform.
 * @brief displayLogIntensity
 * @param fourierImage
 */
void displayFourierImageLog(cv::Mat &fourierImage);

/**
 * Reorder the frequencies of a Fourier transform so that the 0 frequency is at the center.
 * Apply it twice to get back to the original coordinate system.
 * @brief centerFrequencies
 * @param output
 */
void centerFrequencies2D(cufftComplex *FT, int width, int height, cufftComplex *centeredFT);

#endif // FOURIER_H

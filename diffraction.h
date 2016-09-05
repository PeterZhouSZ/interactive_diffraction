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
 *  This program is an implementation of the
 *  Interactive Diffraction from Biological Nanostructures
 *  by Dhillon et al. [2014].
 *
 * \file diffraction.h
 * \brief Implementation of the diffraction algorithm.
 * \author Antoine Toisoul
 * \date August, 24th, 2016
 **/

#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#define M_PI 3.14159265358979323846

//Define TIME to time the computation of a look up table
#define TIME

#define NUMBER_OF_WAVELENGTHS 81 //81 wavelengths between 0.38 microns and 0.78 microns at a sampling of 0.005 microns.
#define WAVELENGTH_SAMPLING 0.005 //81 wavelengths between 0.38 microns and 0.78 microns at a sampling of 0.005 microns.

#include "fourier.h"
#include "mathfunctions.h"
#include "coloursystem.h"

/*---- Standard library ----*/
#include <iostream>
#include <cmath>
#include <time.h>
#include <complex>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>
#include <opencv/highgui.h>

#include <omp.h>

/**
 * Example on how to compute the diffraction look up tables given a image of a grating.
 * @brief computeExample
 * @param filePath
 */
void computeExample(string filePath);

/**
 * Compute the diffraction look up tables with all the optimisations described in the README.
 * The sampling distances have to be in micrometers.
 * @brief computeLookUpTables
 * @param N : number of Taylor terms.
 * @param lookUpTableSize : size of each lookup table (a lookup table is a square image of lookUpTableSize*lookUpTableSize).
 * @param heightMap : grayscale image representing the height map.
 * @param samplingSizeMap : size of a pixel of the height map in micrometers.
 * @param powerSpectrumDistribution : value of the power spectrum distribution at each wavlength.
 * @param numberOfWavelengths : Number of wavelengths
 * @param samplingSizeWavelength : sampling distance between each wavelength.
 * @param Ip
 * @param lambdaMin : lower bound of the wavelengths (0.38 microns for CIE color matching function).
 * @param nonLinearSamplingPower : Power of the non linear sampling. MUST BE ODD. (e.g 1.0 or 3.0 or 5.0 depending on the sample).
 */
void computeLookUpTables(const int N, const int lookUpTableSize, cv::Mat& heightMap, const float samplingSizeMap,
                         float spectrum[],  const int numberOfWavelengths, const float samplingSizeWavelength,
                         cv::Mat Ip[], const float lambdaMin = 0.38, const float nonLinearSamplingPower = 1.0);

/**
 * Computes the 3 integrands for all values of (u,v,lambda).
 * @brief computeIntegrand_forAllUVLambda
 * @param p
 * @param heightMap
 * @param FnsDFT
 * @param FmsDFT
 * @param gaussianU
 * @param gaussianV
 * @param lambdaMin
 * @param powerSpectrumDistribution
 * @param integrand_forAllUVLambdaX
 * @param integrand_forAllUVLambdaY
 * @param integrand_forAllUVLambdaZ
 * @param lookUpTableSize
 * @param numberOfWavelengths
 * @param samplingSizeWavelength
 */
void computeIntegrand_forAllUVLambda(int p, cv::Mat& heightMap, cufftComplex *FnsDFT, cufftComplex *FmsDFT, std::vector<float *> &gaussianU, std::vector<float *> &gaussianV,
                               float lambdaMin, float powerSpectrumDistribution[], std::vector<float *> &integrand_forAllUVLambdaX, std::vector<float *> &integrand_forAllUVLambdaY, std::vector<float *> &integrand_forAllUVLambdaZ,
                                     int lookUpTableSize, int numberOfWavelengths, float samplingSizeWavelength);

/**
 * Given FinalResultFn and FinalResultFm that contains the values of Fn anf Fm for all u and a specific v=vFixed and a specific lambda=lambdaFixed, the function
 * computes the 3 integrands at the points (all u, vFixed, lambdaFixed).
 * @brief computeIntegrand_vLambdaFixed
 * @param p
 * @param finalResultFn
 * @param finalResultFm
 * @param lookUpTableSize
 * @param vIndex
 * @param currentLambda
 * @param lambdaIndex
 * @param powerSpectrumDistribution
 * @param integrand_forAllUVLambdaX
 * @param integrand_forAllUVLambdaY
 * @param integrand_forAllUVLambdaZ
 */
void computeIntegrand_vLambdaFixed(int p, std::complex<float> *finalResultFn, std::complex<float> *finalResultFm, int lookUpTableSize,
                int vIndex, float currentLambda, int lambdaIndex, float powerSpectrumDistribution[], std::vector<float *> &integrand_forAllUVLambdaX, std::vector<float *> &integrand_forAllUVLambdaY, std::vector<float *> &integrand_forAllUVLambdaZ);

/**
 * Given n and the height map computes the DFT needed for the computation of the Fn Function
 * @brief computeDFT_n
 * @param n
 * @param heightMap
 * @param DFTResult
 */
void computeDFT_n(const int n, cv::Mat& heightMap, cufftComplex *DFTResult);

/**
* Function to calculate the value of i^n*heightMap^n. The result is stored in cufftComplex to directly compute the FFT with CUDA.
* @brief complexPowerHeightMap
* @param n
* @param heightMap
* @param resultFn
* @param size
*/
void complexPowerHeightMap(int n, cv::Mat& heightMap, cufftComplex *resultFn);

#endif // DIFFRACTION_H


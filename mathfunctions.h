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
  * \file mathfunctions.h
  * \brief Implementation of several useful mathematical functions.
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#ifndef MATHFUNCTIONS
#define MATHFUNCTIONS

#define M_PI 3.14159265358979323846

#include <iostream>
#include <cmath>
#include <vector>

/*---- OpenCV ----*/
#include <opencv2/core/core.hpp>

#include <omp.h>

using namespace std;
using namespace cv;

/**
 * Function to calculate the integral of a sampled function. size is the number of elements in the array. The sampling is done every samplingDistance.
 * @brief numbericalIntegration
 * @param functionSamples
 * @param samplingDistance
 * @param inf
 * @param sup
 * @return
 */
float numericalIntegration(float functionSamples[], int size, float samplingDistance);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param colorChannel
 * @param integrand_forAllUVLambda
 * @param numberOfUV
 * @param numberOfWavelengths
 * @param samplingDistanceWavelength
 * @param integrationResult
 */
void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param colorChannel
 * @param integrand_forAllUVLambda
 * @param width
 * @param height
 * @param numberOfWavelengths
 * @param samplingDistanceWavelength
 * @param integrationResult
 */
void numericalIntegration_onLambda(int colorChannel, vector<float *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult);
/**
 * Numerical integration over the visible spectrum using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param colorChannel
 * @param integrand_forAllUVLambda
 * @param width
 * @param height
 * @param numberOfWavelengths
 * @param samplingDistanceWavelength
 * @param integrationResult
 */
void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int width, int height,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule.
 * @brief numericalIntegration_onLambda
 * @param colorChannel
 * @param integrand_forAllUVLambda
 * @param numberOfUV
 * @param numberOfWavelengths
 * @param samplingDistanceWavelength
 * @param integrationResult
 */
void numericalIntegration_onLambda(int colorChannel, vector<double *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult);

/**
 * Numerical integration over the visible spectrum using the trapezoidal rule. No optimisations
 * @brief numericalIntegration_onLambda_non_optimised
 * @param p
 * @param colorChannel
 * @param factN
 * @param factM
 * @param integrand_forAllUVLambda
 * @param numberOfUV
 * @param numberOfWavelengths
 * @param samplingDistance
 * @param Ip
 */
void numericalIntegration_onLambda_non_optimised(int p, int colorChannel, double factN, double factM, vector<float *> &integrand_forAllUVLambda, int numberOfUV,
                                   int numberOfWavelengths, float samplingDistance, Mat Ip[]);

/**
 * Function to compute the factorial of n. No overflow for n<100 (at least).
 * @brief factorial
 * @param n
 * @return
 */
double factorial(double n);

/**
 * Given a 2D discrete function described as a matrix, returns function2D^power
 * @brief pow2D
 * @param function2D
 * @param power
 * @return
 */
Mat pow2D(Mat& function2D, unsigned int power);

/**
 * Returns the sign of the value x.
 * @brief sign
 * @param x
 * @return
 */
float sign(float x);

/**
 * Precomputes the gaussian coefficients needed for the calculation of the convolution.
 * @brief preComputeGaussian
 * @param lookUpTableSize
 * @param numberOfCoefficients
 * @param samplingSizeMap
 * @param lambdaMin
 * @param numberOfWavelength
 * @param samplingDistanceWavelength
 * @param result
 */
void preComputeGaussian(const int lookUpTableSize, const int numberOfCoefficients, const float samplingSizeMap, const float lambdaMin,
             const int numberOfWavelength, const float samplingDistanceWavelength, vector<float *> &result, const float nonLinearSamplingPower);


#endif // MATHFUNCTIONS


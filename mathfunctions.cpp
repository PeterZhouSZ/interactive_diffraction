/*
*     Interactive Diffraction
*
*     Author:  Antoine TOISOUL LE CANN, Imperial College London
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
  * \file mathfunctions.cpp
  * \brief Implementation of several useful mathematical functions.
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#include "mathfunctions.h"

/**
 * Function to calculate the integral of a sampled function. size is the number of elements in the array. The sampling is done every samplingDistance.
 * @brief numbericalIntegration
 * @param functionSamples
 * @param samplingDistance
 * @param inf
 * @param sup
 * @return
 */
float numericalIntegration(float functionSamples[], int size, float samplingDistance)
{
    //Assumption (sup-inf)/samplingDistance) is an integer
    float integral = 0.0;

    //Size-1 trapezes
    for(int i = 0 ; i<(size-1) ; i++)
    {
        //Add the area of a trapeze
        integral += samplingDistance*(functionSamples[i+1]+functionSamples[i])/2.0;
    }

    return integral;
}

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
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<numberOfUV ; u++)
    {

        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

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
                                   int numberOfWavelengths, float samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<height ; u++)
    {

        for(int v = 0 ; v<width ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*width+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*width+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

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
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult)
{
   // double currentValue = 0.0;
    //double nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<height ; u++)
    {

        for(int v = 0 ; v<width ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                double currentValue = integrand_forAllUVLambda[w][u*width+v];
                double nextValue = integrand_forAllUVLambda[w+1][u*width+v];

                integrationResult.at<Vec3f>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

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
                                   int numberOfWavelengths, double samplingDistanceWavelength, Mat& integrationResult)
{
   // float currentValue = 0.0;
    //float nextValue = 0.0;

    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int u = 0 ; u<numberOfUV ; u++)
    {
        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                double currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                double nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                integrationResult.at<Vec3d>(u,v).val[colorChannel] += samplingDistanceWavelength*(currentValue+nextValue)/2.0;
            }
        }
    }
}

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
                                   int numberOfWavelengths, float samplingDistance, Mat Ip[])
{
    for(int u = 0 ; u<numberOfUV ; u++)
    {
        for(int v = 0 ; v<numberOfUV ; v++)
        {
            //Integration over lambda with the trapezoidal rule
            for(int w = 0 ; w<numberOfWavelengths-1 ; w++)
            {
                float currentValue = integrand_forAllUVLambda[w][u*numberOfUV+v];
                float nextValue = integrand_forAllUVLambda[w+1][u*numberOfUV+v];

                Ip[p].at<Vec3f>(u,v).val[colorChannel] += static_cast<float>(1.0/(factN*factM))*samplingDistance*(currentValue+nextValue)/2.0;
            }
        }
    }
}

/**
 * Function to compute the factorial of n. No overflow for n<100 (at least).
 * @brief factorial
 * @param n
 * @return
 */
double factorial(double n)
{
    if(n==0)
        return 1;
    else
        return n*factorial(n-1);
}

/**
 * Given a 2D discrete function described as a matrix, returns function2D^power
 * @brief pow2D
 * @param function2D
 * @param power
 * @return
 */
Mat pow2D(Mat& function2D, unsigned int power)
{
    int width = function2D.cols;
    int height = function2D.rows;
    Mat result = Mat::ones(width, height, CV_32FC1);

    function2D.convertTo(function2D, CV_32FC1);

    if(power >0)
    {
        result = function2D.clone();
        for(unsigned int k = 2 ; k<=power ; k++)
        {
           result = result.mul(function2D);
        }
    }

    return result;
}

/**
 * Returns the sign of the value x.
 * @brief sign
 * @param x
 * @return
 */
float sign(float x)
{
    if(x<0)
        return -1.0;
    else
        return 1.0;
}

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
void preComputeGaussian(const int lookUpTableSize, const int heightMapSize, const float samplingSizeMap, const float lambdaMin,
              const int numberOfWavelength, const float samplingDistanceWavelength, vector<float*>& result, const float nonLinearSamplingPower)
{
    int halfSizeUV = (lookUpTableSize-1)/2;

    //The height map represents a size of 65 micrometer^2 in real life
    //Change this value for a different height map
    float sigmaS = 65.0/4.0;
    float sigmaF = 1.0/(2.0*M_PI*sigmaS);

    //Sampling
    float x = 0.0;
    float u = 0.0;
    float ul = 0.0;

    float currentLambda = 0.0;

    for(int k = -halfSizeUV ; k<=halfSizeUV ; k++)
    {
        ul = static_cast<float>(k)/static_cast<float>(halfSizeUV);
        u = 2.0*pow(ul, nonLinearSamplingPower);

        for(int t = 0 ; t<heightMapSize ; t++)
        {
            for(int l = 0 ; l<numberOfWavelength ; l++)
            {
                currentLambda = lambdaMin + (float)l*samplingDistanceWavelength;
                x = u/currentLambda-(t-(float)heightMapSize/2.0)/(float)(heightMapSize*samplingSizeMap);

                result[l][(k+halfSizeUV)*heightMapSize+t] = static_cast<float>(exp(-x*x/(2.0*sigmaF*sigmaF)));
            }
        }
    }
}

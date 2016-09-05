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
  * \file fourier.cpp
  * \brief Implementation of Fourier transforms with CUDA.
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#include "fourier.h"

using namespace std;
using namespace cv;

/**
 * Compute the the 2D FFT using CUDA.
 * Use Square images to calculate Fourier transforms otherwise it does not work.
 * @brief cudaFFT2D_C2C
 * @param input
 * @param width
 * @param height
 * @param output
 */
void cudaFFT2D_C2C(cufftComplex *input, int width, int height, cufftComplex *output)
{
    //input and ouput variables are on the host
    cufftHandle plan;

    /* Create a 2D FFT plan. Complex to Complex*/
    cufftPlan2d(&plan, width, height, CUFFT_C2C);


    //Memory allocation on the GPU
    cufftComplex *input_device;
    cufftComplex *output_device;
    cudaError errorMalloc1 = cudaMalloc((void**)&input_device, sizeof(cufftComplex)*width*height);
    cudaError errorMalloc2 = cudaMalloc((void**)&output_device, sizeof(cufftComplex)*width*height);

    if(errorMalloc1 != cudaSuccess)
    {
        cout << "Memory allcation error" << cudaGetErrorString(errorMalloc1) << endl;
    }
    else if(errorMalloc2 != cudaSuccess)
    {
        cout << "Memory allcation error" << cudaGetErrorString(errorMalloc2) << endl;
    }

    //Copy to GPU memory

    cudaError error = cudaMemcpy(input_device, input, width*height*sizeof(cufftComplex), cudaMemcpyHostToDevice);

    if(error != cudaSuccess)
    {
        cout << "Error when sending data to memory device : " << cudaGetErrorString(error) << endl;
    }

    /* Transform the first signal in place. Complex to Complex*/
    cufftResult cudaFFTResult = cufftExecC2C(plan, input_device, output_device,  CUFFT_FORWARD);

    if(cudaFFTResult != CUFFT_SUCCESS)
    {
        cout << "Error when calculating the FFT : " << cudaFFTResult << endl;
    }

    //The result is on the GPU in the variable output_device : transfer it back to the host CPU.
    cufftComplex *host_complexResult = new cufftComplex[width*height];
    error = cudaMemcpy(host_complexResult, output_device, width*height* sizeof(cufftComplex), cudaMemcpyDeviceToHost);

    if(error != cudaSuccess)
    {
        cout << "Error when sending data to memory host : " << cudaGetErrorString(error) << endl;
    }



    /* Center the frequency in the result output*/
    centerFrequencies2D(host_complexResult, width, height, output);

    //Free memory on the GPU and the CPU
    cufftDestroy(plan);
    cudaFree(input_device);
    cudaFree(output_device);

    delete[] host_complexResult;
}

/**
 * Compute the inverse 2D FFT using CUDA.
 * Use Square images to calculate Fourier transforms otherwise it does not work.
 * @brief cudaInvFFT2D_C2C
 * @param input
 * @param width
 * @param height
 * @param output
 */
void cudaInvFFT2D_C2C(cufftComplex *input, int width, int height, cufftComplex *output)
{
    cufftHandle plan;

    /* Create a 2D FFT plan. Complex to Complex*/
    cufftPlan2d(&plan, width, height, CUFFT_C2C);

    cufftComplex *inputDevice;
    cufftComplex *outputDevice;

    cudaError errorMalloc1 = cudaMalloc((void**)&inputDevice, sizeof(cufftComplex)*width*height);
    cudaError errorMalloc2 = cudaMalloc((void**)&outputDevice, sizeof(cufftComplex)*width*height);

    if(errorMalloc1 != cudaSuccess)
    {
        cout << "Memory allcation error" << cudaGetErrorString(errorMalloc1) << endl;
    }
    else if(errorMalloc2 != cudaSuccess)
    {
        cout << "Memory allcation error" << cudaGetErrorString(errorMalloc2) << endl;
    }

    //Copy to GPU memory
    cudaError error = cudaMemcpy(inputDevice, input, width*height*sizeof(cufftComplex), cudaMemcpyHostToDevice);

    if(error != cudaSuccess)
    {
        cout << "Error when sending data to memory device : " << cudaGetErrorString(error) << endl;
    }

    /* Calculate the inverse FFT on the device */
    cufftResult cudaFFTResult = cufftExecC2C(plan, inputDevice, outputDevice, CUFFT_INVERSE);

    if(cudaFFTResult != CUFFT_SUCCESS)
    {
        cout << "Error when calculating the FFT : " << cudaFFTResult << endl;
    }

    //Copy the memory back to the host

    cufftComplex *host_complexResult = new cufftComplex[width*height];
    error = cudaMemcpy(host_complexResult, outputDevice, width*height* sizeof(cufftComplex), cudaMemcpyDeviceToHost);

    if(error != cudaSuccess)
    {
        cout << "Error when sending data to memory host : " << cudaGetErrorString(error) << endl;
    }

    centerFrequencies2D(host_complexResult, width, height, output);

    //Free memory on the GPU and the CPU
    cufftDestroy(plan);
    cudaFree(inputDevice);
    cudaFree(outputDevice);

    delete[] host_complexResult;
}

/**
 * Display the module of the Fourier transform using a log scale.
 * fourierImage has to be a CV_32C1 image that contains the module of the Fourier transform.
 * @brief displayLogIntensity
 * @param fourierImage
 */
void displayFourierImageLog(Mat &fourierImage)
{
    //Calculate the maximum value
    float maximum = 0.0;
    int width = fourierImage.cols;
    int height = fourierImage.rows;

    for(int i = 0 ; i<height ; i++)
    {
        for(int j = 0 ; j<width ; j++)
        {
            if(fourierImage.at<float>(i,j)>maximum)
            {
                maximum = fourierImage.at<float>(i,j);
            }
        }
    }

    //Display the log scale
    Mat display = Mat::zeros(height, width, CV_32FC1);

    for(int i = 0 ; i<height ; i++)
    {
        for(int j = 0 ; j<width ; j++)
        {
            display.at<float>(i,j) = (float) floor(255.0*log(1.0+fourierImage.at<float>(i,j))/log(1.0+maximum));
        }
    }

    display.convertTo(display, CV_8UC1);
    imshow("Fourier transform", display);
    waitKey(0);
}


/**
 * Reorder the frequencies of a Fourier transform so that the 0 frequency is at the center.
 * Apply it twice to get back to the original coordinate system.
 * @brief centerFrequencies
 * @param output
 */
void centerFrequencies2D(cufftComplex *FT, int width, int height, cufftComplex *centeredFT)
{
    int iCentered = 0;
    int jCentered = 0;

    for (int i=0 ; i< height ; i++)
    {
       for (int j=0 ; j< width ; j++)
       {
           if(i<height/2 && j<width/2)
           {
               iCentered = i+height/2;
               jCentered = j+width/2;

               centeredFT[iCentered*width+jCentered].x = FT[i*width+j].x;
               centeredFT[iCentered*width+jCentered].y = FT[i*width+j].y;
           }
           else if(i<height/2 && j>=width/2)
           {             
               iCentered = i+height/2;
               jCentered = j-width/2;

               centeredFT[iCentered*width+jCentered].x = FT[i*width+j].x;
               centeredFT[iCentered*width+jCentered].y = FT[i*width+j].y;
           }
           else if(i>=height/2 && j<width/2)
           {
               iCentered = i-height/2;
               jCentered = j+width/2;

               centeredFT[iCentered*width+jCentered].x = FT[i*width+j].x;
               centeredFT[iCentered*width+jCentered].y = FT[i*width+j].y;
           }
           else if(i>=height/2 && j>=width/2)
           {
               iCentered = i-height/2;
               jCentered = j-width/2;

               centeredFT[iCentered*width+jCentered].x = FT[i*width+j].x;
               centeredFT[iCentered*width+jCentered].y = FT[i*width+j].y;
           }
       }
    }
}

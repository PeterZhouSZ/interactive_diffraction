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
  * \file diffraction.h
  * \brief Implementation of the diffraction algorithm.
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#include "diffraction.h"

using namespace std;
using namespace cv;

/**
 * Example on how to compute the diffraction look up tables given a image of a grating.
 * @brief computeExample
 * @param filePath
 */
void computeExample(string filePath)
{
    //Open a height map as a grayscale image
    Mat heightMap = imread(filePath.c_str(), CV_LOAD_IMAGE_GRAYSCALE);

    imshow("height map",heightMap);
    waitKey(0);

    //Put its value to the range [0;1]
    heightMap.convertTo(heightMap, CV_32FC1);
    heightMap /= 255.0;

    //N Taylor terms = 2N+1 look up tables
    int sizeTables = 1001;
    int N = 15;
    Mat *Ip = new Mat[2*N+1];

    //Initialisation
    for(int p = 0 ; p<=2*N ; p++)
    {
       Ip[p] = Mat::zeros(sizeTables,sizeTables, CV_32FC3);
    }


    //81 wavelengths between 0.38 microns and 0.78 microns at a sampling of 0.005 microns.
    float lambdaMin = (float) 0.38;
    float powerSpectrum[81];

    //Light with uniform white power spectrum
    for(int i=0 ; i<81 ; i++)
    {
        powerSpectrum[i] = 1.0;
    }

    float samplingDistanceInGrating = (float) 0.083;//micrometers
    computeLookUpTables(30, sizeTables, heightMap, samplingDistanceInGrating, powerSpectrum, 81, (float) 0.005, Ip, lambdaMin, (float)3.0);

}

/**
 * Compute the diffraction look up tables with all the optimisations described in the README.
 * The sampling distances have to be in micrometers.
 * @brief computeLookUpTables
 * @param N : number of Taylor terms.
 * @param lookUpTableSize : size of each lookup table (a lookup table is a square image of lookUpTableSize*lookUpTableSize).
 * @param heightMap : grayscale image representing the height map.
 * @param samplingSizeMap : size of a pixel of the height map in micrometers.
 * @param lambdaMin : lower bound of the wavelengths (0.38 microns for CIE color matching function).
 * @param powerSpectrumDistribution : value of the power spectrum distribution at each wavlength.
 * @param numberOfWavelengths : Number of wavelengths
 * @param samplingSizeWavelength : sampling distance between each wavelength.
 * @param Ip
 */
void computeLookUpTables(const int N, const int lookUpTableSize, Mat& heightMap, const float samplingSizeMap,
                         float powerSpectrumDistribution[],
                         const int numberOfWavelengths, const float samplingSizeWavelength, Mat Ip[], const float lambdaMin, const float nonLinearSamplingPower)
{

    //width height of the heightMap
    int width = heightMap.cols;
    int height = heightMap.rows;

    //Min and max of n in the sum
    int minimumValue = 0;//DO NOT USE unsigned int : p-N can be negative...
    int maximumValue = 0;
    int m = 0;

    /*----- I - Precompute the values of the Gaussians -----*/
    vector<float*> gaussiansU;
    vector<float*> gaussiansV;

    gaussiansU.resize(numberOfWavelengths);
    gaussiansV.resize(numberOfWavelengths);

    for(int i = 0 ; i<numberOfWavelengths ; i++)
    {
        gaussiansU[i] = new float[lookUpTableSize*height];
        gaussiansV[i] = new float[lookUpTableSize*width];
    }

    preComputeGaussian(lookUpTableSize, height, samplingSizeMap, lambdaMin, numberOfWavelengths, samplingSizeWavelength, gaussiansU, nonLinearSamplingPower);
    preComputeGaussian(lookUpTableSize, width, samplingSizeMap, lambdaMin, numberOfWavelengths, samplingSizeWavelength, gaussiansV, nonLinearSamplingPower);


    /*----- II - Main loop : Compute the 2N+1 look up tables -----*/
    //There are (lookUpTableSize)*(lookUpTableSize) possible values of (u,v) for each wavelength
    cufftComplex* FnsDFT;
    cufftComplex* FmsDFT;
    vector<float*> integrand_forAllUVLambdaX;
    vector<float*> integrand_forAllUVLambdaY;
    vector<float*> integrand_forAllUVLambdaZ;

    for(int p = 2 ; p<=2*N ; p++)
    {

#ifdef TIME
     clock_t start = clock();
     cout << "Value of p : " << p << endl;
#endif

         /*---- II-1 Find the lower and upper limits of the sum ----*/

         //Here p is fixed : the possible values of n and m can be known.
         //Hence the values of the DFT can already be computed
         //minimumValue = max(0, p-N);//DO NOT USE unsigned int p-N can be negative...
         //maximumValue = min(p, N);
        if(p<N)
        {
            if(p%2 == 0)
            {

                minimumValue = 0;
                maximumValue = p/2;
            }
            else
            {
                minimumValue = 0;
                maximumValue = (p-1)/2;

            }
        }
        else
        {
            if(p%2 == 0)
            {
                minimumValue = p-N;
                maximumValue = p/2;
            }
            else
            {
                 minimumValue = p-N;
                 maximumValue = (p-1)/2;
            }

        }


        /*---- II-2 Compute the look up table p (calculate the sum) ----*/
        //For each n generate the Fns and the Fms
        for(int n = minimumValue ; n<= maximumValue ; n++)
        {
            /*---- II-2-a Precompute the DFTs that are needed for the current term in the sum----*/
            FnsDFT = new cufftComplex[width*height];
            FmsDFT = new cufftComplex[width*height];

            m = p-n;

            //Compute Fn and Fm
            //Can be improved : if n==m only one DFT to compute
            computeDFT_n(n, heightMap, FnsDFT);
            computeDFT_n(m, heightMap, FmsDFT);


            /*---- II-2-b Compute the integrands for every u,v,lambda----*/
            //There are 3 integrands : one for each channel of the color space CIE XYZ.
            //Each integrand is a function of the variables (u,v,lambda)
            //The function (integrand) is stored in the 3 variables integrand_forAllUVLambdaX/Y/Z
            integrand_forAllUVLambdaX.resize(numberOfWavelengths);
            integrand_forAllUVLambdaY.resize(numberOfWavelengths);
            integrand_forAllUVLambdaZ.resize(numberOfWavelengths);

            for(int k = 0 ; k<numberOfWavelengths ; k++)
            {
                integrand_forAllUVLambdaX[k] = new float[(lookUpTableSize)*(lookUpTableSize)];
                integrand_forAllUVLambdaY[k] = new float[(lookUpTableSize)*(lookUpTableSize)];
                integrand_forAllUVLambdaZ[k] = new float[(lookUpTableSize)*(lookUpTableSize)];
            }

            computeIntegrand_forAllUVLambda(p, heightMap, FnsDFT, FmsDFT, gaussiansU, gaussiansV, lambdaMin, powerSpectrumDistribution, integrand_forAllUVLambdaX,
                                            integrand_forAllUVLambdaY, integrand_forAllUVLambdaZ,
                                            lookUpTableSize, numberOfWavelengths, samplingSizeWavelength);


            /*---- II-2-c Integrate each integrand over the variable lambda----*/
            //The result is stored in an image called integrationResult
            //The channels are stored as XYZ = RGB
            // /!\ Be careful : color channels follows OpenCV rules : 012 = BGR = ZYX
            //Integrate over lambda for every u,v
            Mat integrationResult = Mat::zeros(Ip[p].rows, Ip[p].cols, CV_32FC3);

            numericalIntegration_onLambda(2, integrand_forAllUVLambdaX, lookUpTableSize, numberOfWavelengths, samplingSizeWavelength, integrationResult);
            numericalIntegration_onLambda(1, integrand_forAllUVLambdaY, lookUpTableSize, numberOfWavelengths, samplingSizeWavelength, integrationResult);
            numericalIntegration_onLambda(0, integrand_forAllUVLambdaZ, lookUpTableSize, numberOfWavelengths, samplingSizeWavelength, integrationResult);


            /*---- II-2-d Add the current term to the sum (look up table) ----*/
            //Note that the sum is symmetric due to the commutativity of Dnm = Dmn.
            //As a result the integrand for (n,m) is the same as the one for (m,n)
            // /!\ Be careful if n=m

            //Example :
            //p = 0 D00
            //p = 1 D01 D10 = 2*D01
            //p = 2 D02 D11 D20 = D11+ 2D02
            //p = 3 D03 D12 D21 D30 = 2D03 + 2*D12
            //p = 4 D04 D13 D22 D31 D40 = D22 + 2D04 + 2D13
            //if p is even : then the sum is D(p/2,p/2) + 2*sum(D(i,p/2-i), i=0..p/2)
            //if p is odd : then the sum is 2*sum(D(i,(p-1)/2-i), i=0..(p-1)/2)

            //If N = 15, p= 16 : D115 D214 D313 D412 D511 D610 D79 | D88 | D97 D106 D115 D124 D133
            //If N = 16, p= 17 : D116 D215 D314 D413 D512 D611 D710 D89 | D98 D107 ...

            double factN = factorial(n);
            double factM = factorial(m);

            if(n == m)
            {
                Ip[p] += static_cast<float>(1.0/(factN*factM))*integrationResult;
            }
            else
            {
                Ip[p] += 2.0*static_cast<float>(1.0/(factN*factM))*integrationResult;
            }

            /*---- II-2-e Free the memory for the next look up table----*/
            for(int k = 0 ; k<numberOfWavelengths ; k++)
            {
                delete[] integrand_forAllUVLambdaX[k];
                delete[] integrand_forAllUVLambdaY[k];
                delete[] integrand_forAllUVLambdaZ[k];
            }



            delete[] FnsDFT;
            delete[] FmsDFT;
        }

        #ifdef TIME
            cout << "Time for Ip  " << p << " : " << (clock() - start) / (float)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
        #endif

        //Normalize by the area of the FT
        Ip[p] /= (float)width*(float)height;

        savePFM(Ip[p], string("I2.pfm"));

    }//End loop p

    //The look up tables are now stored in the array Ip

    for(int i = 0 ; i<numberOfWavelengths ; i++)
    {
        delete[] gaussiansU[i];
        delete[] gaussiansV[i];
    }

}

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
void computeIntegrand_forAllUVLambda(int p, Mat& heightMap, cufftComplex* FnsDFT, cufftComplex* FmsDFT,
                                     vector<float*>& gaussianU, vector<float*>& gaussianV,
                               float lambdaMin, float powerSpectrumDistribution[], vector<float*>& integrand_forAllUVLambdaX, vector<float*>& integrand_forAllUVLambdaY,
                                     vector<float*>& integrand_forAllUVLambdaZ, int lookUpTableSize,
                                     int numberOfWavelengths,float samplingSizeWavelength)
{

    //Notations are the same as the paper
    int width = heightMap.cols;
    int height = heightMap.rows;

    //Generate of uniform sampling of [-1 ; 1]x[-1 ; 1]

    //lookUpTableSize is the number of uniformly generated u,v
    int halfSizeUV = (lookUpTableSize-1)/2;
    //Variable to store the intermediate result of the convolution
    complex<float> *intermediateResultFn = new complex<float>[lookUpTableSize*width];
    complex<float> *intermediateResultFm = new complex<float>[lookUpTableSize*width];
    complex<float> *finalResultFn = new complex<float>[lookUpTableSize];
    complex<float> *finalResultFm = new complex<float>[lookUpTableSize];

    //For a given wavelength w
    float currentLambda = 0.0;

    for(int w = 0 ; w<numberOfWavelengths ; w++)
    {
        /*-- Initialise the Fn and Fms to zeros to compute the convolution--*/
        for(int i = 0 ; i<lookUpTableSize ; i++)
        {
            finalResultFn[i] = complex<float>();
            finalResultFm[i] = complex<float>();

            for(int j = 0 ; j<width ; j++)
            {
                intermediateResultFn[i*width+j] = complex<float>();
                intermediateResultFm[i*width+j] = complex<float>();
            }
        }

        //Sampling distance is in micrometers and lambda in micrometers
        currentLambda = lambdaMin+samplingSizeWavelength*(float)w;

        /*-- Compute the first step of the separable convolution on the rows--*/
        for(int k = -halfSizeUV ; k<=halfSizeUV ; k++)
        {
            //For (u,v), the sum is computed on the rows
            //OpenMP faster is each thread has its own const indices
            const int kConst = k+halfSizeUV;

            #pragma omp parallel for num_threads(omp_get_max_threads())
            for(int t = 0 ; t<width ; t++)
            {
                const int tConst = t;
                for(int s = 0 ; s<height ; s++)
                {
                    const int sConst = s;

                    //Filtering with the mask has to be calculated here
                    const float gaussian = gaussianU[w][kConst*height+sConst];

                    //intermediateResultFn is a lookUpTableSize*width matrix
                    intermediateResultFn[kConst*width+tConst] += complex<float>(FnsDFT[sConst*width+tConst].x*gaussian, FnsDFT[sConst*width+tConst].y*gaussian);
                    intermediateResultFm[kConst*width+tConst] += complex<float>(FmsDFT[sConst*width+tConst].x*gaussian, FmsDFT[sConst*width+tConst].y*gaussian);
                }
            }
        }

        //Repeat the procedure on the columns
        //fixes the value of v
        for(int v = -halfSizeUV ; v<=halfSizeUV ; v++)
        {
            /*-- Second step of computing the convolution for Fn and Fm for a specific v --*/
            //For (u,v), the sum is computed on the rows
            #pragma omp parallel for num_threads(omp_get_max_threads())
            for(int k = -halfSizeUV ; k<=halfSizeUV ; k++)
            {
                for(int t = 0 ; t<width ; t++)
                {
                    const float gaussian = gaussianV[w][(v+halfSizeUV)*width+t];

                    //FinalResultFn is a lookUpTableSize*1 matrix : contains the values of Fn(u/lambda,v/lambda) for the current v and all the u.
                    finalResultFn[k+halfSizeUV] += complex<float>(intermediateResultFn[(k+halfSizeUV)*width+t].real()*gaussian, intermediateResultFn[(k+halfSizeUV)*width+t].imag()*gaussian);
                    finalResultFm[k+halfSizeUV] += complex<float>(intermediateResultFm[(k+halfSizeUV)*width+t].real()*gaussian, intermediateResultFm[(k+halfSizeUV)*width+t].imag()*gaussian);
                }
            }

            /*-- Computing the integrands for this specific v --*/
            //Compute the integrand once Fn and Fm have been evaluated
            computeIntegrand_vLambdaFixed(p, finalResultFn, finalResultFm, lookUpTableSize, v+halfSizeUV, currentLambda, w, powerSpectrumDistribution,
                                            integrand_forAllUVLambdaX, integrand_forAllUVLambdaY, integrand_forAllUVLambdaZ);
            //Reinitialise for the next v
            for(int k = -halfSizeUV ; k<=halfSizeUV ; k++)
            {
                  finalResultFn[k+halfSizeUV] = complex<float>();
                  finalResultFm[k+halfSizeUV] = complex<float>();
            }
       }
    }//End of loop over wavelengths

    //Free the memory
    delete[] intermediateResultFn;
    delete[] intermediateResultFm;

    delete[] finalResultFn;
    delete[] finalResultFm;
}


/**
 * Given FinalResultFn and FinalResultFm that contains the values of Fn anf Fm for all u and a specific v=vFixed and a specific lambda=lambdaFixed, the function
 * computes the 3 integrands at the points (all u, vFixed, lambdaFixed).
 * @brief computeIntegrand_vLambdaFixed
 * @param p
 * @param finalResultFn
 * @param finalResultFm
 * @param numberOfUV
 * @param vIndex
 * @param currentLambda
 * @param lambdaIndex
 * @param powerSpectrumDistribution
 * @param integrand_forAllUVLambdaX
 * @param integrand_forAllUVLambdaY
 * @param integrand_forAllUVLambdaZ
 */
void computeIntegrand_vLambdaFixed(int p, complex<float> *finalResultFn, complex<float> *finalResultFm, int lookUpTableSize,
                int vIndex, float currentLambda, int lambdaIndex, float powerSpectrumDistribution[],
                             vector<float *> &integrand_forAllUVLambdaX,
                             vector<float *> &integrand_forAllUVLambdaY,
                             vector<float *> &integrand_forAllUVLambdaZ)
{

    //In this function v=vFixed and lambda=lambdaFixed are fixed
    //The values of Fn and Fm are known for all u and these v and lambda
    //Hence the integrands can be computed (for all u, vFixed, lambdaFixed)
    float Dnm = 0.0;
    float power = (float) pow(2.0*M_PI/currentLambda, p);

    //v is fixed : compute the value of the integrand for each u.
    for(int u = 0 ;u<lookUpTableSize ; u++)
    {

        //Compute Dnm
        Dnm =  finalResultFn[u].real()*finalResultFm[u].real()+finalResultFn[u].imag()*finalResultFm[u].imag();

        //Finish computing the integrand for each u
        integrand_forAllUVLambdaX[lambdaIndex][u*lookUpTableSize+vIndex] = power*Dnm*cie_colour_matching_function[lambdaIndex][0]*powerSpectrumDistribution[lambdaIndex];
        integrand_forAllUVLambdaY[lambdaIndex][u*lookUpTableSize+vIndex] = power*Dnm*cie_colour_matching_function[lambdaIndex][1]*powerSpectrumDistribution[lambdaIndex];
        integrand_forAllUVLambdaZ[lambdaIndex][u*lookUpTableSize+vIndex] = power*Dnm*cie_colour_matching_function[lambdaIndex][2]*powerSpectrumDistribution[lambdaIndex];
    }
}

/**
 * Given n and the height map computes the DFT needed for the computation of the Fn Function
 * @brief computeDFT_n
 * @param n
 * @param heightMap
 * @param DFTResult
 */
void computeDFT_n(const int n, Mat& heightMap, cufftComplex *DFTResult)
{
    int width = heightMap.cols;
    int height = heightMap.rows;
    Mat heightMapPower = Mat(heightMap.cols, heightMap.rows, CV_32FC1);

    /*-------------Compute the n-th power of the height map----------------*/
    heightMap.convertTo(heightMap, CV_32FC1);
    heightMapPower = pow2D(heightMap, n);

    /*-------------Compute i^n*heightMap^n----------------*/
    //Result stored in DFTparameter
    cufftComplex *DFTparameter = new cufftComplex[width*height];
    complexPowerHeightMap(n, heightMapPower, DFTparameter);

    /*-------------Compute the DFT with CUDA----------------*/
    cudaFFT2D_C2C(DFTparameter, width, height, DFTResult);

    delete[] DFTparameter;
}

/**
* Function to calculate the value of i^n*heightMap^n. The result is stored in cufftComplex to directly compute the FFT with CUDA.
* @brief complexPowerHeightMap
* @param n
* @param heightMap
* @param resultFn
* @param size
*/
void complexPowerHeightMap(int n, Mat& heightMap, cufftComplex *resultFn)
{
    int width = heightMap.cols;
    int height = heightMap.rows;

    //Compute i^n. It can take 4 values : 1, i, -1, -i (4 periodic)
    if(n%4 == 0)
    {
        //i^n = 1 : the parameter is real
        for(int i = 0 ; i<height ; i++)
        {
            for(int j = 0 ; j<width ; j++)
            {
                resultFn[i*width+j].x = heightMap.at<float>(i,j);
                resultFn[i*width+j].y = 0.0;
            }
        }
    }
    else if(n%4 == 1)
    {
        //i^n = i : the parameter is complex
        for(int i = 0 ; i<height ; i++)
        {
            for(int j = 0 ; j<width ; j++)
            {
                resultFn[i*width+j].x = 0.0;
                resultFn[i*width+j].y = heightMap.at<float>(i,j);
            }
        }
    }
    else if(n%4 == 2)
    {
        //i^n = -1 : the parameter is real
        for(int i = 0 ; i<height ; i++)
        {
            for(int j = 0 ; j<width ; j++)
            {
                resultFn[i*width+j].x = -heightMap.at<float>(i,j);
                resultFn[i*width+j].y = 0.0;
            }
        }
    }
    else
    {
        //i^n = -i : the parameter is complex
        for(int i = 0 ; i<height ; i++)
        {
            for(int j = 0 ; j<width ; j++)
            {
                resultFn[i*width+j].x = 0.0;
                resultFn[i*width+j].y = -heightMap.at<float>(i,j);
            }
        }
    }
}

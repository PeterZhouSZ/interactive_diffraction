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
  * \file main.cpp
  * \brief main.cpp
  * \author Antoine Toisoul
  * \date August, 24th, 2016
  **/

#include <iostream>

#include"diffraction.h"

using namespace std;

int main(int argc, char *argv[])
{

    computeExample(string("example.jpg"));

    return 0;
}


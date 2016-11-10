/*
GraDe-A: Grain Detection Algorithm.
Copyright (C) 2016 Paul Hoffrogge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ATOMDEFS_H_
#define ATOMDEFS_H_

#define NEAREST_ATOMNEIGHBORHOOD

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <vector>
#include <array>
#include <climits>
#include <string>

#define VOLUMEUNIT "1000nm^3"
#define DIM 3
#define MOORE 26
#define RADTODEG (57.29577951308232087679815481)
#define ONETHIRD (3.33333333333333333333333e-1)
#define SQRT2 	(1.41421356237309504880168872)
#define HALFSQRT2	(7.0710678118654752440084436e-1)
#define SQRT2_3		(8.16496580927726032732428e-1)
#define THIRD2 		(0.666666666666666666666666667)
#define THIRDSQRT2 (0.47140452079103168293389624140)
#define FOURTHIRDPI (4.188790204786390984616857844373)
#define ATOMALLOC 5
#ifdef _MSC_VER
#define DIRCHAR '\\'
#else
#define DIRCHAR '/'
#endif
#define BOX_ID_NAME "boxId"
#define ATOM_ID_NAME "atomId"
#define GRAIN_ID_NAME "grainId"
#define ORIENTATION_ID_NAME "oriId"
#define NO_ORIENTATION -1
#define NO_GRAIN -1
#define TIMEEVOSUBDIR "TimeEvo"
#define FANCYLINE 	"-############################################################################-"
#define LINE 		"------------------------------------------------------------------------------"

#include "OrientationMath.h"

typedef long oID;
typedef long gID;

typedef struct{
	long iB;
	long iA;
} AtomID;

namespace ori  {
	std::string to_string(double val, int precision = 10, bool forceScientific = false);
}
#endif /* ATOMDEFS_H_ */

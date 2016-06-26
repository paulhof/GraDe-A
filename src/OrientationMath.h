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

#ifndef SRC_ORIENTATIONMATH_H_
#define SRC_ORIENTATIONMATH_H_
#include <math.h>

#ifndef PI
#define PI 3.141592653589793238462643383
#endif

#ifndef HALFSQRT2
#define HALFSQRT2 7.0710678118654752440084436e-1
#endif

#ifndef SIGN
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)
#endif

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

#ifndef CUBE
#define CUBE(A) ((A)*(A)*(A))
#endif

#ifndef DIM
#define DIM 3
#endif

namespace ori{
	void inverseDirection(double * direct);
	void axisAngleToMatrix(const double * axis, double angle, double * M );
	void bunge2Quaternion( const double * euler, double * q );
	void quaternion2Bunge( const double * quat, double * euler );
	double misOrientation(const double * q1, const double * q2);
	void unitizeVectors(double * vects, unsigned long n);
	void crossProduct ( const double * v1, const double * v2, double * result);
	double scalarProduct(const double * v1, const double * v2);
	double radFromCosHalf(double cosHalf);
	double cosHalfFromRad(double radAngle);
	double cosHalfMisOrientation(const double * q1, const double * q2);
	void misOrientationQuaternion(const double  *q1, const double * q2, double * qMis);
	double cubicCosHalfMisOrientation(const double *q1, const double *q2);
	double cubicMisOrientation(const double * q1, const double * q2);
	bool haveCloseOrientations(const double * q1,const double * q2, double cosHalfThreshold);
	void bungeToMatrix(const double * euler, double * M );
	void rotationMatrixToQuaternion(const double * r, double * q);
	void uniqueCubicRotationQuaternion(const double *qIn, double *qOut);
}
#endif /* SRC_ORIENTATIONMATH_H_ */

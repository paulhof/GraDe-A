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

#include "OrientationMath.h"
#define QUAT2EUL_ETA 1e-20
double ori::radFromCosHalf(double cosHalf) {
	if(fabs(cosHalf) >= 1.){
		//1->2*0degree --> 0degree
		//-1->2*180degree --> 360degree-->0degree
		return 0.;
	}
	return 2*acos(cosHalf);
}

double ori::cosHalfFromRad(double radAngle){
	return cos(0.5*radAngle);
}

double ori::misOrientation(const double* q1, const double* q2){
	double cosHalfMis = cosHalfMisOrientation(q1, q2);
	return radFromCosHalf(cosHalfMis);
}

void ori::bungeToMatrix(const double * euler, double * M ){
	//See: Gottstein-Shvindlerman "Grain-Boundary-Migration in Metals" for reference
	double s1 = sin(euler[0]);
	double c1 = cos(euler[0]);
	double s2 = sin(euler[2]);
	double c2 = cos(euler[2]);
	double S = sin(euler[1]);
	double C = cos(euler[1]);
	//
	M[0] = c1*c2 - s1*s2*C;
	M[1] = s1*c2 + c1*s2*C;
	M[2] = s2*S;
	//
	M[3] = -c1*s2 - s1*c2*C;
	M[4] = -s1*s2 + c1*c2*C;
	M[5] = c2*S;
	//
	M[6] = s1*S;
	M[7] = -c1*S;
	M[8] = C;
}

void ori::axisAngleToMatrix(const double * axis, double angle, double * M ){
	//see Gottstein : "Physical Foundations of Materials Science"
	double unitAxis [3] = {axis[0],axis[1],axis[2]};
	ori::unitizeVectors(unitAxis,1);
	double c1 = (1-cos(angle));
	double c  = cos(angle);
	double s = sin(angle);
	//diag
	M[0] = SQR(unitAxis[0]) * c1 + c;//11
	M[4] = SQR(unitAxis[1]) * c1 + c;//22
	M[8] = SQR(unitAxis[2]) * c1 + c;//33
	//
	M[1] = unitAxis[0] * unitAxis[1] * c1 + unitAxis[2] * s; //12
	M[3] = M[1] - 2 * unitAxis[2] * s; //21
	//
	M[2] = unitAxis[0] * unitAxis[2] * c1 - unitAxis[1] * s; //13
	M[6] = M[2] + 2 * unitAxis[1] * s; //31
	//
	M[5] = unitAxis[1] * unitAxis[2] * c1 + unitAxis[0] * s; //23
	M[7] = M[5] - 2 * unitAxis[0] * s; //32
}

void ori::inverseDirection(double * direct){
	direct[0] *= -1;
	direct[1] *= -1;
	direct[2] *= -1;
}

void ori::unitizeVectors(double * vects, unsigned long n){
	double * v = vects;
	double l;
	for (unsigned long i = 0; i < n; i ++){
		l = sqrt(SQR(v[0]) + SQR(v[1]) + SQR(v[2]));
		v[0] /= l;
		v[1] /= l;
		v[2] /= l;
		v += DIM;
	}
}

void ori::crossProduct ( const double * v1, const double * v2, double * result){
	result [0] = v1[1] * v2[2] - v1[2] * v2[1];
	result [1] = v1[2] * v2[0] - v1[0] * v2[2];
	result [2] = v1[0] * v2[1] - v1[1] * v2[0];
}
double ori::scalarProduct(const double * v1, const double * v2){
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double ori::cosHalfMisOrientation(const double *q1, const double *q2){
	return fabs(q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3]);
}

bool ori::haveCloseOrientations(const double *q1, const double *q2, double cosHalfThreshold){
	//return cubicCosHalfMisOrientation(q1, q2) > cosHalfThreshold;
	return cosHalfMisOrientation(q1, q2) > cosHalfThreshold;
}

double ori::cubicMisOrientation(const double* q1, const double* q2) {
	double cosHalfMis = cubicCosHalfMisOrientation(q1, q2);
	//resolve numerical errors leading to a cosine value greater 1.
	return radFromCosHalf(cosHalfMis);
}

void ori::misOrientationQuaternion(const double  *q1, const double * q2, double * qMis) {
	//qMis = q1*q2^-1
	//00 + 11 + 22 + 33
	qMis[0] = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];
	//10 - 01 + 32 - 23
	qMis[1] = q1[1]*q2[0] - q1[0]*q2[1] + q1[3]*q2[2] - q1[2]*q2[3];
	//20 - 02 + 13 - 31
	qMis[2] = q1[2]*q2[0] - q1[0]*q2[2] + q1[1]*q2[3] - q1[3]*q2[1];
	//21 - 12 + 30 - 03
	qMis[3] = q1[2]*q2[1] - q1[1]*q2[2] + q1[3]*q2[0] - q1[0]*q2[3];
}

double ori::cubicCosHalfMisOrientation(const double *q1, const double *q2) {
	double qMis[4];
	misOrientationQuaternion(q1,q2,qMis);
	double a = qMis[0];
	double b = qMis[1];
	double c = qMis[2];
	double d = qMis[3];
	//all allowed (cosine half) rotation angles in cubic system with 4-fold symmetry
	double equivCosHalfMisorientations[24] =
	{
			//12 <100> 90 degree rotations
			//a+-b, c+-d
			(a+b)*HALFSQRT2,	(a-b)*HALFSQRT2,	(c+d)*HALFSQRT2,	(c-d)*HALFSQRT2,
			//a+-c, b+-d
			(a+c)*HALFSQRT2,	(a-c)*HALFSQRT2,	(b+d)*HALFSQRT2,	(b-d)*HALFSQRT2,
			//a+-d, b+-c
			(a+d)*HALFSQRT2,	(a-d)*HALFSQRT2,	(b+c)*HALFSQRT2,	(b-c)*HALFSQRT2,

			//8 <111> rotations
			(a+b+c+d)*0.5, (a+b-c-d)*0.5, (a-b+c-d)*0.5, (a-b-c+d)*0.5,
			(a+b+c-d)*0.5, (a+b-c+d)*0.5, (a-b+c+d)*0.5, (a-b-c-d)*0.5,

			//4 <100> 180degree rotations
			a,	b,	c,	d
	};
	//find maximum value
	double maxCosHalfMisori = 0.0;
	double absCurVal;
	for(char i=0; i < 24; i++){
		absCurVal = fabs(equivCosHalfMisorientations[i]);
		if (absCurVal > maxCosHalfMisori){
			maxCosHalfMisori = absCurVal;
		}
	}
	return maxCosHalfMisori;
}

void ori::rotationMatrixToQuaternion(const double * r, double * q) {
	if ( r[4] > -r[8] || r[0] > -r[4] || r[0] > -r[8]){
		q[0] = .5 * sqrt(1+r[0]+r[4]+r[8]);
		q[1] = .25 * (r[5] - r[7]) / q[0];
		q[2] = .25 * (r[6] - r[2]) / q[0];
		q[3] = .25 * (r[1] - r[3]) / q[0];
		return;
	}
	if ( r[4] < -r[8] || r[0] > r[4] || r[0] > r[8]){
		q[1] = .5 * sqrt(1+r[0]-r[4]-r[8]);
		q[0] = .25 * (r[5] - r[7]) / q[1];
		q[2] = .25 * (r[1] + r[3]) / q[1];
		q[3] = .25 * (r[6] + r[2]) / q[1];
		return;
	}
	if ( r[4] > r[8] || r[0] < r[4] || r[0] < -r[8]){
		q[2] = .5 * sqrt(1-r[0]+r[4]-r[8]);
		q[0] = .25 * (r[6] - r[2]) / q[2];
		q[1] = .25 * (r[1] + r[3]) / q[2];
		q[3] = .25 * (r[5] + r[7]) / q[2];
		return;
	}
	if ( r[4] < r[8] || r[0] < -r[4] || r[0] < r[8]){
		q[3] = .5 * sqrt(1-r[0]-r[4]+r[8]);
		q[0] = .25 * (r[1] - r[3]) / q[3];
		q[1] = .25 * (r[6] + r[2]) / q[3];
		q[2] = .25 * (r[5] + r[7]) / q[3];
		return;
	}
}

void ori::bunge2Quaternion( const double * euler, double * q )
{
	/*20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006
	Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	//cos(a+b) = c(a+b) = cacb-sasb
	//cos(a-b) = c(a-b) = cacb+sasb
	//sin(a+b) = s(a+b) = sacb+casb
	//sin(a-b) = s(a-b) = sacb-casb
	 */

	double p1 = euler[0];
	double t  = euler[1];
	double p2 = euler[2];

	double co1 = cos(t/2);
	double s1 = sin(t/2);

	double p[4] = {co1*cos((p1+p2)/2),s1*cos((p1-p2)/2),s1*sin((p1-p2)/2),co1*sin((p1+p2)/2)}; //applying sin, cos addition theorems

	//double test[4]={0};
	//quaternion2Euler( p,test);

	q[0] = p[0];
	q[1] = p[1];
	q[2] = p[2];
	q[3] = p[3];
}

void ori::quaternion2Bunge( const double * quat, double * euler )
{
	//convention: Bunge, ZXZ, equal to case (3,1,3) as
	//  analyzed in Diebel, James, 2006:
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401 - 413
	//dont forget to define QUAT2EUL_ETA 1e-20

	double q0 = quat[0];
	double q1 = quat[1];
	double q2 = quat[2];
	double q3 = quat[3];
	double PHI, sP, phi1, phi2;

	double cosPHI = SQR(q3)-SQR(q2)-SQR(q1)+SQR(q0);

	double y0 =	2*q1*q3	-	2*q0*q2;
	double x0 =	2*q2*q3	+	2*q0*q1;
	double y1 =	2*q1*q3	+	2*q0*q2;
	double x1 = -	2*q2*q3	+	2*q0*q1;

	if( cosPHI > 1. ) cosPHI = 1.;

	if( SQR( 1. - cosPHI ) <= QUAT2EUL_ETA )
		PHI = 0.;
	else
		PHI = acos( cosPHI );

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a quarternion does not define a Bunge Euler angle uniquely

	if( sP != 0 )
	{
		phi2 = atan2( y0 / sP, x0 / sP );
		phi1 = atan2( y1 / sP, x1 / sP );
	}else
	{
		phi1 = atan2( (2*q1*q2+2*q0*q3),SQR(q0)+SQR(q1)-SQR(q2)-SQR(q3) );
		phi2 = 0.;
	}

	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * PI;
	if (phi2 < 0.0)
		phi2 += 2 * PI;

	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
}

void ori::uniqueCubicRotationQuaternion(const double *qIn, double *qOut){
	double maxCos=0.0;
	double curAbsVal;
	char maxQuatID;
	char negator;
	//q1-q4
	curAbsVal = fabs(qIn[0]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 0;}
	curAbsVal = fabs(qIn[1]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 1;}
	curAbsVal = fabs(qIn[2]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 2;}
	curAbsVal = fabs(qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 3;}
	//q5-q8
	curAbsVal = 0.5 * fabs(qIn[0] - qIn[1] - qIn[2] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 4;}
	curAbsVal = 0.5 * fabs(qIn[0] + qIn[1] + qIn[2] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 5;}
	curAbsVal = 0.5 * fabs(qIn[0] - qIn[1] + qIn[2] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 6;}
	curAbsVal = 0.5 * fabs(qIn[0] + qIn[1] - qIn[2] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 7;}
	//q9-q12
	curAbsVal = 0.5 * fabs(qIn[0] + qIn[1] - qIn[2] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 8;}
	curAbsVal = 0.5 * fabs(qIn[0] - qIn[1] + qIn[2] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 9;}
	curAbsVal = 0.5 * fabs(qIn[0] + qIn[1] + qIn[2] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 10;}
	curAbsVal = 0.5 * fabs(qIn[0] - qIn[1] - qIn[2] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 11;}
	//q13-q16
	curAbsVal = HALFSQRT2 * fabs(qIn[0] - qIn[1]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 12;}
	curAbsVal = HALFSQRT2 * fabs(qIn[0] - qIn[2]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 13;}
	curAbsVal = HALFSQRT2 * fabs(qIn[0] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 14;}
	curAbsVal = HALFSQRT2 * fabs(qIn[1] + qIn[2]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 15;}
	//q17-q20
	curAbsVal = HALFSQRT2 * fabs(qIn[2] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 16;}
	curAbsVal = HALFSQRT2 * fabs(qIn[1] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 17;}
	curAbsVal = HALFSQRT2 * fabs(qIn[0] + qIn[1]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 18;}
	curAbsVal = HALFSQRT2 * fabs(qIn[0] + qIn[2]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 19;}
	//q21-q24
	curAbsVal = HALFSQRT2 * fabs(qIn[0] + qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 20;}
	curAbsVal = HALFSQRT2 * fabs(qIn[1] - qIn[2]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 21;}
	curAbsVal = HALFSQRT2 * fabs(qIn[2] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 22;}
	curAbsVal = HALFSQRT2 * fabs(qIn[1] - qIn[3]);
	if( curAbsVal > maxCos) { maxCos = curAbsVal; maxQuatID = 23;}
	char qSign;
	switch (maxQuatID) {
	//q1
	case 0: qSign = SIGN(qIn[0]);
			qOut[0] = maxCos;
			qOut[1] = qSign * (qIn[1]);
			qOut[2] = qSign * (qIn[2]);
			qOut[3] = qSign * (qIn[3]); break;
	//q2
	case 1: qSign = SIGN(-qIn[1]);
			qOut[0] = maxCos;
			qOut[1] = qSign * (qIn[0]);
			qOut[2] = qSign * (qIn[3]);
			qOut[3] = qSign * (-qIn[2]); break;
	//q3
	case 2: qSign = SIGN(-qIn[2]);
			qOut[0] = maxCos;
			qOut[1] = qSign * (-qIn[3]);
			qOut[2] = qSign * (qIn[0]);
			qOut[3] = qSign * (qIn[1]); break;
	//q4
	case 3: qSign = SIGN(-qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * (qIn[2]);
			qOut[2] = qSign * (-qIn[1]);
			qOut[3] = qSign * (qIn[0]); break;
	//q5 +--- +++- +-++ ++-+
	case 4: qSign = SIGN(qIn[0] - qIn[1] - qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (qIn[0] + qIn[1] + qIn[2] - qIn[3]);
			qOut[2] = qSign * 0.5 * (qIn[0] - qIn[1] + qIn[2] + qIn[3]);
			qOut[3] = qSign * 0.5 * (qIn[0] + qIn[1] - qIn[2] + qIn[3]); break;
	//q6 ++++ -+-+ -++- --++
	case 5: qSign = SIGN(qIn[0] + qIn[1] + qIn[2] + qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (-qIn[0] + qIn[1] - qIn[2] + qIn[3]);
			qOut[2] = qSign * 0.5 * (-qIn[0] + qIn[1] + qIn[2] - qIn[3]);
			qOut[3] = qSign * 0.5 * (-qIn[0] - qIn[1] + qIn[2] + qIn[3]); break;
	//q7 +-+- ++++ --++ +--+
	case 6: qSign = SIGN(qIn[0] - qIn[1] + qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (qIn[0] + qIn[1] + qIn[2] + qIn[3]);
			qOut[2] = qSign * 0.5 * (- qIn[0] - qIn[1] + qIn[2] + qIn[3]);
			qOut[3] = qSign * 0.5 * (qIn[0] - qIn[1] - qIn[2] + qIn[3]); break;
	//q8 ++-+ -+-- +++- -+++
	case 7: qSign = SIGN(qIn[0] + qIn[1] - qIn[2] + qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (-qIn[0] + qIn[1] - qIn[2] - qIn[3]);
			qOut[2] = qSign * 0.5 * (qIn[0] + qIn[1] + qIn[2] - qIn[3]);
			qOut[3] = qSign * 0.5 * (-qIn[0] + qIn[1] + qIn[2] + qIn[3]); break;
	//q9 ++-- -++- +-+- ++++
	case 8: qSign = SIGN(qIn[0] + qIn[1] - qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (-qIn[0] + qIn[1] + qIn[2] - qIn[3]);
			qOut[2] = qSign * 0.5 * (qIn[0] - qIn[1] + qIn[2] - qIn[3]);
			qOut[3] = qSign * 0.5 * (qIn[0] + qIn[1] + qIn[2] + qIn[3]); break;
	//q10 +-++ ++-+ -+++ ---+
	case 9: qSign = SIGN(qIn[0] - qIn[1] + qIn[2] + qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (qIn[0] + qIn[1] - qIn[2] + qIn[3]);
			qOut[2] = qSign * 0.5 * (-qIn[0] + qIn[1] + qIn[2] + qIn[3]);
			qOut[3] = qSign * 0.5 * (-qIn[0] - qIn[1] - qIn[2] + qIn[3]); break;
	//q11 +++- -+++ --+- +-++
	case 10: qSign = SIGN(qIn[0] + qIn[1] + qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (-qIn[0] + qIn[1] + qIn[2] + qIn[3]);
			qOut[2] = qSign * 0.5 * (-qIn[0] - qIn[1] + qIn[2] - qIn[3]);
			qOut[3] = qSign * 0.5 * (qIn[0] - qIn[1] + qIn[2] + qIn[3]); break;
	//q12 +--+ ++-- ++++ -+-+
	case 11: qSign = SIGN(qIn[0] - qIn[1] - qIn[2] + qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * 0.5 * (qIn[0] + qIn[1] - qIn[2] - qIn[3]);
			qOut[2] = qSign * 0.5 * (qIn[0] + qIn[1] + qIn[2] + qIn[3]);
			qOut[3] = qSign * 0.5 * (-qIn[0] + qIn[1] - qIn[2] + qIn[3]); break;
	//---------------------------------
	//q13
	case 12: qSign = SIGN(qIn[0] - qIn[1]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[0] + qIn[1]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[2] + qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (-qIn[2] + qIn[3]); break;
	//q14
	case 13: qSign = SIGN(qIn[0] - qIn[2]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[1] - qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[0] + qIn[2]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[1] + qIn[3]); break;
	//q15
	case 14: qSign = SIGN(qIn[0] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[1] + qIn[2]);
			qOut[2] = qSign * HALFSQRT2 * (-qIn[1] + qIn[2]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[0] + qIn[3]); break;
	//q16
	case 15: qSign = SIGN(-qIn[1] - qIn[2]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[0] - qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[0] + qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[1] - qIn[2]); break;
	//q17
	case 16: qSign = SIGN(- qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[2] - qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[0] - qIn[1]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[0] + qIn[1]); break;
	//q18
	case 17: qSign = SIGN(- qIn[1] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[0] + qIn[2]);
			qOut[2] = qSign * HALFSQRT2 * (-qIn[1] + qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[0] - qIn[2]); break;
	//q19
	case 18: qSign = SIGN(qIn[0] + qIn[1]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (-qIn[0] + qIn[1]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[2] - qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[2] + qIn[3]); break;
	//q20
	case 19: qSign = SIGN(qIn[0] + qIn[2]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[1] + qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (-qIn[0] + qIn[2]);
			qOut[3] = qSign * HALFSQRT2 * (-qIn[1] + qIn[3]); break;
	//q21
	case 20: qSign = SIGN(qIn[0] + qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[1] - qIn[2]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[1] + qIn[2]);
			qOut[3] = qSign * HALFSQRT2 * (-qIn[0] + qIn[3]); break;
	//q22
	case 21: qSign = SIGN(qIn[1] - qIn[2]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (-qIn[0] - qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (qIn[0] - qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[1] + qIn[2]); break;
	//q23
	case 22: qSign = SIGN(qIn[2] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (qIn[2] + qIn[3]);
			qOut[2] = qSign * HALFSQRT2 * (-qIn[0] - qIn[1]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[0] - qIn[1]); break;
	//q24
	case 23: qSign = SIGN(qIn[1] - qIn[3]);
			qOut[0] = maxCos;
			qOut[1] = qSign * HALFSQRT2 * (-qIn[0] + qIn[2]);
			qOut[2] = qSign * HALFSQRT2 * (-qIn[1] - qIn[3]);
			qOut[3] = qSign * HALFSQRT2 * (qIn[0] + qIn[2]); break;
	}

}

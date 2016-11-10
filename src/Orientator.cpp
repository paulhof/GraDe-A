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

#include "Orientator.h"
#define COS60DEG .5
#define COS120DEG -.5
#define COS90DEG .0
#define COS35DEG SQRT2_3 //35.26 Degree, angle between 110 and 111
#define COS45DEG HALFSQRT2
#define COS135DEG -HALFSQRT2

#define COS1DEGPREC (1.5230484360876084299E-4)
#define COS2DEGPREC (6.0917298090426999376E-4)
#define COS3DEGPREC (1.3704652454261262155E-3)
#define COS4DEGPREC (2.4359497401757523868E-3)
#define COS5DEGPREC (3.8053019082544677050E-3)
#define COS10DEGPREC (1.52E-2)
//
#define SINSMALLPREC (1.7452406437283512819E-3)
#define SIN1DEGPREC (1.7452406437283512819E-2)
#define SIN2DEGPREC (3.4899496702500971646E-2)
#define SIN3DEGPREC (5.2335956242943832722E-2)
#define SIN4DEGPREC (6.9756473744125300776E-2)
#define SIN5DEGPREC (8.7155742747658173558E-2)
//
#define ST_1DEGPREC (0.0151904)
#define ST_2DEGPREC (0.0305284)
#define ST_3DEGPREC (0.0460095)
#define ST_4DEGPREC (0.0616289)
#define ST_5DEGPREC (0.0773817)
//
#define F4TY5_5DEGPREC (0.06432)
#define TH3TY5_5DEGPREC (0.05343)
#define COSTRESHOLD COS5DEGPREC
#define SINTRESHOLD 2.e-1
#define SIXTYDEGTRESHH ST_1DEGPREC
#define DEFAULTLEAFSIZE 10

Orientator::Orientator(AtomBox * inBoxes) {
	init(inBoxes, ORIENTALLOC);
}

Orientator::Orientator(AtomBox * inBoxes, unsigned long initCapacity){
	init(inBoxes, initCapacity);
}

void Orientator::init(AtomBox * inBoxes, unsigned long initCapacity){
	nOrientations = 0;
	orientAlloc = initCapacity;
	orientSize = 0;
	orientations = nullptr;
	boxes = inBoxes;
}

Orientator::~Orientator() {
	if(orientSize > 0){
		delete [] orientations;
	}
}

oID Orientator::orientateFCC(double * neighborPositions, unsigned char nNextNeighbors){
	if (nNextNeighbors > 12) {
		return NO_ORIENTATION;
	}
	double * unitVectors = new double [nNextNeighbors * DIM];
	for(unsigned char i = 0; i < nNextNeighbors * DIM; i ++){
		unitVectors[i] = neighborPositions[i];
	}
	ori::unitizeVectors(unitVectors, nNextNeighbors);
	unsigned char nVects;
	nVects = reduceAntiparallelVectors(unitVectors, nNextNeighbors);
	double * nextNborVects = new double [nVects*DIM];
	for(long i = 0; i < nVects*DIM; i++){
		nextNborVects[i] = unitVectors[i];
	}
	delete [] unitVectors;
	//nextNborVects has now size nVects
#ifndef DEBUGMODE
	if (nVects < 6) {
		delete [] nextNborVects;
		return NO_ORIENTATION;
	}
#else
	if (nVects < 6) {
	delete [] nextNborVects;
	std::cout << "NO ORI BECAUSE OF nVects < 6" << std::endl;
	return NO_ORIENTATION;
	}
#endif
	//each vector pair:
	char iVec2;
	char nPerpend = 0;
	double normal100Vects[3*DIM];
	//obtain <100> directions
	//always 2 next-neighbor pairs which are perpendicular lay inside a {100} plane
	//the normal vector is obtained by calculating the cross product
	for (char iVec = 0; iVec < nVects; iVec++){
		if (3 == nPerpend) break;
		for (iVec2 = iVec + 1; iVec2 < nVects; iVec2++){
			if( vectsArePerpend(nextNborVects + DIM * iVec, nextNborVects + DIM * iVec2) )
			{
				ori::crossProduct(nextNborVects + DIM * iVec, nextNborVects + DIM * iVec2, normal100Vects + DIM * nPerpend);
				if(3 == ++nPerpend){
					break;
				}
			}
		}
	}
	delete [] nextNborVects;
	if (nPerpend < 3) {
		//not enough perpend directions found
		return NO_ORIENTATION;
	}
	return calcOrientationFromThree100Directions(normal100Vects, normal100Vects + DIM, normal100Vects + 2*DIM);
}

oID Orientator::calcOrientationFromThree100Directions(double * v100, double * v010, double * v001){
#ifdef USE_ARMADILLO
	arma::mat m(DIM,DIM);
#else
	Eigen::Matrix3d m;
#endif
	double det;
	//row 0
	m(0,0) = v100[0]; m(0,1) = v100[1]; m(0,2) = v100[2];
	//row 1
	m(1,0) = v010[0]; m(1,1) = v010[1]; m(1,2) = v010[2];
	//row 2
	m(2,0) = v001[0]; m(2,1) = v001[1];	m(2,2) = v001[2];

#ifdef USE_ARMADILLO
	det = arma::det(m);
#else
	det = m.determinant();
#endif
	//check if determinant is negative (matrix is left-handed)
	if (det < 0.){
		m(0,0) = -v100[0];
		m(0,1) = -v100[1];
		m(0,2) = -v100[2];
	}
#ifdef USE_ARMADILLO
	det = arma::det(m);
#else
	det = m.determinant();
#endif
	if (det < 0.){
		std::cout << "ERROR: determinant of matrix negative" << std::endl;
		return NO_ORIENTATION;
	}
	return closestOrientation(m);
}

#ifdef USE_ARMADILLO
oID Orientator::closestOrientation (arma::mat  & m3x3){
#else
oID Orientator::closestOrientation (Eigen::Matrix3d  & m3x3){
#endif
	double M[9] = {
	m3x3(0,0), m3x3(0,1), m3x3(0,2),
	m3x3(1,0), m3x3(1,1), m3x3(1,2),
	m3x3(2,0), m3x3(2,1), m3x3(2,2)
	};
	return closestOrientation(M);
}

oID Orientator::calcFCCOrientation_90Deg(double * v110, double * vm110){
	double v001[3];
	ori::crossProduct(v110,vm110,v001);
	double M [9] = {
	//row 0
	HALFSQRT2 * (v110[0] - vm110[0]) ,//11-0
	HALFSQRT2 * (v110[0] + vm110[0]) ,//12-1
	v001[0] ,		//13-2
	//row 1
	HALFSQRT2 * (v110[1] - vm110[1]) ,//21-3
	HALFSQRT2 * (v110[1] + vm110[1]) ,//22-4
	v001[1],		//23-5
	//row 2
	HALFSQRT2 * (v110[2] - vm110[2]) ,//31-6
	HALFSQRT2 * (v110[2] + vm110[2]) ,//32-7
	v001[2]//33-8
	};
	return closestOrientation(M);
}

oID Orientator::calcFCCOrientation_60Deg(double * v110, double * v101){
	double v1m1m1[3];
	ori::crossProduct(v110,v101,v1m1m1);
	for(char i = 0; i < DIM; i++){
		v110[i] *= THIRDSQRT2;
		v101[i] *= THIRDSQRT2;
		v1m1m1[i] *= THIRD2;
	}
	double M[9] = {
	//row 0
	v110[0] + v101[0] + v1m1m1[0],
	2 * v110[0] - v101[0] - v1m1m1[0],
	- v110[0] + 2* v101[0] - v1m1m1[0],
	//row 1
	v110[1] + v101[1] + v1m1m1[1],
	2 * v110[1] - v101[1] - v1m1m1[1],
	- v110[1] + 2* v101[1] - v1m1m1[1],
	//row 2
	v110[2] + v101[2] + v1m1m1[2],
	2 * v110[2] - v101[2] - v1m1m1[2],
	- v110[2] + 2* v101[2] - v1m1m1[2]
	};
	return closestOrientation(M);
}

void Orientator::matrixToClosestQuaternion(const double * M, double * q){
	//see Itzhack 2000 " New Method for Extracting the Quaternion from a Rotation Matrix "
#ifdef USE_ARMADILLO
	arma::mat m(4,4);
	arma::vec mEigval;
	arma::mat mEigvec;
#else
	Eigen::Matrix4d m;
#endif
	m(0,0) = ONETHIRD * (M[0] - M[4] - M[8]);//d11 - d22 - d33 : 0
	m(0,1) = ONETHIRD * (M[3] + M[1]);		//d21 + d12 : 1
	m(0,2) = ONETHIRD * (M[6] + M[2]);		//d31 + d13 : 2
	m(0,3) = ONETHIRD * (M[5] - M[7]);		//d23 - d32 : 3
	//row 1
	m(1,0) = m(0,1);//d21 + d12 : 4
	m(1,1) = ONETHIRD * (M[4] - M[0] - M[8]);//d22 - d11 - d33 : 5
	m(1,2) = ONETHIRD * (M[7] + M[5]);		//d32 + d23 : 6
	m(1,3) = ONETHIRD * (M[6] - M[2]); 		//d31 - d13 : 7
	//row 2
	m(2,0) = m(0,2);//d31 + d13 : 8
	m(2,1) = m(1,2);//d32 + d23 : 9
	m(2,2) = ONETHIRD * (M[8] - M[0] - M[4]);//d33 - d11 - d22 : 10
	m(2,3) = ONETHIRD * (M[1] - M[3]);//d12 - d21 : 11
	//row 3
	m(3,0) = m(0,3);//d23 - d32 : 12
	m(3,1) = m(1,3);//d31 - d13 : 13
	m(3,2) = m(2,3);//d12 - d21 : 14
	m(3,3) = ONETHIRD * (M[0] + M[4] + M[8]);//d11 + d22 + d33 : 15

	char nEigvals;
#ifndef USE_ARMADILLO
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigenSolver;
#endif
	try{
#ifdef USE_ARMADILLO
		arma::eig_sym(mEigval, mEigvec, m);
		nEigvals = mEigval.n_elem;
#else//Use Eigen solver
		eigenSolver.compute(m);
		nEigvals = 4;

#endif
	} catch(...){
		return;
	}

	if( nEigvals > 0){
#ifndef USE_ARMADILLO
	const Eigen::Matrix4d & mEigvec = eigenSolver.eigenvectors();
#endif
	q[0] = mEigvec(3,nEigvals - 1);
	q[1] = mEigvec(0,nEigvals - 1);
	q[2] = mEigvec(1,nEigvals - 1);
	q[3] = mEigvec(2,nEigvals - 1);
	}
}

oID Orientator::closeOrientbyQuaternion(double * q){
	double cubicQuat[4];
	ori::uniqueCubicRotationQuaternion(q,cubicQuat);
	return newOrientbyQuaternion(cubicQuat);
}

double Orientator::cubicCosHalfMisOrientation(const Atom* atom1, const Atom* atom2) const {
	const double  * q1, * q2;
		q1= orientations[atom1->getOrientationId()].getQuaternion();
		q2= orientations[atom2->getOrientationId()].getQuaternion();
		return ori::cubicCosHalfMisOrientation(q1,q2);
}

bool Orientator::haveCloseOrientations(const Atom *atom1, const Atom *atom2, double cosHalfThreshold) const{
	//return cubicCosHalfMisOrientation(atom1, atom2) > cosHalfThreshold;
	//cosHalfMis seems to be sufficient for small angles
	return cosHalfMisOrientation(atom1, atom2) > cosHalfThreshold;
}

double Orientator::cosHalfMisOrientation(const Atom *atom1, const Atom *atom2) const{
	const double  * q1, * q2;
	q1= orientations[atom1->getOrientationId()].getQuaternion();
	q2= orientations[atom2->getOrientationId()].getQuaternion();
	return ori::cosHalfMisOrientation(q1,q2);
}

const Orientation * Orientator::getOrientation(oID inOriId) const
{
	return orientations + inOriId;
}

oID Orientator::newOrientbyQuaternion(double * q){
	if(nOrientations < orientSize){//space is sufficient
		orientations[nOrientations].initbyQuaternion(q);
		nOrientations ++;
		return nOrientations - 1;
	}
	Orientation * oldOrients;
	if (0 != nOrientations){//reallocation
		 oldOrients = new Orientation [nOrientations];
		for(long i = 0; i < nOrientations; i ++){
				oldOrients[i] = orientations[i];
		}
		delete [] orientations;

	}
	orientSize += orientAlloc;
    orientations = new Orientation[orientSize];

    if(0 != nOrientations){
		for (long i = 0; i < nOrientations; i ++){
			orientations[i] = oldOrients[i];
		}
		delete [] oldOrients;
    }

	orientations[nOrientations].initbyQuaternion(q);
	nOrientations ++;
	return nOrientations - 1;
}

AtomBox * Orientator::getBoxes(){
	return boxes;
}

oID Orientator::closestOrientation (const double * M){
	double q[4];
	try{
		matrixToClosestQuaternion(M, q);
	} catch (...){
		return NO_ORIENTATION;
	}
	oID ret;
	 ret = closeOrientbyQuaternion(q);
	return ret;
}

unsigned char Orientator::reduceAntiparallelVectors(doubleP &vects, unsigned char n){
	double * v1, * v2;
	bool * redundant = new bool [n];
	for(unsigned char i = 0; i < n; i++ ){
		redundant[i] = false;
	}
	unsigned char ii;
	std::vector<double> v;

	for(unsigned char i = 0; i < n; i++ ){
			v1 = vects + i*DIM;
			//check whether there is any vector to v1 with a near 180Deg relationship
			for(ii = i+1; ii < n; ii++){
				//both pair atoms should be unmarried
					if(!redundant[i] && !redundant[ii]){
						v2 = vects + ii*DIM;
						//check the angle (is close to 180Deg? (scalar-product close to -1?))
						if (fabs(fabs(ori::scalarProduct(v1, v2)) - 1 )< COSTRESHOLD ){
							//now both are forced to wear wedding rings
							redundant[i] = true;
							redundant[ii] = true;
							//save the mean direction as vector
							v.push_back(.5*(v1[0]-v2[0]));
							v.push_back(.5*(v1[1]-v2[1]));
							v.push_back(.5*(v1[2]-v2[2]));
						}
					}
			}
			//if no partner was found, remain single forever!
			if(!redundant[i]){
				v.push_back(v1[0]);
				v.push_back(v1[1]);
				v.push_back(v1[2]);
			}
	}
	//redundant information is no longer needed
	delete [] redundant;
	//initial vectors are no longer needed
	delete [] vects;
	//allocate for reduced vector list
	vects = new double [v.size()];
	for (unsigned char i = 0; i < v.size(); i ++){
		vects[i] = v[i];
	}
	return static_cast<unsigned char> (v.size()/DIM);
}

bool Orientator::vectsArePerpend(double * v1, double * v2){
	if ( fabs(ori::scalarProduct(v1, v2)) < SINTRESHOLD) return true;
	return false;
}

bool Orientator::findBestPerpendPair(double * directs, unsigned char nDirects, unsigned char &v110Id, unsigned char &vm110Id){
	double * vTest1;
	double * vTest2;
	double minDev = 1.;
	double curDev = 1.;
	for (unsigned char i = 0; i < nDirects; i ++){
		vTest1 = directs + i*DIM;
		for ( unsigned char ii = i+1; ii < nDirects; ii++ ){
			vTest2 = directs + ii*DIM;
			curDev = fabs(ori::scalarProduct(vTest1, vTest2));
			if ( curDev < minDev ){
				v110Id = i;
				vm110Id = ii;
				minDev = curDev;
			}
		}
	}
	if (minDev < SINTRESHOLD) return true;
	return false;
}

unsigned char Orientator::findPerpendDirect(double * directs, unsigned char nDirects, unsigned char me){
	double * vTest = directs;
	double * vMe = &vTest[me];
	for (unsigned char i = 0; i < nDirects; i ++){
		if( i != me){
			if(fabs(ori::scalarProduct(vTest, vMe)) < SINTRESHOLD){
				return i;
			}
		}
		vTest += DIM;
	}
	return me;
}

unsigned char Orientator::find60DegDirect(double * directs, unsigned char nDirects, unsigned char me){
	double * vTest = directs;
	double * vMe = &vTest[me];
	for (unsigned char i = 0; i < nDirects; i ++){
		if( i != me){
			if(fabs(ori::scalarProduct(vTest, vMe) - .5) < SIXTYDEGTRESHH)
				//60Degree found
				return i;
			if(fabs(ori::scalarProduct(vTest, vMe) + .5) < SIXTYDEGTRESHH){
				//120Degree found
				ori::inverseDirection(vTest);
				return i;
			}
		}
		vTest += DIM;
	}
	return me;
}

bool sortVectsInSpace(double * vA, double  * vB){ return (7*vA[0] + 41*vA[1] + vA[2]) > (7*vB[0] + 41*vB[1] + vB[2]); }

void Orientator::sortVectsList(double * vList, long nVects){
	double ** vects = new double * [nVects];
	double * sortedValues = new double [nVects * DIM];
	for (long i = 0; i < nVects; i++){
		vects[i] = vList + i*DIM;
	}
	std::sort(vects,vects + nVects,sortVectsInSpace);
	long iP;
	for (long i = 0; i < nVects; i++){
	iP = i*DIM;
	sortedValues[iP] = vects[i][0];
	sortedValues[iP+1] = vects[i][1];
	sortedValues[iP+2] = vects[i][2];
	}
	delete [] vects;
	for (long i = 0; i < nVects; i++){
	iP = i*DIM;
	vList[iP] = sortedValues[iP];
	iP++;
	vList[iP] = sortedValues[iP];
	iP++;
	vList[iP] = sortedValues[iP];
	}
	delete [] sortedValues;
}

long Orientator::getNumOrientations() const{
	return nOrientations;
}

Orientation * Orientator::getOrientations(){
	return orientations;
}

OrientatorPrinter::OrientatorPrinter(const Orientator* inOrient,
		const std::string inFileName) {
	orient = inOrient;
	writer.setFileName(inFileName);
	writer.addNewColumn("phi_1 [degree]");
	writer.addNewColumn("PHI [degree]");
	writer.addNewColumn("phi_2 [degree]");
	writer.addNewColumn("q0");
	writer.addNewColumn("q1");
	writer.addNewColumn("q2");
	writer.addNewColumn("q3");
}

void OrientatorPrinter::print() {
	for (oID iO = 0; iO < orient->getNumOrientations(); iO++){
		printOrientation(orient->getOrientation(iO));
	}
	writer.write();
}

void OrientatorPrinter::printOrientation(const Orientation* ori) {
	writer.addNewLine();
	const double * q = ori->getQuaternion();
	double bunge[3];
	ori->obtainBungeAngles(bunge);
	for (char i = 0; i < 3; i++){
		writer.addEntryToCurrentLine(ori::to_string(bunge[i]*RADTODEG),i);
	}
	for (char i = 0; i < 4; i++){
		writer.addEntryToCurrentLine(ori::to_string(q[i]),i+3);
	}
}

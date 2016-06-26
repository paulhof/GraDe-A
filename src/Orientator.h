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

#ifndef ORIENTATOR_H_
#define ORIENTATOR_H_

#define ORIENTALLOC 10000
#include "GradeA_Defs.h"
#include "Atom.h"
#include "Orientation.h"
#include <armadillo>
#include "io/CSVTableWriter.h"
typedef double * doubleP;
typedef struct{
	long iB;
	long iA;
}AtomID;
class Atom;
class AtomBox;
class MeanOrientation;

class Orientator {
public:
	Orientator(){};
	Orientator(AtomBox * boxes);
	Orientator(AtomBox * boxes, unsigned long initCapacity);
	virtual ~Orientator();
	oID orientateFCC(double * neighborPositions, unsigned char nNextNeighbors);
	oID closestOrientation (arma::mat  & m3x3);
	oID closestOrientation (const double * M);
	long getNumOrientations() const;
	void matrixToClosestQuaternion (const double * M, double * q);
	void sortVectsList(double * vList, long nVects);
	Orientation * getOrientations();
	const Orientation * getOrientation(oID inOriId) const;
	AtomBox * getBoxes();
	double cosHalfMisOrientation(const Atom * atom1, const Atom * atom2) const;
	double cubicCosHalfMisOrientation(const Atom *atom1, const Atom *atom2) const;
	bool haveCloseOrientations(const Atom * atom1, const Atom * atom2, double cosHalfThreshold) const;

private:
	void init(AtomBox * boxes, unsigned long initCapacity);
	oID calcOrientationFromThree100Directions(double * v100, double * v010, double * v001);
	oID calcFCCOrientation_90Deg(double * v110, double * vm110);
	oID calcFCCOrientation_60Deg(double * v110, double * v101);
	//void inverseDirection(double * direct);
	unsigned char reduceAntiparallelVectors(doubleP &vects, unsigned char n);
	unsigned char find60DegDirect(double * directs, unsigned char nDirects, unsigned char me);
	unsigned char findPerpendDirect(double * directs, unsigned char nDirects, unsigned char me);
	bool vectsArePerpend(double * v1, double * v2);
	bool findBestPerpendPair(double * directs, unsigned char nDirects, unsigned char &v110Id, unsigned char &vm110Id);
	oID closeOrientbyQuaternion(double * q);
	oID newOrientbyQuaternion(double * q);
	AtomBox * boxes = nullptr;
	Orientation * orientations = nullptr;

	unsigned long orientAlloc = ORIENTALLOC;
	long nOrientations = 0;
	long orientSize = 0;
};

class OrientatorPrinter {
public:
	OrientatorPrinter(const Orientator * inOrient, const std::string inFileName);
	void print();
private:
	void printOrientation(const Orientation * ori);
	const Orientator * orient;
	CSVTableWriter writer;
};

#endif /* ORIENTATOR_H_ */

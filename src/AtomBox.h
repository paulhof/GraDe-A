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

#ifndef ATOMBOX_H_
#define ATOMBOX_H_
#include "Atom.h"
#include "GradeA_Defs.h"
#include "Orientator.h"
struct ABoxNeighbor;
class AtomBox {
public:
	AtomBox();
	virtual ~AtomBox();
	void init(double * origin, double * size, ABoxNeighbor * neighbors, long nNeighbors, bool sizeLinked = true);
	void srtNeighbors();
	void addAtom(double * pos);
	void calculateAtomOrientationsFCC(long boxId, double angleThreshold, Orientator  * orient, const double rSqrMin, const double rSqrMax);
	unsigned char atomNeighbors(const long atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, double * nborPositions);
	unsigned char atomNeighbors(const long atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, AtomBox ** outNborBoxesList, long * outNborAtomIdList, double * outNborPosList);
	unsigned char nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, double * nborPositions);
	unsigned char nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, Atom ** nborAtoms);
	Atom * getAtom(long atomId);
	void obtainGlobalAtomPos(long atomId, double * outPos) const;
	void printAtoms();
	double * getOrigin();
	double * getSize();
	long getNumNeighbors();
	long getNumAtoms();
	ABoxNeighbor * getNeighbors();
	long mostFrequentGrainId() const;
private:
	inline double sqrDist(const double * p1, const double * p2) const;
	bool sizeLinked;
	double origin[DIM];
	double * size;
	Atom * atoms = nullptr;
	long nAtoms;
	long atomArrSize;
	ABoxNeighbor * neighbors = nullptr;
	long nNeighbors;
};
typedef AtomBox* AtomBoxP;


struct ABoxNeighbor{
	AtomBoxP box;
	char coord[DIM];
};

#endif /* ATOMBOX_H_ */

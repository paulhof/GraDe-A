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

#include "AtomBox.h"
AtomBox::AtomBox() {
	nAtoms = 0;
	atomArrSize = 0;
	nNeighbors = 0;
	atoms = nullptr;
	neighbors = nullptr;
	size = nullptr;
	sizeLinked = false;
}

void AtomBox::init(const double * inOrigin, const double * inSize, const ABoxNeighbor * inNeighbors, long inNumNeighbors, bool inSizeLinked){
	sizeLinked = inSizeLinked;
	if (sizeLinked){
		size = (double*) inSize;
	} else {
		size = new double [DIM];
		for (long i = 0; i < DIM; i++){
		size[i] = inSize[i];
		}
	}

	for (long i = 0; i < DIM; i++){
		origin[i] = inOrigin[i];
	}

	nNeighbors = inNumNeighbors;
	neighbors = new ABoxNeighbor[nNeighbors];
	for(long iNbor = 0; iNbor < nNeighbors; iNbor ++){
		neighbors[iNbor] = inNeighbors[iNbor];
	}
}

AtomBox::~AtomBox() {
	if (!sizeLinked) delete [] size;
	delete [] neighbors;
	delete [] atoms;
}

void AtomBox::addAtom(double * inPos ){
	if(0 == nAtoms){
		atomArrSize = ATOMALLOC;
		atoms = new Atom[atomArrSize];
	} else  if ( nAtoms == atomArrSize){
		long i;
		Atom * oldAtoms = new Atom[nAtoms];

		for (i = 0; i < nAtoms; i ++){
			oldAtoms[i] = atoms[i];
		}

		delete [] atoms;
		atomArrSize += ATOMALLOC;
		atoms = new Atom[atomArrSize];

		for (i = 0; i < nAtoms; i ++){
			atoms[i] = oldAtoms[i];
		}

		delete [] oldAtoms;
	}
	atoms[nAtoms].init(inPos, NO_ORIENTATION, NO_GRAIN);
	nAtoms ++;
}

void AtomBox::obtainGlobalAtomPos(long atomId, double * outPos) const{
	double * relPos = atoms[atomId].getPos();
	outPos[0] = origin[0] + relPos[0];
	outPos[1] = origin[1] + relPos[1];
	outPos[2] = origin[2] + relPos[2];
}

void AtomBox::calculateAtomOrientationsFCC(double angleThreshold, Orientator * orient, const double rSqrMin, const double rSqrMax){
	unsigned char nAtomNeighbors;
	double atomNborPositions[12*DIM];
	Atom * atom;
	for (long iA = 0; iA < nAtoms; iA++){
		atom = atoms + iA;
#ifndef NEAREST_ATOMNEIGHBORHOOD
		nAtomNeighbors = atomNeighbors(iA, rSqrMin , rSqrMax , 12, atomNborPositions);
#else
		nAtomNeighbors = nearestAtomNeighbors(iA, 12, atomNborPositions);
#endif
		long n = orient->orientateFCC(atomNborPositions, nAtomNeighbors);
		atom->setOrientationId(n);
	}
}

unsigned char AtomBox::atomNeighbors(const long atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, double * outNborPositions){
	double * atomPos = atoms[atomId].getPos();
	double * nborAtomPos = nullptr;
	unsigned char nFoundNeighbors = 0;
	double sqrDistance; //distance atom-atom squared
	int pos;
	//check own AtomBox first
	long iA;
	for ( iA = 0; iA < nAtoms; iA++){
		if( iA != atomId){
			nborAtomPos = atoms[iA].getPos();
			sqrDistance = sqrDist(atomPos, nborAtomPos);
			if( sqrDistance < rSqrMax && sqrDistance > rSqrMin){
				if(nFoundNeighbors == nMaxAtomNeighbors) return 0;
				pos = nFoundNeighbors * DIM;
				outNborPositions[pos] = nborAtomPos[0] - atomPos[0];
				outNborPositions[pos+1] = nborAtomPos[1] - atomPos[1];
				outNborPositions[pos+2] = nborAtomPos[2] - atomPos[2];
				nFoundNeighbors ++;
			}
		}
	}
	//check neighboring AtomBoxes next
	if(nFoundNeighbors < nMaxAtomNeighbors){
	AtomBoxP nborBox = nullptr;
	Atom * nborAtom;
	double relAtomPos [DIM];//relative Position of the atom in neighbor's box perspective
	for (long iBoxes = 0; iBoxes < nNeighbors; iBoxes ++ ){
		nborBox = neighbors[iBoxes].box;
		relAtomPos[0] = atomPos[0] - neighbors[iBoxes].coord[0] * size[0];
		relAtomPos[1] = atomPos[1] - neighbors[iBoxes].coord[1] * size[1];
		relAtomPos[2] = atomPos[2] - neighbors[iBoxes].coord[2] * size[2];
		for(iA = 0; iA < nborBox->getNumAtoms(); iA++){
			nborAtom = nborBox->getAtom(iA);
			nborAtomPos = nborAtom->getPos();
			sqrDistance = sqrDist(relAtomPos, nborAtomPos);
			if(sqrDistance < rSqrMax && sqrDistance > rSqrMin){
				if(nFoundNeighbors == nMaxAtomNeighbors) return 0;
				pos = nFoundNeighbors * DIM;
				outNborPositions[pos] = nborAtomPos[0] - relAtomPos[0];
				outNborPositions[pos+1] = nborAtomPos[1] - relAtomPos[1];
				outNborPositions[pos+2] = nborAtomPos[2] - relAtomPos[2];
				nFoundNeighbors ++;
			}
		}
	}
	}
	return nFoundNeighbors;
}

unsigned char AtomBox::atomNeighbors(const long  atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, AtomBoxP *outNborBoxesList, long * outNborAtomIdList, double * outNborPosList)
{
	double * atomPos = atoms[atomId].getPos();//relative position of the atom in the box
	double * nborAtomPos = nullptr;
	unsigned char nFoundNeighbors = 0;
	double sqrDistance; //distance atom-atom squared
	int pos;
		//check own AtomBox first
		long iA;
		for ( iA = 0; iA < nAtoms; iA++){
			if( iA != atomId){
				nborAtomPos = atoms[iA].getPos();//relative position of the neighbor atom in the box
				sqrDistance = sqrDist(atomPos, nborAtomPos);
				if( sqrDistance < rSqrMax && sqrDistance > rSqrMin){
					if(nFoundNeighbors == nMaxAtomNeighbors) return 0;
					outNborBoxesList[nFoundNeighbors] = this;
					outNborAtomIdList[nFoundNeighbors] = iA;
					//
					pos = nFoundNeighbors * DIM;
					outNborPosList[pos] = nborAtomPos[0] - atomPos[0];
					outNborPosList[pos+1] = nborAtomPos[1] - atomPos[1];
					outNborPosList[pos+2] = nborAtomPos[2] - atomPos[2];
					nFoundNeighbors ++;
				}
			}
		}
		//check neighboring AtomBoxes next
		if(nFoundNeighbors < nMaxAtomNeighbors){
		AtomBoxP nborBox = nullptr;
		Atom * nborAtom;
		double relAtomPos [DIM];//relative Position of the atom in neighbor's box perspective
		for (long iBoxes = 0; iBoxes < nNeighbors; iBoxes ++ ){
			nborBox = neighbors[iBoxes].box;
			relAtomPos[0] = atomPos[0] - neighbors[iBoxes].coord[0] * size[0];
			relAtomPos[1] = atomPos[1] - neighbors[iBoxes].coord[1] * size[1];
			relAtomPos[2] = atomPos[2] - neighbors[iBoxes].coord[2] * size[2];
			for(iA = 0; iA < nborBox->getNumAtoms(); iA++){
				nborAtom = nborBox->getAtom(iA);
				nborAtomPos = nborAtom->getPos();
				sqrDistance = sqrDist(relAtomPos, nborAtomPos);
				if(sqrDistance < rSqrMax && sqrDistance > rSqrMin){
					if(nFoundNeighbors == nMaxAtomNeighbors) return 0;
					outNborBoxesList[nFoundNeighbors] = nborBox;
					outNborAtomIdList[nFoundNeighbors] = iA;
					//
					pos = nFoundNeighbors * DIM;
					outNborPosList[pos] = nborAtomPos[0] - relAtomPos[0];
					outNborPosList[pos+1] = nborAtomPos[1] - relAtomPos[1];
					outNborPosList[pos+2] = nborAtomPos[2] - relAtomPos[2];
					nFoundNeighbors++;
				}
			}
		}
		}
	return nFoundNeighbors;
}

bool sortLengthAscending(double * vA, double * vB) { return SQR(vA[0])+SQR(vA[1])+SQR(vA[2]) < SQR(vB[0])+SQR(vB[1])+SQR(vB[2]); }

unsigned char AtomBox::nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, double * outNborPositions){
	//finds the nAtomNeighbors next neighbors of a given atom
	double * atomPos = atoms[atomId].getPos();
	double * nborAtomPos = nullptr;
	//list to sort afterwards
	std::vector<double> nborAtomPosList;

	double sqrDistance; //distance atom-atom squared
	int pos;
	//check own AtomBox first
	long iA;
	for ( iA = 0; iA < nAtoms; iA++){
		if( iA != atomId){
			nborAtomPos = atoms[iA].getPos();
			nborAtomPosList.push_back(nborAtomPos[0] - atomPos[0]);
			nborAtomPosList.push_back(nborAtomPos[1] - atomPos[1]);
			nborAtomPosList.push_back(nborAtomPos[2] - atomPos[2]);
		}
	}
	AtomBoxP nborBox = nullptr;
	Atom * nborAtom;
	double relAtomPos [DIM];//relative position of the atom in neighbor's box perspective
	for (long iBoxes = 0; iBoxes < nNeighbors; iBoxes ++ ){
		nborBox = neighbors[iBoxes].box;
		relAtomPos[0] = atomPos[0] - neighbors[iBoxes].coord[0] * size[0];
		relAtomPos[1] = atomPos[1] - neighbors[iBoxes].coord[1] * size[1];
		relAtomPos[2] = atomPos[2] - neighbors[iBoxes].coord[2] * size[2];
		for(iA = 0; iA < nborBox->getNumAtoms(); iA++){
			nborAtom = nborBox->getAtom(iA);
			nborAtomPos = nborAtom->getPos();
			nborAtomPosList.push_back( nborAtomPos[0] - relAtomPos[0]);
			nborAtomPosList.push_back( nborAtomPos[1] - relAtomPos[1]);
			nborAtomPosList.push_back( nborAtomPos[2] - relAtomPos[2]);
		}
	}
	long nFoundNeighbors = nborAtomPosList.size()/DIM;
	double ** nborAtomList = new double * [nFoundNeighbors];
	double * listBegin = &nborAtomPosList.front();
	//save pointer to vectors in list
	for(long iN = 0; iN < nFoundNeighbors; iN++){
		nborAtomList[iN] = listBegin + iN * DIM;
	}
	//sort neighbor list (neighbors by distance)
	std::sort(nborAtomList, nborAtomList+nFoundNeighbors, sortLengthAscending);

	if (nFoundNeighbors > nAtomNeighbors) nFoundNeighbors = nAtomNeighbors;
	//copy position data in outputarray
	for(unsigned char iN = 0; iN < nFoundNeighbors; iN++){
		pos = iN*DIM;
		outNborPositions[pos++] = nborAtomList[iN][0];
		outNborPositions[pos++] = nborAtomList[iN][1];
		outNborPositions[pos] = nborAtomList[iN][2];
	}
	delete [] nborAtomList;
	return nFoundNeighbors;
}

typedef struct{
	double position[DIM];
	Atom * atom;
} AtomNeighbor;
bool sortLengthAscending2(AtomNeighbor nA, AtomNeighbor nB) { return SQR(nA.position[0])+SQR(nA.position[1])+SQR(nA.position[2]) <
		SQR(nB.position[0])+SQR(nB.position[1])+SQR(nB.position[2]);
}

unsigned char AtomBox::nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, Atom ** outAtoms){
	//finds the nAtomNeighbors next neighbors of a given atom
	double * atomPos = atoms[atomId].getPos();
	double * nborAtomPos = nullptr;
	//list to sort afterwards
	std::vector<AtomNeighbor> nborAtomList;
	//check own AtomBox first
	long iA;
	AtomNeighbor curNbor;
	for ( iA = 0; iA < nAtoms; iA++){
		if( iA != atomId){
			nborAtomPos = atoms[iA].getPos();
			curNbor.atom = atoms +iA;
			curNbor.position[0] = nborAtomPos[0] - atomPos[0];
			curNbor.position[1] = nborAtomPos[1] - atomPos[1];
			curNbor.position[2] = nborAtomPos[2] - atomPos[2];
			nborAtomList.push_back(curNbor);
		}
	}
	AtomBoxP nborBox = nullptr;
	Atom * nborAtom;
	double relAtomPos [DIM];//relative position of the atom in neighbor's box perspective
	for (long iBoxes = 0; iBoxes < nNeighbors; iBoxes ++ ){
		nborBox = neighbors[iBoxes].box;
		relAtomPos[0] = atomPos[0] - neighbors[iBoxes].coord[0] * size[0];
		relAtomPos[1] = atomPos[1] - neighbors[iBoxes].coord[1] * size[1];
		relAtomPos[2] = atomPos[2] - neighbors[iBoxes].coord[2] * size[2];
		for(iA = 0; iA < nborBox->getNumAtoms(); iA++){
			nborAtom = nborBox->getAtom(iA);
			curNbor.atom = nborAtom;
			nborAtomPos = nborAtom->getPos();
			curNbor.position[0] = nborAtomPos[0] - relAtomPos[0];
			curNbor.position[1] = nborAtomPos[1] - relAtomPos[1];
			curNbor.position[2] = nborAtomPos[2] - relAtomPos[2];
			nborAtomList.push_back(curNbor);
		}
	}
	long nFoundNeighbors = nborAtomList.size();
	//sort neighbor list (neighbors by distance)
	std::sort(nborAtomList.begin(), nborAtomList.end(), sortLengthAscending2);
	if (nFoundNeighbors > nAtomNeighbors) nFoundNeighbors = nAtomNeighbors;
	for(unsigned char iN = 0; iN < nFoundNeighbors; iN++){
		outAtoms[iN] = nborAtomList[iN].atom;
	}
	return nFoundNeighbors;
}

void AtomBox::printAtoms(){
	double * p;
	std::cout << "AtomPrint:: " << std::endl;
	for ( long i = 0; i < nAtoms ; i ++){
		p = atoms[i].getPos();
	std::cout <<  " atom " << i << " : Pos(x,y,z) " << p[0]/size[0] << " " << p[1]/size[1] << " " << p[2]/size[2] << std::endl;

	}

}

Atom * AtomBox::getAtom(long i){
	return (atoms + i);
}

bool sortDistAscending(ABoxNeighbor a, ABoxNeighbor b) { return SQR(a.coord[0])+SQR(a.coord[1])+SQR(a.coord[2]) < SQR(b.coord[0])+SQR(b.coord[1])+SQR(b.coord[2]); }

void AtomBox::srtNeighbors(){
	std::sort(neighbors, neighbors + nNeighbors, sortDistAscending);
}

const double * AtomBox::getOrigin() const{
	return origin;
}

const double * AtomBox::getSize() const{
	return size;
}

long AtomBox::getNumAtoms(){
	return nAtoms;
}

long AtomBox::getNumNeighbors(){
	return nNeighbors;
}

ABoxNeighbor * AtomBox::getNeighbors(){
	return neighbors;
}

double AtomBox::sqrDist(const double * p1,const double * p2) const{
	return SQR(p1[0] - p2[0]) + SQR(p1[1] - p2[1]) + SQR(p1[2] - p2[2]);
}

long AtomBox::getAtomNum(const Atom* atom) const {
	return atom - atoms;
}

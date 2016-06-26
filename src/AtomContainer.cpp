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

#include "AtomContainer.h"
#define DEFAULT_ANGULARTHRESHOLD 0.5e-2//1degree
#define DEFAULT_ORICAPACITY 10000 //320KB reserved as default for orientations to reduce the frequency of reallocations - critical, slow operation
#define DEFAULT_MEANORILEAFSIZE 10

AtomContainer::AtomContainer(double inMinBoxSize){
	numberAtoms = 0;
	minBoxSize = inMinBoxSize;
	nY = 0;
	nX = 0;
	nZ = 0;
	nXY = 0;
	nBoxes = 0;
}

AtomContainer::AtomContainer(double  * inCenter, double * inSize, unsigned long * inFragmentation){
	init(inCenter, inSize, inFragmentation, DEFAULT_ORICAPACITY);
}

AtomContainer::AtomContainer(double * inCenter, double * inSize, unsigned long * inFragmentation, unsigned long initOriCapacity){
	init(inCenter, inSize, inFragmentation, initOriCapacity);
}


void AtomContainer::init(double * inCenter, double * inSize, unsigned long * inFragmentation, unsigned long initOriCapacity ){
	numberAtoms = 0;
	nX = inFragmentation[0];//number of boxes on a line
	nY = inFragmentation[1];
	nZ = inFragmentation[2];
	nXY = nX * nY;//number of boxes on a face
	nBoxes = nXY * nZ;
	for (char i = 0; i < DIM; i++){
			size[i] = inSize[i];
			origin[i] = inCenter[i] - .5 * size[i];

	}
	boxes = new AtomBox[nBoxes];
	initBoxes();
	orient = new Orientator(boxes, initOriCapacity);
	grains = new GrainIdentificator(orient,boxes,nBoxes,DEFAULT_ANGULARTHRESHOLD,12);
}

void AtomContainer::destruct() {

}

AtomContainer::~AtomContainer() {
	delete orient;
	delete grains;
	if ( nBoxes > 0) delete [] boxes;
}

void AtomContainer::initBoxes(){
	ABoxNeighbor * neighbors;
	long ixyz, nxyz;
	long nid;
	long nx, ny, nz;
	boxSize[0] = size[0]/nX;
	boxSize[1] = size[1]/nY;
	boxSize[2] = size[2]/nZ;
//#pragma omp parallel private(neighbors, ixyz, nxyz, nid, nx, ny, nz) default(none)
{
	neighbors = new ABoxNeighbor[MOORE];
	double boxOrigin[DIM];
//#pragma omp for collapse(3)
	for(long iz = 0; iz < nZ; iz ++){
		for(long iy = 0; iy < nY; iy ++){
			for(long ix = 0; ix < nX; ix ++){
			ixyz = id(ix,iy,iz);
			nid = 0;
			//moore-neighborhood
			for (long nz = iz-1; nz <= iz+1; nz ++){
			for (long ny = iy-1; ny <= iy+1; ny ++){
			for (long nx = ix-1; nx <= ix+1; nx ++){
				//the box itself is not considered in moore-neighborhood
				if (nx != ix || ny != iy || nz != iz){
					//proof that neighbor lays inside the container
					if( valid(nx, ny, nz) ) {
					nxyz = id(nx, ny, nz);
					neighbors[nid].box = boxes + nxyz ;
					neighbors[nid].coord[0] = nx - ix;
					neighbors[nid].coord[1] = ny - iy;
					neighbors[nid].coord[2] = nz - iz;
					nid ++;
					}
				}
			}}}

			//all neighbors defined
			boxOrigin[0] = origin[0] + boxSize[0] * ix;
			boxOrigin[1] = origin[1] + boxSize[1] * iy;
			boxOrigin[2] = origin[2] + boxSize[2] * iz;
			boxes[ixyz].init(boxOrigin, boxSize, neighbors, nid);
			boxes[ixyz].srtNeighbors();
			}
		}
	}

	delete [] neighbors;
}//end parallel
}
void AtomContainer::setSize(double * inSize){
	size[0] = inSize[0];
	size[1] = inSize[1];
	size[2] = inSize[2];
}

void AtomContainer::setOrigin(double * inOrigin){
	origin[0] = inOrigin[0];
	origin[1] = inOrigin[1];
	origin[2] = inOrigin[2];
}
void AtomContainer::reducedPoint(const double* p, double* redP) const{
	redP[0] = reducedCoordinate(p[0],0);
	redP[1] = reducedCoordinate(p[1],1);
	redP[2] = reducedCoordinate(p[2],2);
}

double AtomContainer::reducedCoordinate(double pos, unsigned char dimension) const{
	if (!isPeriodic() || dimension > 2){
		return pos;
	}
	double relPos = pos - origin[dimension];
	//rest = position in [-1,1]*size
	double rest = fmod(relPos, size[dimension]);
	//result is position in [0,1]*size + origin
	if (rest >= 0.){
		return origin[dimension] + rest;
	}
	return origin[dimension] + rest + size[dimension];
}

void AtomContainer::calculateAtomOrientations(const double rSqrMin, const double rSqrMax){
	std::cout << "Thread " << omp_get_thread_num() << ": Orientation Calculation started. " << std::endl;
	AtomBoxP box;
//#pragma omp parallel default(shared) private(box)
//#pragma omp for
	long tenPercentNum = nBoxes/10;
	for (long iBox = 0; iBox < nBoxes; iBox ++){
		box = boxes + iBox;
		box->calculateAtomOrientationsFCC(iBox, DEFAULT_ANGULARTHRESHOLD,orient, rSqrMin, rSqrMax);
		if(iBox % tenPercentNum == 0){
#pragma omp critical
{
	std::cout << "Thread " << omp_get_thread_num() << ": Orientation Calculation finished " << iBox/tenPercentNum *10 << " % " << std::endl;
}
		}
	}
#pragma omp critical
{
	std::cout << LINE << "\n"
	<<"Thread " << omp_get_thread_num() <<": Orientation Calculation Done (1/3)\n"
	<< LINE <<std::endl;
}
}

const AtomBox * AtomContainer::getBoxes() const{
	return boxes;
}

long AtomContainer::getNumBoxes() const{
	return nBoxes;
}

void AtomContainer::addAtoms(double * inPos, long nAtoms){
	double * iPos;
//#pragma omp parallel shared(cout, inPos, nAtoms) private(iPos) default(none)
//#pragma omp for
	for (long i  = 0; i < nAtoms; i++){
		iPos = inPos + i*DIM ;
		if(!addAtom(iPos)) {
			std::cout << "WARNING: Adding atom " << i <<" failed (placed outside container)" << std::endl;
		}
	}
}

bool AtomContainer::addAtom(double * inPos){
	long ix, iy, iz;
	if(inPos[0] < origin[0] ) return false;
	if(inPos[1] < origin[1] ) return false;
	if(inPos[2] < origin[2] ) return false;
	ix = (inPos[0]-origin[0])/boxSize[0];
	iy = (inPos[1]-origin[1])/boxSize[1];
	iz = (inPos[2]-origin[2])/boxSize[2];
	if (!valid(ix,iy,iz)) return false;
	AtomBoxP box = &boxes[id(ix,iy,iz)];
	double * boxOrigin = box->getOrigin();
	double boxPos[DIM] =	{ inPos[0] - boxOrigin[0], inPos[1] - boxOrigin[1], inPos[2] - boxOrigin[2] };
	box->addAtom(boxPos);
	numberAtoms ++;
	return true;
}



void AtomContainer::outBox(long id){
	std::cout << "BOX " << id << " with address = " << &boxes[id] << std::endl;
	std::cout << " has " << boxes[id].getNumNeighbors() << " Neighbors and " << boxes[id].getNumAtoms() << " Atoms" << std::endl;
	double * boxSize = boxes[id].getSize();
	double * boxOrigin = boxes[id].getOrigin();
	for (char i=0; i < DIM; i++){
		std::cout << "boxSize[" << (int)i << "] = " << boxSize[i] << "   ";
		std::cout << "boxOrigin[" << (int)i << "] = " << boxOrigin[i] << std::endl;
	}
}

void AtomContainer::receiveAtomData(double * pos, long * boxIds, long * atomIds, oID * orientIDs, gID * grainIDs){
	receiveAtomData(pos,orientIDs,grainIDs);
	AtomBox * box;
	long iA;
	for (long iB = 0; iB < nBoxes; iB++){
			box = boxes + iB;
			for(iA = 0; iA < box->getNumAtoms(); iA++){
			*atomIds = iA;
			*boxIds = iB;
			atomIds ++;
			boxIds ++;
			}
		}
}

void AtomContainer::receiveAtomData(double * pos, oID * orientIDs, gID * grainIDs){
	long iA;
	AtomBox * box;
	Atom * a;
	double * posCur;
	double * boxPos;
	double * p = pos;
	oID * o = orientIDs;
	gID * g = grainIDs;
	for (long iB = 0; iB < nBoxes; iB++){
		box = boxes + iB;
		boxPos = box->getOrigin();
		for(iA = 0; iA < box->getNumAtoms(); iA++){
			a = box->getAtom(iA);
			posCur = a->getPos();
			p[0] = boxPos[0] + posCur[0];
			p[1] = boxPos[1] + posCur[1];
			p[2] = boxPos[2] + posCur[2];
			*o = a->getOrientationId();
			*g = a->getGrainId();
			p+=DIM;
			o++;
			g++;
		}
	}

}

long AtomContainer::getNumOrientations() const{
	return orient->getNumOrientations();
}

const Orientation * AtomContainer::getOrientations() const{
	return orient->getOrientations();
}
long AtomContainer::getNumGrains() const{
	return grains->getNumGrains();
}

const Grain * AtomContainer::getGrain(long  grainNum) const{
	return grains->getGrain(grainNum);
}

void AtomContainer::identifyGrains(double angularThreshold, double rSqrMin, double rSqrMax)
{
	grains->setSearchRadiiSquared(rSqrMin,rSqrMax);
	grains->run(angularThreshold);
	grains->assignOrphanAtoms();
	grains->sort();
}

const Orientator * AtomContainer::getOrientator() const {
	return orient;
}

void AtomContainer::generate(long nAtoms)
{
	nX = size[0]/minBoxSize;
	nY = size[1]/minBoxSize;
	nZ = size[2]/minBoxSize;
	if(nX > MAXFRAGMENT) nX = MAXFRAGMENT;
		else if(nX <= 0) nX = 1;
	if(nY > MAXFRAGMENT) nY = MAXFRAGMENT;
		else if(nY <= 0) nY = 1;
	if(nZ > MAXFRAGMENT) nZ = MAXFRAGMENT;
		else if(nZ <= 0) nZ = 1;
	nXY = nX*nY;
	nBoxes = nXY*nZ;
	boxes = new AtomBox[nBoxes];
	initBoxes();
	orient = new Orientator(boxes, nAtoms);
	grains = new GrainIdentificator(orient,boxes,nBoxes,DEFAULT_ANGULARTHRESHOLD,12);
}

const double * AtomContainer::getSize() const {
	return size;
}

long AtomContainer::getNumAtoms() const{
	return numberAtoms;
}

const double * AtomContainer::getOrigin() const{
	return origin;
}

long AtomContainer::id(long ix, long iy, long iz) const{
	return iz * nXY + iy * nX + ix;

}

bool AtomContainer::valid(long ix, long iy, long iz) const{
	if (ix < 0) return false;
	if (iy < 0) return false;
	if (iz < 0) return false;
	if (ix >= nX) return false;
	if (iy >= nY) return false;
	if (iz >= nZ) return false;
	return true;
}
PeriodicAtomContainer::PeriodicAtomContainer(double * inCenter, double * inSize, unsigned long * inFragmentation, unsigned long initOriCapacity){
	init(inCenter, inSize, inFragmentation, initOriCapacity);
}

PeriodicAtomContainer::PeriodicAtomContainer(double * inCenter, double * inSize, unsigned long * inFragmentation){
	init(inCenter, inSize, inFragmentation, DEFAULT_ORICAPACITY);
}

PeriodicAtomContainer::~PeriodicAtomContainer() {
}

void PeriodicAtomContainer::addAtoms(double * inPos, long nAtoms){
	double * iPos;
//#pragma omp parallel shared(cout, inPos, nAtoms) private(iPos) default(none)
//#pragma omp for
	for (long i  = 0; i < nAtoms; i++){
		iPos = inPos +i*DIM;
		addAtom(iPos);
	}
}

void PeriodicAtomContainer::initBoxes(){
	long ixyz, nxyz;
	long nid;
	long nx, ny, nz;//neighbor ids
	ABoxNeighbor * neighbors;
	boxSize[0] = size[0]/nX;
	boxSize[1] = size[1]/nY;
	boxSize[2] = size[2]/nZ;
//#pragma omp parallel private(neighbors, ixyz, nxyz, nid, nx, ny, nz) shared(cout) default(none)
{
	neighbors = new ABoxNeighbor[MOORE];
	double boxOrigin[DIM];
	//loop over all boxes
//#pragma omp for collapse(3)
		for(long iz = 0; iz < nZ; iz ++){
			for(long iy = 0; iy < nY; iy ++){
				for(long ix = 0; ix < nX; ix ++){
				ixyz = id(ix,iy,iz);
				nid = 0;
				//moore-neighborhood
				for (long nz = iz-1; nz <= iz+1; nz ++){
				for (long ny = iy-1; ny <= iy+1; ny ++){
				for (long nx = ix-1; nx <= ix+1; nx ++){
					//the box itself is not considered in moore-neighborhood
					if (nx != ix || ny != iy || nz != iz){
						//for periodic case all neighbors are forced to be accessible
						nxyz = id(nx, ny, nz);
						neighbors[nid].box = &boxes[nxyz];
						neighbors[nid].coord[0] = nx - ix;
						neighbors[nid].coord[1] = ny - iy;
						neighbors[nid].coord[2] = nz - iz;
						nid ++;
					}
				}}}
				//all neighbors defined
				boxOrigin[0] = origin[0] + boxSize[0] * ix;
				boxOrigin[1] = origin[1] + boxSize[1] * iy;
				boxOrigin[2] = origin[2] + boxSize[2] * iz;
				boxes[ixyz].init(boxOrigin, boxSize, neighbors, nid);
				boxes[ixyz].srtNeighbors();
				}
			}
		}
	delete [] neighbors;
}
}

bool PeriodicAtomContainer::addAtom(double * inPos){
	long ix, iy, iz;
	long iX, iY, iZ;
	//relative Position of the atom to the container origin
	double relPos[DIM] = {inPos[0]-origin[0], inPos[1]-origin[1], inPos[2]-origin[2]};
	ix = relPos[0]/boxSize[0];
	if(relPos[0] < 0.) ix --;
	iy = relPos[1]/boxSize[1];
	if(relPos[1] < 0.) iy --;
	iz = relPos[2]/boxSize[2];
	if(relPos[2] < 0.) iz --;
	//
	iX = relPos[0]/size[0];
	iY = relPos[1]/size[1];
	iZ = relPos[2]/size[2];

	AtomBoxP box = boxes +id(ix,iy,iz);
	double * boxOrigin = box->getOrigin();
	double boxPos[DIM] = {inPos[0] - boxOrigin[0], inPos[1] - boxOrigin[1], inPos[2] - boxOrigin[2]};
	//atom position inside the boxes are defined in the box's local coordinate system
	//hence, we need to translate position by periodic translation vectors in order to obtain
	//true relative coordinates (smaller as the box size)
	if( iX > 0 ) boxPos[0] -= iX * size[0];
	else if (relPos[0] < 0. ) boxPos[0] -= (iX - 1) * size[0];
	//
	if( iY > 0 ) boxPos[1] -= iY * size[1];
	else if (relPos[1] < 0. ) boxPos[1] -= (iY - 1) * size[1];
	//
	if( iZ > 0 ) boxPos[2] -= iZ * size[2];
	else if (relPos[2] < 0. ) boxPos[2] -= (iZ - 1) * size[2];
	//add atom to box
	box->addAtom(boxPos);
	numberAtoms++;
	return true;
}

long PeriodicAtomContainer::id(long ix, long iy, long iz) const{
	if(ix >= nX) ix = ix % nX;
	else if (ix < 0) ix = nX + ((ix+1) % nX - 1);

	if(iy >= nY) iy = iy % nY;
	else if (iy < 0) iy = nY + ((iy+1) % nY - 1);

	if(iz >= nZ) iz = iz % nZ;
	else if (iz < 0) iz = nZ + ((iz+1) % nZ - 1);

	return iz * nXY + iy * nX + ix;
}

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
	capacity = initOriCapacity;
	boxes = new AtomBox[nBoxes];
	initBoxes();
	orient = new Orientator(boxes, capacity);
	grains = new GrainIdentificator(orient,boxes,nBoxes,DEFAULT_ANGULARTHRESHOLD,12);
	atomInputOrder.reserve(capacity, nBoxes);
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
	neighbors = new ABoxNeighbor[MOORE];
	double boxOrigin[DIM];
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
}
void AtomContainer::setSize(const double * inSize){
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
	AtomBoxP box;
	long tenPercentNum = nBoxes/10;
	for (long iBox = 0; iBox < nBoxes; iBox ++){
		box = boxes + iBox;
		box->calculateAtomOrientationsFCC(DEFAULT_ANGULARTHRESHOLD,orient, rSqrMin, rSqrMax);
		if(iBox % tenPercentNum == 0){
#pragma omp critical
{
	std::cout << "Thread " << omp_get_thread_num() << ": Orientation Calculation finished " << iBox/tenPercentNum *10 << " % " << std::endl;
}
		}
	}
}

const AtomBox * AtomContainer::getBoxes() const{
	return boxes;
}

long AtomContainer::getNumBoxes() const{
	return nBoxes;
}

bool isFloat(const std::string & value){
	//returns true if value contains one of '.','e' or 'E'
	return (value.find_last_of(".eE") != std::string::npos);
}

void AtomContainer::addAtomProperty(const std::string& name) {
	atomPropertyList.addProperty(name, capacity);
}
void AtomContainer::addAtom(const double * pos, const std::vector<std::string>& atomProperties){
	int iProperty;
	//after adding the atom was successful, continue with adding additional properties
	if(addAtom(pos)) {
		//save additional properties
		for (iProperty = 0; iProperty < atomProperties.size(); iProperty++){
			if (!atomPropertyList.isPropertyInt(iProperty)) {
				//property has float type
				atomPropertyList.addFloatPropertyValue(iProperty, atof(atomProperties[iProperty].c_str()));
			} else {
				//property has integer type
				if (isFloat(atomProperties[iProperty])) {
					//but current value is of float type -> conversion
					atomPropertyList.convertPropertyToFloat(iProperty);
					atomPropertyList.addFloatPropertyValue(iProperty, atof(atomProperties[iProperty].c_str()));
				} else {
					//value has also integer type -> add value
					atomPropertyList.addIntPropertyValue(iProperty,atoi(atomProperties[iProperty].c_str()));
				}
			}
		}
	} else {
		std::cout << "WARNING: Adding atom " << getNumAtoms() <<" failed (placed outside container)" << std::endl;
	}
}

void AtomContainer::addAtoms(double * inPos, long nAtoms){
	double * iPos;
	for (long i  = 0; i < nAtoms; i++){
		iPos = inPos + i*DIM ;
		if(!addAtom(iPos)) {
			std::cout << "WARNING: Adding atom " << i <<" failed (placed outside container)" << std::endl;
		}
	}
}

long AtomContainer::getBoxId(const AtomBox* box) const {
	return box - boxes;
}

double AtomContainer::getAtomsProperty(int propertyNum, const AtomID& atomId) {
	long atomNum = atomInputOrder.getAtomNum(atomId);
	if (atomNum >= 0) {
		if (atomPropertyList.isPropertyInt(propertyNum)) {
			//is int
			return static_cast<double> (atomPropertyList.getIntPropertyValue(propertyNum, atomNum));
		} else {
			//is float
			return atomPropertyList.getFloatPropertyValue(propertyNum, atomNum);
		}
	}
	return 0.0;
}

void AtomContainer::calculateGrainProperties() {
	AtomID id;
	std::vector<std::vector<double>> grainProperties(grains->getNumGrains());
	std::vector<int > nums (grains->getNumGrains(),0);
	//initialize with zeroes
	for (long iG = 0; iG < grains->getNumGrains(); iG++){
		grainProperties[iG].resize(atomPropertyList.getNumProperties(),0.0);
	}
	AtomBox * box;
	gID grainId;
	int iP;
	//for each grain sum up all values for each corresponding atom
	for (id.iB = 0; id.iB < nBoxes; id.iB++){
		box = boxes + id.iB;
		for (id.iA = 0; id.iA < box->getNumAtoms(); id.iA++){
			grainId = box->getAtom(id.iA)->getGrainId();
			if (grainId >= 0) {
				for (iP = 0; iP < atomPropertyList.getNumProperties(); iP++){
					grainProperties[grainId][iP] += getAtomsProperty(iP,id);
				}
			}
		}
	}
	//calculate the average value
	Grain * grain;
	for (long iG = 0; iG < grains->getNumGrains(); iG++){
		grain = grains->getGrain(iG);
		for (iP = 0; iP < atomPropertyList.getNumProperties(); iP++){
				grainProperties[iG][iP] /= grain->getNumberOfAtoms();
		}
		grain->setProperties(grainProperties[iG]);
	}
}

bool AtomContainer::addAtom(const double * inPos){
	long ix, iy, iz;
	AtomID atomId;
	if(inPos[0] < origin[0] ) return false;
	if(inPos[1] < origin[1] ) return false;
	if(inPos[2] < origin[2] ) return false;
	ix = (inPos[0]-origin[0])/boxSize[0];
	iy = (inPos[1]-origin[1])/boxSize[1];
	iz = (inPos[2]-origin[2])/boxSize[2];
	if (!valid(ix,iy,iz)) return false;
	//retrieve and save the box-id
	atomId.iB = id(ix,iy,iz);
	//get the pointer to the box
	AtomBoxP box = boxes + atomId.iB;
	//obtain read access to the box's origin
	const double * boxOrigin = box->getOrigin();
	//calculate the position relative to the box's origin
	double boxPos[DIM] =	{ inPos[0] - boxOrigin[0], inPos[1] - boxOrigin[1], inPos[2] - boxOrigin[2] };
	//the corresponding atom-id is equal to the number of atoms
	//before adding the atom to the box
	atomId.iA = box->getNumAtoms();
	box->addAtom(boxPos);
	//remember the order the atoms were put into the container
	atomInputOrder.add(atomId);
	numberAtoms ++;
	return true;
}



void AtomContainer::outBox(long id){
	std::cout << "BOX " << id << " with address = " << &boxes[id] << std::endl;
	std::cout << " has " << boxes[id].getNumNeighbors() << " Neighbors and " << boxes[id].getNumAtoms() << " Atoms" << std::endl;
	const double * boxSize = boxes[id].getSize();
	const double * boxOrigin = boxes[id].getOrigin();
	for (char i=0; i < DIM; i++){
		std::cout << "boxSize[" << (int)i << "] = " << boxSize[i] << "   ";
		std::cout << "boxOrigin[" << (int)i << "] = " << boxOrigin[i] << std::endl;
	}
}

void AtomContainer::getAtomPropertyNames(std::vector<std::string> & outProperties, bool withDefaults) const{
	int i;
	outProperties.clear();
	if(withDefaults){
		//add default properties
		for ( i = 0; i < numDefaultProperties; i++){
			outProperties.push_back(defaultAtomProperties[i]);
		}
	}
	//add all remaining properties
	for ( i = 0; i < atomPropertyList.getNumProperties(); i++){
		outProperties.push_back(atomPropertyList.getPropertyName(i));
	}
}

void AtomContainer::getGrainAvgPropertyNames(std::vector<std::string> & outProperties) const{
	//add all grain average properties
	for (int i = 0; i < atomPropertyList.getNumProperties(); i++){
		outProperties.push_back(atomPropertyList.getPropertyName(i) + "_grainAvg");
	}
}

void AtomContainer::getAtomsPosition(long atomNum, double * outPos) const{
	const AtomID * id;
	AtomBox * box;
	Atom * a;
	const double * boxPos;
	const double * curPos;
	if(atomNum >=  0 && atomNum < numberAtoms){
		id = atomInputOrder.getAtomId(atomNum);
		box = boxes + id->iB;
		a =  box->getAtom( id->iA );
		boxPos = box->getOrigin();
		curPos = a->getPos();
		//calculate global coordinates from local coordinates and origin
		outPos[0] = boxPos[0] + curPos[0];
		outPos[1] = boxPos[1] + curPos[1];
		outPos[2] = boxPos[2] + curPos[2];
	}
}

void AtomContainer::getAtomsProperties(long atomNum, std::vector<std::string>& outProperties, bool withDefaults, bool withGrainAvg) const{
	outProperties.clear();
	const AtomID * id = atomInputOrder.getAtomId(atomNum);
	AtomBox * box = boxes + id->iB;
	Atom * a =  box->getAtom( id->iA );
	if(atomNum >=  0 && atomNum < numberAtoms){
		if (withDefaults){
			//add default properties
			outProperties = {
				std::to_string(id->iB),
				std::to_string(id->iA),
				std::to_string(a->getOrientationId()),
				std::to_string(a->getGrainId())
			};
		}
		for (int i = 0; i < atomPropertyList.getNumProperties(); i++){
			//add all remaining properties
			if (atomPropertyList.isPropertyInt(i)) {
				//property is int
				outProperties.push_back(std::to_string(atomPropertyList.getIntPropertyValue(i,atomNum)));
			} else	{ // property is float
				outProperties.push_back(ori::to_string(atomPropertyList.getFloatPropertyValue(i,atomNum)));
			}
		}
		if(withGrainAvg) {
			//add all grain avg properties
			if (a->getGrainId() != NO_GRAIN){
				const std::vector<double> & grainAvgProperties =
						grains->getGrain(a->getGrainId())->getProperties();
				for (int i = 0; i < grainAvgProperties.size(); i++){
					outProperties.push_back(ori::to_string(grainAvgProperties[i]));
				}
			} else {
				//for atoms not belonging to any grain set the grain-avg-property equal to zero
				for (int i = 0; i < atomPropertyList.getNumProperties(); i++){
					outProperties.push_back("0");
				}
			}
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
	calculateGrainProperties();
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

PeriodicAtomContainer::~PeriodicAtomContainer() {
}

void PeriodicAtomContainer::addAtoms(double * inPos, long nAtoms){
	double * iPos;
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
	neighbors = new ABoxNeighbor[MOORE];
	double boxOrigin[DIM];
	//loop over all boxes
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

bool PeriodicAtomContainer::addAtom(const double * inPos){
	long ix, iy, iz;
	long iX, iY, iZ;
	AtomID atomId;
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
	//obtain and save the id of the box
	atomId.iB = id(ix,iy,iz);
	//get the corresponding pointer to the box
	AtomBoxP box = boxes + atomId.iB;
	const double * boxOrigin = box->getOrigin();
	//calculate the relative coordinates to the box's origin
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
	//the corresponding atom-id is equal to the number of atoms
	//before adding the atom to the box
	atomId.iA = box->getNumAtoms();
	//add atom to box
	box->addAtom(boxPos);
	//remember the order the atoms were put into the container
	atomInputOrder.add(atomId);
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

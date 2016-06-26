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

#ifndef ATOMCONTAINER_H_
#define ATOMCONTAINER_H_
#include "AtomBox.h"
#include "Orientator.h"
#include "Grain.h"
#include "GrainIdentificator.h"
#define MAXFRAGMENT 80

class AtomContainer {
public:
	AtomContainer(){}
	AtomContainer(double inMinBoxSize);
	AtomContainer(double * center, double * size, unsigned long * fragmentation);
	AtomContainer(double * center, double * size, unsigned long * fragmentation,  unsigned long initCapacity);
	virtual ~AtomContainer();
	void destruct();
	virtual void generate(long nAtoms);
	void setSize(double * inSize);
	const double * getSize() const;
	const double * getOrigin() const;
	void setOrigin(double * inOrigin);
	void calculateAtomOrientations(const double rSqrMin, const double rSqrMax);
	void identifyGrains(double angularThreshold, double rSqrMin, double rSqrMax);
	virtual void addAtoms(double * atomPos, long nAtoms);
	void outBox(long id);
	void reducedPoint(const double* p, double* redP) const;
	long getNumAtoms() const;
	long getNumBoxes() const;
	long getNumBoxesInX() const { return nX;};
	long getNumBoxesInY() const { return nY;};
	long getNumBoxesInZ() const { return nZ;};
	const AtomBox * getBoxes() const;
	long getNumOrientations() const;
	const Orientation * getOrientations() const;
	const Orientator * getOrientator() const;
	long getNumGrains() const;
	const Grain * getGrain(long grainNum) const;
	virtual bool isPeriodic() const {return false;}
	void receiveAtomData(double * pos, long * boxIds, long * atomIds, oID * orientIDs, gID * grainIDs);
	void receiveAtomData(double * pos, oID * orientIDs, oID * meanOrientIDs);
protected:
	void init(double * center, double * size, unsigned long * fragmentation, unsigned long initCapacity);
	double reducedCoordinate(double pos, unsigned char dimension) const;
	double origin[DIM];
	double boxSize[DIM];
	double size[DIM];
	long nX = 0, nY = 0, nZ = 0;
	long nBoxes = 0, nXY = 0;
	long numberAtoms = 0;
	double minBoxSize;
	AtomBox * boxes = nullptr;
	Orientator * orient = nullptr;
	GrainIdentificator * grains = nullptr;
private:
	virtual void initBoxes();
	virtual inline bool addAtom(double * pos);
	virtual inline bool valid(long ix, long iy, long iz) const;
	virtual inline long id(long ix, long iy, long iz) const;
};
//
class PeriodicAtomContainer : public AtomContainer{
public:
	PeriodicAtomContainer(double inMinBoxSize) : AtomContainer(inMinBoxSize){}
	PeriodicAtomContainer(double * center, double * size, unsigned long * fragmentation);
	PeriodicAtomContainer(double * center, double * size, unsigned long * fragmentation,  unsigned long initCapacity);
	~PeriodicAtomContainer();
	void addAtoms(double * atomPos, long nAtoms);
	bool isPeriodic() const {return true;}
private:
	void initBoxes();
	inline bool addAtom(double * pos);
	inline long id(long ix, long iy, long iz) const;
};

#endif /* ATOMCONTAINER_H_ */

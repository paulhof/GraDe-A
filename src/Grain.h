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

#ifndef GRAIN_H_
#define GRAIN_H_
#include "GradeA_Defs.h"
#include "Orientator.h"
#include "MeanOrientation.h"
#include "CubicLattices.h"
class Grain {
public:
	Grain();
	void add(Atom * atom, const Orientator * orient);
	void addOrphan(Atom * oAtom);
	void addToCenter(const double * vec);
	void calcAverageCenter();
	void calcTotalOrientationSpread(const Orientator * orient);
	void calcRegularOrientationSpread(const Orientator * orient);
	void calcOrphanOrientationSpread(const Orientator * orient);
	void setProperties(const std::vector<double> & properties);
	const std::vector<double>& getProperties() const;
	void recalculateMeanOrientation();
	double orientationSpread() const;
	double cosHalfOrientationSpread() const;
	void unAssign();
	void assign(gID grainId);
	void assignRegularAtoms(gID grainId);
	void assignOrphanAtoms(gID grainId);
	gID getAssignedId() const;
	double getVolume(const CubicLattice & material) const;
	double getVolumeInLatticeUnit(const CubicLattice & material) const;
	long getNumberOfAtoms() const;
	long getNumberOfRegularAtoms() const;
	long getNumberOfOrphanAtoms() const;
	const double * getPosition() const;
	const Orientation * getOrientation() const;
	virtual ~Grain();
private:
	std::vector<Atom*> atoms;
	std::vector<Atom*> orphanAtoms;
	std::vector<double> meanProperties;
	double center[DIM];
	MeanOrientation meanOrient;
	double oriSpread = 0., cosHalfOriSpread = 1.;//in rad
	double orphanOriSpread = 0., orphanCosHalfOriSpread = 1.;//in rad
	double totalOriSpread = 0., totalCosHalfOriSpread = 1.;
	long assignedId = NO_GRAIN;
	long nCenter = 0;
};

#endif /* GRAIN_H_ */

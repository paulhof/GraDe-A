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

#ifndef GRAINDATA_H_
#define GRAINDATA_H_
#include "GradeA_Defs.h"
#include "Orientator.h"
#include "Grain.h"
//! \brief A class used to store information usually referring to a Grain object.
//! Does not include information about corresponding atoms.
class GrainData {
public:
	GrainData(){};
	GrainData(double * inCenter, double * inQuaternion, long inNumAtoms,
			gID inAssignedId, long inNumRegularAtoms = 0, long inNumOprhanAtoms = 0,
			double inOrientationSpread = 0., const std::vector<double> & properties = {});
	GrainData(const Grain * inGrain);
	virtual ~GrainData();
	//!\brief returns a pointer to access the entries of center.
	const double * getCenter() const;
	//!\brief returns a pointer to meanOrient.
	const Orientation * getOrientation() const ;
	//!\brief returns the volume in A^3
	double getVolume(const CubicLattice& material) const;
	//!\brief returns the volume in the unit specified in material
	double getVolumeInLatticeUnit(const CubicLattice& material) const;
	//!\brief returns the number of atoms.
	long getNumberOfAtoms() const;
	//!\brief returns the number of regular atoms.
	long getNumberOfRegularAtoms() const;
	//!\brief returns the number of orphan atoms.
	long getNumberOfOrphanAtoms() const;
	//!\brief returns the orientation spread.
	double getOrientationSpread() const;
	//!\brief returns the assigned id.
	gID getAssignedId() const;
	//!\brief sets the assigned id.
	void setAssignedId(gID inAssignedId);
	//!\brief calculates misOrientation initOri and saves it into misOrientationToInit.
	void setMisOrientation(const Orientation * initOri);
	//!\brief sets distanceToInit
	void setDistanceToInit(double inDistance);

	double getDistanceToInit() const {
		return distanceToInit;
	}

	double getMisOrientationToInit() const {
		return misOrientationToInit;
	}

	double getRedMisOrientationToInit() const {
		return redMisOrientationToInit;
	}
	const std::vector<double> & getProperties() const;
private:
	double center[DIM] = {0.,0.,0.};
	Orientation meanOrient;
	long numAtoms = 0;
	long numOrphanAtoms = 0;
	long numRegularAtoms = 0;
	double orientationSpread = 0.;
	//
	double misOrientationToInit = 0.;
	double redMisOrientationToInit = 0.;
	double distanceToInit = 0.;
	gID assignedId = NO_GRAIN;
	std::vector<double> meanProperties;
};

#endif /* GRAINDATA_H_ */

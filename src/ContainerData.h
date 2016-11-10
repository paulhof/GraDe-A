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

#ifndef CONTAINERDATA_H_
#define CONTAINERDATA_H_
#include "AtomContainer.h"
#include "GradeA_Defs.h"
#include "GrainData.h"
//!\brief A class which contains important information usually referring to a  AtomContainer object.
class ContainerData {
public:
	ContainerData(bool inPeriodic = true){
		periodic = inPeriodic;
	};
	ContainerData(AtomContainer * con);
	virtual ~ContainerData();
	void setSize(double * inSize);
	void setPeriodic(bool inPeriodic);
	void setOrigin(double * inOrigin);
	//!\brief Calculates the squared distance of two points p1 and p2 in this container.
	double sqrDistance(const double * p1, const double * p2) const;
	//!\brief Returns the number of grains stored in grains.
	long getNumberOfGrains() const;
	//!\brief Adds a grain to grains.
	void addGrain (GrainData &grain);
	//!\brief If valid, returns a pointer to the grainId's entry inside grains.
	//! Returns nullptr if grainId greater or equal to size of grains.
	const GrainData * getGrain(gID grainId) const;
	const double * getSize() const;
	const double * getOrigin() const;
	void setAtomPropertyNames(const std::vector<std::string> & inProperties);
	void getAtomPropertyNames(std::vector<std::string> & outProperties) const;
	bool isPeriodic() const {
		return periodic;
	}
	gID maxGrainId() const;
	//!\brief calculates the reduced position to a given point p.
	//\param[in] p Pointer to input point data.
	//\param[out] redP Pointer to data, where result is written into.
	void reducedPoint(const double * p, double *redP) const;
private:
	//!\brief Calculates the reduced coordinate of a given dimension to pos.
	//! If the container is periodic, returns a value inside the containers boundaries.
	double reducedCoordinate(double pos, unsigned char dimension) const;
	std::vector<GrainData>  grains;
	std::vector<std::string> propertyNames;
	double origin[DIM] {
		0., 0., 0.
	};
	double size[DIM];
	bool periodic;
};

#endif /* CONTAINERDATA_H_ */

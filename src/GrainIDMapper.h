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

#ifndef GRAINIDMAPPER_H_
#define GRAINIDMAPPER_H_

#define NO_GRAIN -1
#define SMALLGRAINVOL 20000
#include "StopWatch.h"
#include "GradeA_Defs.h"
#include "ContainerData.h"
#include "io/CFGEditor.h"
#include "io/GrainCSVFileFormat.h"
#include "io/CSVTableWriter.h"

//!\brief A class which maps similar grains from two different datasets
class GrainIDMapper {
public:
	//! \brief A constructor, which setups a link to two already existing ContainerData datasets.
	GrainIDMapper(const CubicLattice & material, ContainerData * inOldContainer, ContainerData * inCurContainer, long inNewGrainIdBegin = 0);
	virtual ~GrainIDMapper();
	//! \brief Initializes values of the object.
	void init(double inMaxCosHalfMisOri, double inMaxVolFraction, long inNewGrainsIdBegin);
	//! \brief This function calculates the mapping.
	void map();
	//! \brief This method assigns the current grainIds to the given container object.
	//! Does not do anything if the AtomContainer object does not comprise enough grains.
	void assign(AtomContainer * inCurContainer);
	//
	void edit(CFGEditor * editor);
	//
	void print(CSVTableWriter * csvWriter);
	//! \brief Returns the mapped grain id to the current grain identified by curGrainNum.
	gID getMappedId(long curGrainNum) const;
	//! \brief Returns the mapped grain id to the current grain identified by its assigned-id before mapping.
	gID getMappedIdToOldAssignedId(gID oldAssignedId) const;
	long getNumCurGrains() const;
	long getNumPrevGrains() const;
	long getNewGrainId() const;

	const ContainerData * getCurContainer() const;
	const std::vector<gID> & getCurMappedIds() const;
	const std::vector<gID> & getCurOldAssignedIds() const;

	const CubicLattice * getMaterial() const {
		return material;
	}

private:
	//!\brief returns the grain id of a grain in oldContainer with similar properties to the given grain
	gID correspondingOldGrainId(const GrainData * grain);
	const CubicLattice * material = nullptr;
	ContainerData * prevContainer = nullptr;
	ContainerData * curContainer = nullptr;
	std::vector<gID> curMappedIds;
	std::vector<gID> curOldAssignedIds;
	std::vector<bool> prevGrainIsMapped;
	//! maximum grain-grain distance
	double maxSqrDistance;
	//! maximum grain-grain cosine misorientation
	double maxCosHalfMisOri;
	//! maximum grain-grain relative volume deviation
	double maxVolFraction;
	long numPrevGrains;
	long numCurGrains;
	//
	//! id, where numbering of unmapped grains begins
	long newGrainIdBegin;
	//! id, which will be assigned to any new grain.
	long newGrainId;
};

class GrainIDMapping {
public:
	GrainIDMapping(const GrainIDMapper * mapper);
	//! \brief This method assigns the current grainIds to the given container object.
	//! Does not do anything if the AtomContainer object does not comprise enough grains.
	void assign(AtomContainer * inCurContainer);
	//
	void edit(CFGEditor * editor);
	//
	void print(CSVTableWriter * csvWriter);
	//! \brief Returns the mapped grain id to the current grain identified by curGrainNum.
	gID getMappedId(long curGrainNum) const;
	//! \brief Returns the mapped grain id to the current grain identified by its assigned-id before mapping.
	gID getMappedIdToOldAssignedId(gID oldAssignedId) const;
private:
	const CubicLattice * material = nullptr;
	ContainerData curContainer;
	std::vector<gID> curMappedIds;
	std::vector<gID> curOldAssignedIds;
	long numCurGrains;

};
#endif /* GRAINIDMAPPER_H_ */

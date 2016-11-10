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
//!\brief A class, which allows to store atoms directly.
//!The box is a cuboid cell described by its origin and size.
//!For each box a neighborhood is defined, which is used for the underlying nearest-neighbor search functions.
class AtomBox {
public:
	//!\brief Initializes default values for the AtomBox object.
	AtomBox();
	virtual ~AtomBox();

	//!\brief Initializes the AtomBox object utilizing the given parameters.
	//!\param[in] origin Origin of the box. Must contain at least three accessible elements.
	//!\param[in] size Size of the box. Must contain at least three accessible elements.
	//!\param[in] neighbors List containing the neighbors of the box. Must contain at least nNeighbors elements.
	//!\param[in] nNeighbors Number of elements in the neighbors list.
	//!\param[in] sizeLinked Flag indicating whether to allocate data for the size or not.
	//!\c true indicates no allocation but the data \c size points to must be available during existence.
 	void init(const double * origin, const double * size, const ABoxNeighbor * neighbors,
 			long nNeighbors, bool sizeLinked = true);

 	//!\brief Sorts the neighbors by distance so that the nearest neighbors are listed first.
	void srtNeighbors();

	//!\brief Adds an atom to the box.
	//!\param[in] The position of the atom. At least three elements must be accessible.
	void addAtom(double * pos);

	//!\brief Calculates the orientations of the stored atoms.
	//!\param[in] angleThreshold
	//!\param[in,out] orient Orientator-object used to calculate and store the orientations.
	//!\param[in] rSqrMin Minimum squared radius used for nearest-neighbor search
	//!\param[in] rSqrMax Maximum squared radius used for nearest-neighbor search.
	void calculateAtomOrientationsFCC(double angleThreshold, Orientator  * orient, const double rSqrMin, const double rSqrMax);

	//!\brief Calculates the nearest neighbors to an atom that lay in the sphere-segment defined by an inner and outer radius.
	unsigned char atomNeighbors(const long atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, double * nborPositions);

	//!\brief Calculates the nearest neighbors to an atom that lay in the sphere-segment defined by an inner and outer radius.
	unsigned char atomNeighbors(const long atomId, const double rSqrMin, const double rSqrMax, const unsigned char nMaxAtomNeighbors, AtomBox ** outNborBoxesList, long * outNborAtomIdList, double * outNborPosList);

	//!\brief Tries to find \c nAtomNeighbors nearest neighbors to an atom that are closest to that atom.
	unsigned char nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, double * nborPositions);

	//!\brief Tries to find \c nAtomNeighbors nearest neighbors to an atom that are closest to that atom.
	unsigned char nearestAtomNeighbors(const long atomId, const unsigned char nAtomNeighbors, Atom ** nborAtoms);

	//!\return A pointer to the atom identified by \c atomId
	Atom * getAtom(long atomId);

	//!\return Returns the number an atom is stored.
	long getAtomNum(const Atom * atom) const;

	//!\brief Calculates the global position of an atom.
	//!\param[in] atomId Number used to identify an atom.
	//!\param[out] outPos Data to write the global position into. At least three elements must be accessible.
	void obtainGlobalAtomPos(long atomId, double * outPos) const;

	//!\brief Couts data for each atom in the box.
	void printAtoms();

	//!\return The origin of the box.
	const double * getOrigin() const;

	//!\return The size of the box.
	const double * getSize() const;

	//!\return The number of neighbors of the box.
	long getNumNeighbors();

	//!\return The number of atoms stored in the box.
	long getNumAtoms();

	//!\return The neighbors of the box.
	ABoxNeighbor * getNeighbors();
private:
	//!\return The distance between two points.
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

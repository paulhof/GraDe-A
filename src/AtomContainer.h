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
#include "AtomPropertyList.h"
#include "AtomIdList.h"
#define MAXFRAGMENT 80
//!\brief Container class inside which a whole atom-position configuration is stored.\n
//! An AtomContainer object represents a three-dimensional block which boundaries are defined by its origin and size.\n
//!	The class comprises of a cellular structure (i.e. boxes) which is used to store the atoms.\n
//! It provides methods to calculate the orientations of all atoms and to identify grains and locally stores all of the corresponding information.
//! In order to use a
class AtomContainer {
public:
	//!\brief Constructor which sets the minimum box size.
	//!\param[in] inMinBoxSize Minimum box-size to set.
	AtomContainer(double inMinBoxSize);

	//!\brief Constructor which fully initializes the container object with a default capacity.
	//!\param[in] center The center of the container in Angstrom. Must refer to data with at least three accessible elements.
	//!\param[in] size The size of the container in Angstrom. Must refer to data with at least three accessible elements.
	//!\param[in] fragmentation The number of boxes into which the container is split for each coordinate direction.
	//! Must refer to data with at least three accessible elements.
	AtomContainer(double * center, double * size, unsigned long * fragmentation);

	//!\brief Constructor which fully initializes the container object with an adjustable capacity.
	//!\param[in] center The center of the container in Angstrom. Must refer to data with at least three accessible elements.
	//!\param[in] size The size of the container in Angstrom. Must refer to data with at least three accessible elements.
	//!\param[in] fragmentation The number of boxes into which the container is split for each coordinate direction.
	//! Must refer to data with at least three accessible elements.
	//!\param[in] initCapacity Number of atoms for which memory is allocated.
	AtomContainer(double * center, double * size, unsigned long * fragmentation,  unsigned long initCapacity);
	virtual ~AtomContainer();
	//!\brief Adds an atom to the container.
	//!\param[in] pos The position of the atom. Must refer to data with at least three accessible elements.
	//!\param[in] atomProperties List containing the values of the additional (optional) properties of the atom.
	//! The number of elements of atomProperties must be consistent with the number of properties stored in atomPropertyList.
	void addAtom(const double * pos, const std::vector<std::string> & atomProperties);

	//!\brief Adds an additional user-defined property to the container.
	//! Increases the number of properties stored in atomPropertyList by 1.
	void addAtomProperty(const std::string & name);

	//!\brief Puts out the names of all currently considered atom-properties.
	//!\param[out] outProperties List into which the names of the properties are written.
	//!\param[in] withDefaults Defines whether the default properties should be included at the front of the list.
	void getAtomPropertyNames(std::vector<std::string> & outProperties, bool withDefaults = true) const;

	//!\brief Puts out the names of the grain-average properties.
	//! Each grain-average property corresponds to an user-defined atom-property, but has a different name for distinguishing.
	//!\param[out] outProperties List into which the names of the grain-average properties are written.
	void getGrainAvgPropertyNames(std::vector<std::string> & outProperties) const;

	//!\brief Puts out the property values of a single atom as a string vector.
	//!\param[in] atomNum The atom-number used to identify an atom. Must be inside the range 0 to getNumAtoms()-1.
	//!\param[out] outProperties List into which the property-values are written.
	//!\param[in] withDefaults Defines whether the values of the default properties should be included at the front of the list.
	//!\param[in] withGrainAvg Defines whether the values of the grain-average properties should be included at the end of the list.
	void getAtomsProperties(long atomNum, std::vector<std::string> & outProperties, bool withDefaults = true, bool withGrainAvg = true) const;

	//!\brief Puts out the position of an atom.
	//!\param[in] atomNum The atom-number used to identify an atom. Must be inside the range 0 to getNumAtoms()-1.
	//!\param[out] outPos Data into which the position is written. Must refer to a datablock with at least three accessible elements.
	void getAtomsPosition(long atomNum, double * outPos) const;

	//!\brief Returns a user-defined property-value of an atom.
	//!\param[in] propertyNum The number of the corresponding property. Must be in the range 0 to number of properties stored in atomPropertyList - 1.
	//!\param[in] atomId Structure to identify the atom of interest.
	double getAtomsProperty(int propertyNum, const AtomID & atomId);

	//!\brief Initializes the container object allocating memory for nAtoms atoms.
	//!\param[in] nAtoms Number of atoms to allocate memory for.
	virtual void generate(long nAtoms);

	//!\brief Sets the size of the container object (in Angstrom).
	//!\param[in] inSize Reference to the x,y and z coordinates of the size-values to set. The corresponding data-block must contain at least three accessible elements.
	void setSize(const double * inSize);

	//!\brief Returns a reference to the size coordinates of the container object (in Angstrom).
	const double * getSize() const;

	//!\brief Returns a reference to the origin coordinates of the container object (in Angstrom).
	const double * getOrigin() const;

	//!\brief Sets the origin of the container object (in Angstrom).
	//!\param[in] inOrigin Reference to the x, y and z coordinates of the origin-values to set. The corresponding data-block must contain at least three accessible elements.
	void setOrigin(double * inOrigin);

	//!\brief Tries to calculate an orientation for every atom stored inside the container object.
	//! Uses \c rSqrMin and \c rSqrMax for the nearest-neighbor-search in case \c NEAREST_NEIGHBORHOOD is not defined.
	//!\param[in] rSqrMin Minimum squared radius for the nearest-neighbor identification step (in Angstrom^2).
	//!\param[in] rSqrMax Maximum squared radius for the nearest-neighbor identification step (in Angstrom^2).
	void calculateAtomOrientations(const double rSqrMin, const double rSqrMax);

	//!\brief Determines all grain entities by grouping similar oriented atoms together.
	//! The orientations of the atoms must be calculated beforehand by the use of \c calculateAtomOrientations()
	//!\param[in] angularThreshold Maximum angular misorientation for the local criterion.
	//!The global criterion's misorientation-threshold is three times \c angularThreshold
	//! Uses \c rSqrMin and \c rSqrMax for the nearest-neighbor-search in case \c NEAREST_NEIGHBORHOOD is not defined.
	//!\param[in] rSqrMin Minimum squared radius for the nearest-neighbor identification step (in Angstrom^2).
	//!\param[in] rSqrMax Maximum squared radius for the nearest-neighbor identification step (in Angstrom^2).
	void identifyGrains(double angularThreshold, double rSqrMin, double rSqrMax);

	//!\brief Adds atoms to the container by a list of atom-positions.
	//!\param[in] atomPos Pointer to the beginning of the list. At least 3*nAtoms elements must be accessible.
	//!\param[in] nAtoms Number of atoms contained in the list.
	virtual void addAtoms(double * atomPos, long nAtoms);

	//!\brief Prints out the origins and sizes of all boxes into stdout.
	//! Deprecated / not in use.
	void outBox(long id);

	//!\brief Calculates the reduced position of a position.
	//!\param[in] p The position which reduced position will be calculated. At least 3 elements must be accessible.
	//!\param[out] redP The reduced position to \c p. At least 3 elements must be accessible.
	void reducedPoint(const double* p, double* redP) const;

	//!\return The number of atoms contained in the container object.
	long getNumAtoms() const;

	//!\return The number of boxes that make up the container.
	long getNumBoxes() const;

	//!\return The number of boxes that are arranged along the x-direction.
	long getNumBoxesInX() const { return nX;};
	//!\return The number of boxes that are arranged along the y-direction.
	long getNumBoxesInY() const { return nY;};
	//!\return The number of boxes that are arranged along the z-direction.
	long getNumBoxesInZ() const { return nZ;};

	//!\return A pointer to the underlying box-list. \c getNumBoxes() elements are accessible.
	const AtomBox * getBoxes() const;

	//!\return The corresponding box-index in the underlying box-list to \c box.
	//!\param[in] A pointer to a box. Must be contained in the underlying box-list.
	long getBoxId (const AtomBox * box) const;

	//!\return The number of orientations stored inside the container object.
	long getNumOrientations() const;

	//!\return A pointer to the underlying orientation-list. Contains \c getNumOrientations() elements.
	const Orientation * getOrientations() const;

	//!\return A pointer to the underlying orientator-object.
	const Orientator * getOrientator() const;

	//!\return The number of grains stored inside the container object.
	long getNumGrains() const;

	//!\return A pointer to an element in the underlying grain list.
	//!\param[in] grainNum Grain-number used to identify the underlying grain-object. Must be in the range 0 to \c getNumGrains().
	const Grain * getGrain(long grainNum) const;

	//!\return Whether the container has periodic boundary conditions or not.
	virtual bool isPeriodic() const {return false;}
protected:
	void init(double * center, double * size, unsigned long * fragmentation, unsigned long initCapacity);
	double reducedCoordinate(double pos, unsigned char dimension) const;
	void calculateGrainProperties();
	double origin[DIM];
	double boxSize[DIM];
	double size[DIM];
	long nX = 0, nY = 0, nZ = 0;
	long nBoxes = 0, nXY = 0;
	long numberAtoms = 0;
	double minBoxSize;
	long capacity = 0;
	AtomBox * boxes = nullptr;
	Orientator * orient = nullptr;
	GrainIdentificator * grains = nullptr;
	AtomIdList atomInputOrder;
	AtomPropertyList atomPropertyList;
	const static int numDefaultProperties = 4;
	std::string defaultAtomProperties[numDefaultProperties] =
		{BOX_ID_NAME, ATOM_ID_NAME, ORIENTATION_ID_NAME, GRAIN_ID_NAME};
private:
	virtual void initBoxes();
	virtual inline bool addAtom(const double * pos);
	virtual inline bool valid(long ix, long iy, long iz) const;
	virtual inline long id(long ix, long iy, long iz) const;
};
//
class PeriodicAtomContainer : public AtomContainer{
public:
	PeriodicAtomContainer(double inMinBoxSize) : AtomContainer(inMinBoxSize){};
	PeriodicAtomContainer(double * center, double * size, unsigned long * fragmentation) : AtomContainer(center, size, fragmentation){};
	PeriodicAtomContainer(double * center, double * size, unsigned long * fragmentation,  unsigned long initCapacity) : AtomContainer(center, size, fragmentation, initCapacity){};
	~PeriodicAtomContainer();
	void addAtoms(double * atomPos, long nAtoms);
	bool isPeriodic() const {return true;}
private:
	void initBoxes();
	inline bool addAtom(const double * pos);
	inline long id(long ix, long iy, long iz) const;
};

#endif /* ATOMCONTAINER_H_ */

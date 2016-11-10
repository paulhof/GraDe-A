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

#ifndef GRAINIDENTIFICATOR_H_
#define GRAINIDENTIFICATOR_H_

#define GRAINALLOC 5000

#include "GradeA_Defs.h"
#include "Orientator.h"
#include "AtomBox.h"
#include "Grain.h"
class RecursiveGrainIdentificationEngine;
typedef struct {
	long id;
	long count;
} Occurrence;

class GrainIdentificator {
public:
	GrainIdentificator(){};
	GrainIdentificator(Orientator * inOrient, AtomBox * inBoxes, long inNumBoxes, double angleThreshold, unsigned char nMaxAtomNeighbors);
	long getNumGrains();
	virtual ~GrainIdentificator();
	Grain * getGrain(gID grainID);
	void setSearchRadiiSquared(double inRsqrMin, double inRsqrMax);
	long run(double angularThreshold);//returns the number of found grains
	void calculateOrientationSpread();
	void assignOrphanAtoms(long depth = 0);
	void sort();//sorts grains with decreasing number of atoms (= volume) and calls the assign() routine
	void assign();//sets the grain ids of the atoms
private:
	inline void rebuildOrphanAtomList(std::vector<AtomID> & unassignedAtoms, std::vector<bool> &isUnassigned, long numUnassignedAtoms);
	inline void orphanAtomAssignAttempt(std::vector<AtomID> & unassignedAtoms, std::vector<bool> &isUnassigned, long &numUnAssignedAtoms);
	inline Occurrence sortAndFindMaxGrainOccurrence(Atom ** atomSet, unsigned char setSize) const;
	void newEmptyGrain();
	void deleteLastGrain();
	//extern data
	AtomBox * boxes = nullptr;
	Orientator * orient = nullptr;
	long numBoxes = 0;
	//intern data
	std::vector <Grain * > grains;
	long numGrains = 0;
	long grainCapacity = 0;
	long grainAlloc = 0;
	RecursiveGrainIdentificationEngine * engine = nullptr;
	unsigned char nMaxAtomNeighbors = 0;
};

class GrainCandidate{
public:
	GrainCandidate();
	GrainCandidate(AtomBox * inBox, long inAtomId, double * inPos);
	~GrainCandidate();
	void init(AtomBox * inBox, long inAtomId, double * inPos);
	void activate();
	void setParent(GrainCandidate * inParent);
	bool isActive();
	unsigned char findNeighbors(unsigned char inNumMaxAtomNeighbors, GrainCandidate * outNeighborsList, double rSqrMin, double rSqrMax);
	Atom * getAtom();
	Atom * getParentAtom();
	AtomBox * getBox();
	double * getPosition();
	long getAtomId();
private:
	double position[DIM];
	long atomId;
	bool active;
	AtomBox * box;
	Atom * parent;
};

class RecursiveGrainIdentificationEngine;
class GrainCandidateList{
public:
	GrainCandidateList(RecursiveGrainIdentificationEngine *inOwner, unsigned char inNumMaxAtomNeighbors);
	virtual ~GrainCandidateList();
	void init(RecursiveGrainIdentificationEngine *inOwner,unsigned char inNumMaxAtomNeighbros);
	void add(GrainCandidate * candidate);
	void testAndAddToGrain();
	//void addCandidates(GrainCandidate * inTestCandidates, unsigned char nNeighbors);
	void addCandidates(GrainCandidate * parent, GrainCandidate * inTestCandidates, unsigned char nNeighbors);
	long getNumCandidates();
	GrainCandidate * getCandidate(long iC);
	void buildNextGeneration(GrainCandidateList& nextGenList, double rSqrMin, double rSqrMax);
private:
	//Orientator * orient;
	RecursiveGrainIdentificationEngine * owner;
	std::vector<GrainCandidate> list;
	unsigned char nMaxAtomNeighbors;
};

class RecursiveGrainIdentificationEngine {
	public:
	RecursiveGrainIdentificationEngine(Orientator * inOrient, unsigned char inNumMaxAtomNeighbors, double inAngleThreshold);
	virtual ~RecursiveGrainIdentificationEngine();
	void setSearchRadiiSquared(double inRsqrMin, double inRsqrMax);
	void init(Orientator * orient, unsigned char inNumMaxAtomNeighbors, double inAngleThreshold);
	void setup(gID inGrainId, Grain * inGrain, double angularThreshold);
	long start(Atom * parent, AtomBox * curBox, long curAtomNum, double * inRelAtomPos);
	void addToGrain(GrainCandidate * candidate);
	bool belongsToGrain(GrainCandidate * candidate);
	bool strictBelongsToGrain(GrainCandidate * candidate);
	void test(GrainCandidateList * candidates);
	gID getGrainId();
	long getNumberOfAtoms();
	//Start atom for infection identified by curAtomNum and curBox
	//returns the number of found atoms for that start atom
	private:
	void cleanupCandidatesMem();
	GrainCandidateList * curCandidates = nullptr;
	GrainCandidateList * nextGenCandidates = nullptr;
	double rSqrMin, rSqrMax;
	double cosHalfThreshold, bigCosHalfThreshold;
	unsigned char nMaxAtomNeighbors;
	Grain * grain;
	gID grainId;
	Orientator * orient;
};

#endif /* GRAINIDENTIFICATOR_H_ */

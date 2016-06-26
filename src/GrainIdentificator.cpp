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

#include "GrainIdentificator.h"

GrainIdentificator::GrainIdentificator(Orientator * inOrient, AtomBox * inBoxes, long inNumBoxes, double angleThreshold, unsigned char inNumMaxAtomNeighbors) {
	numGrains = 0;
	grainCapacity = 0;
	grainAlloc = GRAINALLOC;
	boxes = inBoxes;
	numBoxes = inNumBoxes;
	nMaxAtomNeighbors = inNumMaxAtomNeighbors;
	newEmptyGrain();
	orient = inOrient;
	engine = new RecursiveGrainIdentificationEngine(orient, nMaxAtomNeighbors, angleThreshold);
}

GrainIdentificator::~GrainIdentificator() {
	long grainSize = grains.size();
	for (long i = 0; i < grainSize; i++) {
		deleteLastGrain();
	}
	delete engine;
}

long GrainIdentificator::run(double angularThreshold) {
	std::cout << "Thread " << omp_get_thread_num() << ": Grain Identification started." << std::endl;
	//for all atoms run the recursive engine
	AtomBox * box;
	long iA;
	long numAssignedAtoms;
	double atomPos[3];
	MeanOrientation curMeanOri;
	for (long iB = 0; iB < numBoxes; iB++) {
		box = boxes + iB;
		for (iA = 0; iA < box->getNumAtoms(); iA++) {
			box->obtainGlobalAtomPos(iA, atomPos);
			engine->setup(grains.size() - 1, grains.back(), angularThreshold);
			numAssignedAtoms = engine->start(box->getAtom(iA), box, iA, atomPos);
			if (numAssignedAtoms > 200) {
				std::cout << "Thread " << omp_get_thread_num() << ": Found grain " << grains.size() - 1 << " with " << numAssignedAtoms << " atoms" << std::endl;
				newEmptyGrain();
				continue;
			}
			deleteLastGrain();
			newEmptyGrain();
		}
	}
	deleteLastGrain();
	calculateOrientationSpread();
#pragma omp critical
{
	std::cout << LINE << "\n"
	<<"Thread " << omp_get_thread_num() <<": Grain Identification Done (2/3)\n"
	<< LINE <<std::endl;
}
	return grains.size();
}

bool sortGrainsByVolume(Grain* g1, Grain* g2) {
	return g1->getNumberOfAtoms() > g2->getNumberOfAtoms();
}

void GrainIdentificator::sort() {
	std::sort(grains.begin(), grains.end(), sortGrainsByVolume);
	//assign the new defined grain ids to all atoms
	assign();
}

void GrainIdentificator::assign() {
	for (long iG = 0; iG < grains.size(); iG++) {
		//dont save the too small ones
		if (grains[iG]->getNumberOfAtoms() > 1) {
			grains[iG]->assign(iG);
		}
	}
}

void GrainIdentificator::setSearchRadiiSquared(double inRsqrMin, double inRsqrMax) {
	engine->setSearchRadiiSquared(inRsqrMin, inRsqrMax);
}

Grain * GrainIdentificator::getGrain(gID grainID) {
	return grains.at(grainID);
}

long GrainIdentificator::getNumGrains() {
	return grains.size();
}

void GrainIdentificator::newEmptyGrain() {
	grains.push_back(new Grain());
	numGrains++;
}

bool sortByGrainIds(Atom * a1, Atom *a2) {
	return a1->getGrainId() < a2->getGrainId();
};

void GrainIdentificator::assignOrphanAtoms(long depth) {
	std::cout <<"Thread " << omp_get_thread_num() <<": Starting Orphan Atom Adoption." << std::endl;
	std::vector<AtomID> unassignedAtoms;
	//test all atoms and check whether they are assigned to any grain or not
	Atom * atom;
	AtomBox * box;
	AtomID atomId;
	for (atomId.iB = 0; atomId.iB < numBoxes; atomId.iB++) {
		box = boxes + atomId.iB;
		for (atomId.iA = 0; atomId.iA < box->getNumAtoms(); atomId.iA++) {
			atom = box->getAtom(atomId.iA);
			if (atom->getGrainId() == NO_GRAIN) {
				unassignedAtoms.push_back(atomId);
			}
		}
	}
#ifdef DEBUGMODE
	std::cout << "All " << unassignedAtoms.size() << " unassigned atoms found " << std::endl;
#endif
	//now try to assign the most frequent grain id of the neighbors
	long numUnAssignedAtoms = unassignedAtoms.size();
	std::vector<bool> isUnassigned;
	isUnassigned.reserve(numUnAssignedAtoms);
	for (long iA = 0; iA < numUnAssignedAtoms; iA++) {
		isUnassigned.push_back(true);
	}

	long oldNumUnAssignedAtoms;
	long curDepth = 0;
	if (depth == 0) {
		do {
			curDepth++;
			std::cout <<"Thread " << omp_get_thread_num() <<": Orphan Atom Adoption Iteration " << curDepth << " with " << numUnAssignedAtoms << " unassigned atoms" << std::endl;
			oldNumUnAssignedAtoms = numUnAssignedAtoms;
			orphanAtomAssignAttempt(unassignedAtoms, isUnassigned, numUnAssignedAtoms);
			if (unassignedAtoms.size() > 2 * numUnAssignedAtoms && unassignedAtoms.size() > 500) {
				rebuildOrphanAtomList(unassignedAtoms, isUnassigned, numUnAssignedAtoms);
			}
		} while (numUnAssignedAtoms < oldNumUnAssignedAtoms);
	} else {
		for (curDepth = 1; curDepth <= depth; curDepth++) {
			std::cout <<"Thread " << omp_get_thread_num() <<": Orphan Atom Adoption Iteration " << curDepth << " with " << numUnAssignedAtoms << " unassigned atoms" << std::endl;
			oldNumUnAssignedAtoms = numUnAssignedAtoms;
			orphanAtomAssignAttempt(unassignedAtoms, isUnassigned, numUnAssignedAtoms);
			if (unassignedAtoms.size() > 2 * numUnAssignedAtoms && unassignedAtoms.size() > 500) {
				rebuildOrphanAtomList(unassignedAtoms, isUnassigned, numUnAssignedAtoms);
			}
		}
	}
#pragma omp critical
{
	std::cout << LINE << "\n"
	<<"Thread " << omp_get_thread_num() <<": Orphan Atom Adoption Done (3/3)\n"
	<< LINE <<std::endl;
}
}

void GrainIdentificator::rebuildOrphanAtomList(std::vector<AtomID> & unassignedAtoms, std::vector<bool> &isUnassigned, long numUnassignedAtoms) {
	std::vector<AtomID> shortListUnassignedAtoms;
	std::vector<bool> shortIsUnassigned;
	shortListUnassignedAtoms.reserve(numUnassignedAtoms);
	shortIsUnassigned.reserve(numUnassignedAtoms);
	for (long iA = 0; iA < unassignedAtoms.size(); iA++) {
		if (isUnassigned[iA]) {
			shortListUnassignedAtoms.push_back(unassignedAtoms[iA]);
			shortIsUnassigned.push_back(true);
		}
	}
	unassignedAtoms = shortListUnassignedAtoms;
	isUnassigned = shortIsUnassigned;
}

//!\brief Tries to assign the given unassigned (='orphan') atoms to an already existing grain
void GrainIdentificator::orphanAtomAssignAttempt(std::vector<AtomID> & unassignedAtoms, std::vector<bool> &isUnassigned, long &numUnAssignedAtoms) {
	AtomBox * box;
	Atom * atom;
	Atom * nborAtom;
	unsigned char nNeighbors;
	Atom * neighbors[nMaxAtomNeighbors];
	Occurrence maxGrainOcc;
	for (long iA = 0; iA < unassignedAtoms.size(); iA++) {
		//don't check already assigned atoms
		if (!isUnassigned[iA]) {
			continue;
		}
		//create references to the current box and atom
		box = boxes + unassignedAtoms[iA].iB;
		atom = box->getAtom(unassignedAtoms[iA].iA);
		//calculate nearest neighbors
		nNeighbors = box->nearestAtomNeighbors(unassignedAtoms[iA].iA, nMaxAtomNeighbors, neighbors);
		if (nNeighbors <= 0) {
			continue;
		}
		//find the maximal occurrence of grain ids.
		maxGrainOcc = sortAndFindMaxGrainOccurrence(neighbors, nNeighbors);
		//end maximal occurrence found
		if (maxGrainOcc.count < 4) {
			continue;
		}
		if (maxGrainOcc.id == NO_GRAIN) {
			continue;
		}
		//now assign atom to grain with grain id.
		grains[maxGrainOcc.id]->addOrphan(atom);
		atom->setGrainId(maxGrainOcc.id);
		isUnassigned[iA] = false;
		numUnAssignedAtoms--;
	}

}
//!\brief Sorts a given atom-set by ascending grain-id and finds the grain, which occurs most frequently.
//!\param[in,out] atomSet pointer to list of atom-pointers
//!\param[in] setSize size of the atom set (number of entries).
//!\return Occurrence object, containing the most frequent grain-id and corresponding number of occurrence (count).
Occurrence GrainIdentificator::sortAndFindMaxGrainOccurrence(Atom ** atomSet, unsigned char setSize) const{
	Occurrence curOcc;
	Occurrence maxOcc;
	maxOcc.id = NO_GRAIN;
	maxOcc.count = 0;
	curOcc.id = atomSet[0]->getGrainId();
	curOcc.count = 0;
	Atom * atom;
	gID curGrainId;
	std::sort(atomSet, atomSet + setSize, sortByGrainIds);
	for (unsigned char iN = 0; iN < setSize; iN++) {
		atom = atomSet[iN];
		curGrainId = atom->getGrainId();
		if (curGrainId == NO_GRAIN) {
			continue;
		}
		if (curGrainId == curOcc.id) {
			curOcc.count++;
		} else {
			if (curOcc.count > maxOcc.count) {
				maxOcc = curOcc;
			}
			curOcc.id = curGrainId;
			curOcc.count = 1;
		}
	}
	if (curOcc.count > maxOcc.count) {
		maxOcc = curOcc;
	}
	return maxOcc;
}

void GrainIdentificator::deleteLastGrain() {
	delete grains.back();
	grains.pop_back();
	numGrains--;
}

RecursiveGrainIdentificationEngine::RecursiveGrainIdentificationEngine(Orientator *inOrient, unsigned char inNumMaxAtomNeighbors, double inAngleThreshold) {
	init(inOrient, inNumMaxAtomNeighbors, inAngleThreshold);
}

RecursiveGrainIdentificationEngine::~RecursiveGrainIdentificationEngine() {
	cleanupCandidatesMem();
}

void RecursiveGrainIdentificationEngine::init(Orientator *inOrient, unsigned char inNumMaxAtomNeighbors, double inAngleThreshold) {
	orient = inOrient;
	nMaxAtomNeighbors = inNumMaxAtomNeighbors;
	cosHalfThreshold = ori::cosHalfFromRad(inAngleThreshold);
}

void RecursiveGrainIdentificationEngine::setSearchRadiiSquared(double inRsqrMin, double inRsqrMax) {
	rSqrMin = inRsqrMin;
	rSqrMax = inRsqrMax;
}

void RecursiveGrainIdentificationEngine::setup(gID inGrainId, Grain *inGrain, double angularThreshold) {
	grainId = inGrainId;
	grain = inGrain;
	cosHalfThreshold = ori::cosHalfFromRad(angularThreshold);
	bigCosHalfThreshold = ori::cosHalfFromRad(3 * angularThreshold);
}

long RecursiveGrainIdentificationEngine::start(Atom * parent, AtomBox * curBox, long curAtomNum, double * inRelAtomPos) {
	//initialization for the while loop (or while the for loop?) :-)
	cleanupCandidatesMem();
	curCandidates = new GrainCandidateList(this, nMaxAtomNeighbors);
	nextGenCandidates = new GrainCandidateList(this, nMaxAtomNeighbors);
	//setup an active graincandidate object for the parent atom
	GrainCandidate parentCandidate;
	parentCandidate.init(curBox, curAtomNum, inRelAtomPos);
	parentCandidate.setParent(&parentCandidate);
	//add the parent atom to the candidates list
	if (belongsToGrain(&parentCandidate)) {
		curCandidates->add(&parentCandidate);
	}
	//end initialization
	//cout << "Started iteration for grain " << grainId << endl;
	while (curCandidates->getNumCandidates() > 0) {
		curCandidates->testAndAddToGrain();
		curCandidates->buildNextGeneration(*nextGenCandidates, rSqrMin, rSqrMax);
		delete curCandidates; //for the next iteration we dont need the past candidates
		curCandidates = nextGenCandidates; // for the next iteration, next generation becomes adult
		nextGenCandidates = new GrainCandidateList(this, nMaxAtomNeighbors); //birth of a next generation
	}
	if (grain->getNumberOfAtoms() > 0) {
		grain->calcAverageCenter();
		grain->recalculateMeanOrientation();
	}
	cleanupCandidatesMem();
	return (grain->getNumberOfAtoms());
}

void RecursiveGrainIdentificationEngine::addToGrain(GrainCandidate *candidate) {
	Atom * atom;
	double * pos;
	atom = candidate->getAtom();
	pos = candidate->getPosition();
	//
	atom->setGrainId(grainId);
	grain->add(atom, orient);
	grain->addToCenter(pos);
	//recalc the orientation each 100s atom
	if (grain->getNumberOfAtoms() % 100 == 1) {
		grain->recalculateMeanOrientation();
	}
}

bool RecursiveGrainIdentificationEngine::belongsToGrain(GrainCandidate * candidate) {
	Atom * atom = candidate->getAtom();
	Atom * parentAtom;
	//3 exit criteria:
	//1. atom is already assigned to any grain
	//2. atom has no defined orientation (is a GB-atom)
	//3. atom has a "different orientation" than the parent atom
	if (atom->getGrainId() != NO_GRAIN) {
		return false;
	}
	if (atom->getOrientationId() == NO_ORIENTATION) {
		return false;
	}
	parentAtom = candidate->getParentAtom();
	if (!orient->haveCloseOrientations(parentAtom, atom, cosHalfThreshold)) {
		return false;
	}
	return true;
}

bool RecursiveGrainIdentificationEngine::strictBelongsToGrain(GrainCandidate * candidate) {
	Atom * atom = candidate->getAtom();
	Atom * parentAtom;
	//4 exit criteria:
	//1. atom is already assigned to any grain
	//2. atom has no defined orientation (is a GB-atom)
	//3. atom has a "different orientation" than the parent atom
	//4. atom has a huge different orientation from the grains mean orientation
	if (atom->getGrainId() != NO_GRAIN) {
		return false;
	}
	if (atom->getOrientationId() == NO_ORIENTATION) {
		return false;
	}
	parentAtom = candidate->getParentAtom();
	if (!orient->haveCloseOrientations(parentAtom, atom, cosHalfThreshold)) {
		return false;
	}
	if (!ori::haveCloseOrientations(grain->getOrientation()->getQuaternion(), orient->getOrientation(atom->getOrientationId())->getQuaternion(),
			bigCosHalfThreshold)) {
		return false;
	}
	return true;
}

void RecursiveGrainIdentificationEngine::test(GrainCandidateList *candidates) {
	//3 exit criteria:
	//1. atom is already assigned to any grain
	//2. atom has no defined orientation (is a GB-atom)
	//3. atom has a "different orientation" than the parent atom
	Atom * atom;
	Atom * parentAtom;
	GrainCandidate * candidate;
	for (long i = 0; i < candidates->getNumCandidates(); i++) {
		candidate = candidates->getCandidate(i);
		atom = candidate->getAtom();
		if (atom->getGrainId() != NO_GRAIN) {
			continue;
		}
		if (atom->getOrientationId() == NO_ORIENTATION) {
			continue;
		}
		parentAtom = candidate->getParentAtom();
		if (!orient->haveCloseOrientations(parentAtom, atom, cosHalfThreshold)) {
			continue;
		}

		candidate->activate();
	}
}

gID RecursiveGrainIdentificationEngine::getGrainId() {
	return grainId;
}

long RecursiveGrainIdentificationEngine::getNumberOfAtoms() {
	return grain->getNumberOfAtoms();
}

void RecursiveGrainIdentificationEngine::cleanupCandidatesMem() {
	if (curCandidates != nullptr) {
		delete curCandidates;
		curCandidates = nullptr;
	}
	if (nextGenCandidates != nullptr) {
		delete nextGenCandidates;
		nextGenCandidates = nullptr;
	}
}

GrainCandidate::GrainCandidate() {
	active = false;
}

GrainCandidate::GrainCandidate(AtomBox *inBox, long inAtomId, double *inPos) {
	init(inBox, inAtomId, inPos);
}

void GrainCandidate::init(AtomBox *inBox, long inAtomId, double * inPos) {
	box = inBox;
	atomId = inAtomId;
	position[0] = inPos[0];
	position[1] = inPos[1];
	position[2] = inPos[2];
	active = false;
}

void GrainCandidate::activate() {
	active = true;
}

unsigned char GrainCandidate::findNeighbors(unsigned char inNumMaxAtomNeighbors, GrainCandidate *outNeighborsList, double rSqrMin, double rSqrMax) {
	AtomBoxP * nborBoxesList = new AtomBoxP[inNumMaxAtomNeighbors];
	long * nborAtomIdList = new long[inNumMaxAtomNeighbors];
	double * neighborVects = new double[DIM * inNumMaxAtomNeighbors];
	unsigned char nNeighbors = box->atomNeighbors(atomId, rSqrMin, rSqrMax, inNumMaxAtomNeighbors, nborBoxesList, nborAtomIdList, neighborVects);
	double * curVect;
	for (unsigned char iNA = 0; iNA < nNeighbors; iNA++) {
		curVect = neighborVects + iNA * DIM;
		if ((SQR(curVect[0]) + SQR(curVect[1]) + SQR(curVect[2])) > rSqrMax) {
			std::cout << "ERROR: " << "it is too big " << std::endl;
		}
		curVect[0] += position[0]; //add the position of this candidate to the new candidates
		curVect[1] += position[1];
		curVect[2] += position[2];
		outNeighborsList[iNA].init(nborBoxesList[iNA], nborAtomIdList[iNA], curVect);
	}
	delete[] nborBoxesList;
	delete[] nborAtomIdList;
	delete[] neighborVects;
	return nNeighbors;
}

GrainCandidate::~GrainCandidate() {
}

GrainCandidateList::GrainCandidateList(RecursiveGrainIdentificationEngine * inOwner, unsigned char inNumMaxAtomNeighbors) {
	init(inOwner, inNumMaxAtomNeighbors);
}

void GrainCandidateList::init(RecursiveGrainIdentificationEngine * inOwner, unsigned char inNumMaxAtomNeighbors) {
	owner = inOwner;
	nMaxAtomNeighbors = inNumMaxAtomNeighbors;
}

void GrainCandidateList::add(GrainCandidate * candidate) {
	list.push_back(*candidate);
}

void GrainCandidateList::addCandidates(GrainCandidate * parent, GrainCandidate * inTestCandidates, unsigned char nNeighbors) {
	GrainCandidate * curTestCandidate;
	for (unsigned char iA = 0; iA < nNeighbors; iA++) {
		curTestCandidate = inTestCandidates + iA;
		curTestCandidate->setParent(parent);
		add(curTestCandidate);
	}
}

GrainCandidateList::~GrainCandidateList() {
	list.clear();
}

void GrainCandidateList::testAndAddToGrain() {
	GrainCandidate * firstCandidate = &(list[0]);
	GrainCandidate * curCandidate;
	for (long i = 0; i < list.size(); i++) {
		curCandidate = firstCandidate + i;
		if (owner->getNumberOfAtoms() <= 100) {
			if (owner->belongsToGrain(curCandidate)) {
				owner->addToGrain(curCandidate);
				curCandidate->activate();
			}
			continue;
		}
		//restrict the later atoms
		if (owner->strictBelongsToGrain(curCandidate)) {
			owner->addToGrain(curCandidate);
			curCandidate->activate();
		}
	}
}

void GrainCandidateList::buildNextGeneration(GrainCandidateList & nextGenList, double rSqrMin, double rSqrMax) {
	//for all current candidates build the neighbor list
	//and save them in the nextGenList, if they meet all criteria
	GrainCandidate * testCandidates = new GrainCandidate[nMaxAtomNeighbors];
	unsigned char nNeighbors;
	//loop over all candidates
	for (long i = 0; i < list.size(); i++) {
		if (list[i].isActive()) {
			//
			nNeighbors = list[i].findNeighbors(nMaxAtomNeighbors, testCandidates, rSqrMin, rSqrMax);
			if (nNeighbors > nMaxAtomNeighbors) {
				std::cerr << "ERROR too much neighbors" << std::endl;
			}
			nextGenList.addCandidates(&(list[i]), testCandidates, nNeighbors);
		}
	}
	delete[] testCandidates;
}

long GrainCandidateList::getNumCandidates() {
	return list.size();
}

GrainCandidate * GrainCandidateList::getCandidate(long iC) {
	return &(list[iC]);
}

Atom * GrainCandidate::getAtom() {
	return box->getAtom(atomId);
}

AtomBox * GrainCandidate::getBox() {
	return box;
}

bool GrainCandidate::isActive() {
	return active;
}

double *GrainCandidate::getPosition() {
	return position;
}

void GrainCandidate::setParent(GrainCandidate *inParent) {
	parent = inParent->getAtom();
}

Atom * GrainCandidate::getParentAtom() {
	return parent;
}

long GrainCandidate::getAtomId() {
	return atomId;
}

void GrainIdentificator::calculateOrientationSpread() {
	for(long iG = 0; iG < numGrains; iG++){
		grains[iG]->calcTotalOrientationSpread(orient);
	}
}

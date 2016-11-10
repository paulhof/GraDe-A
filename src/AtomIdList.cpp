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

#include "AtomIdList.h"

AtomIdList::AtomIdList() {
}

AtomIdList::~AtomIdList() {
}

void AtomIdList::add(AtomID& id) {
	if (id.iB >= inverseList.size()){
		//list is too small -> reserve
		reserveInverse(id.iB+1);
	}
	if (id.iA >= inverseList[id.iB]->size()) {
		//leaf is too small -> resize
		inverseList[id.iB]->resize(id.iA+1);
	}
	//save corresponding list index
	(*inverseList[id.iB])[id.iA] = list.size();
	//and save the id
	list.push_back(id);
}

void AtomIdList::reserve(long nAtoms, long nBoxes) {
	if (nAtoms >= 0) {
		list.reserve(nAtoms);
		reserveInverse(nBoxes);
	}
}

const AtomID* AtomIdList::getAtomId(long iAtom) const{
	if (iAtom >= 0 && iAtom < list.size()) {
		return list.data() + iAtom;
	}
	return nullptr;
}

void AtomIdList::reserveInverse(long nBoxes) {
	if (inverseList.size() < nBoxes) {
		inverseList.reserve(nBoxes);
		for (long i = inverseList.size() - 1; i < nBoxes; i++ ){
			inverseList.push_back(new std::vector<long>);
		}
	}
}

long AtomIdList::getAtomNum(const AtomID& id) const {
	if (id.iB < inverseList.size()){
		if(id.iA < inverseList[id.iB]->size()){
			return (*inverseList[id.iB])[id.iA];
		}
	}
	return -1;
}

void AtomIdList::deleteInverseList() {
	for (long i = 0; i < inverseList.size(); i++) {
		delete inverseList[i];
	}
}

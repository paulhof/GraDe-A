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

#ifndef ATOMIDLIST_H_
#define ATOMIDLIST_H_
#include "GradeA_Defs.h"

class AtomIdList {
public:
	AtomIdList();
	virtual ~AtomIdList();
	void reserve(long nAtoms, long nBoxes);
	void add(AtomID & id);
	const AtomID * getAtomId(long iAtom) const;
	long getAtomNum(const AtomID & id ) const;
private:
	void reserveInverse(long nBoxes);
	void deleteInverseList();
	std::vector<AtomID> list;
	std::vector<std::vector<long>* > inverseList;
};


#endif /* ATOMIDLIST_H_ */

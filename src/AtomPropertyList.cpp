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

#include "GradeA_Defs.h"
#include "AtomPropertyList.h"

AtomPropertyList::AtomPropertyList() {

}

AtomPropertyList::~AtomPropertyList() {
	int i;
	for(i = 0; i < integerProperties.size(); i++){
		delete integerProperties[i];
	}
	for(i = 0; i < floatProperties.size(); i++){
		delete floatProperties[i];
	}
}

int AtomPropertyList::addProperty(const std::string & name, long reservedSize, bool isInteger) {
	propertyNames.push_back(name);
	propertyIsInt.push_back(isInteger);
	if (isInteger) {
		integerProperties.push_back(new std::vector<int>);
		integerProperties.back()->reserve(reservedSize);
		propertyId.push_back(integerProperties.size() - 1);
	} else {
		floatProperties.push_back(new std::vector<double>);
		floatProperties.back()->reserve(reservedSize);
		propertyId.push_back(floatProperties.size() - 1);
	}
	return getNumProperties() - 1;
}

void AtomPropertyList::addIntPropertyValue(int propertyNum, int value) {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
		if  (propertyIsInt[propertyNum]) {
				//property is int
				int intId = propertyId[propertyNum];
				integerProperties[intId]->push_back(value);
		}
	}
}

void AtomPropertyList::addFloatPropertyValue(int propertyNum, double value) {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
		if  (!propertyIsInt[propertyNum]) {
			//property is float
			int floatId = propertyId[propertyNum];
			floatProperties[floatId]->push_back(value);
		}
	}
}

int AtomPropertyList::getNumProperties() const {
	return propertyNames.size();
}

std::string AtomPropertyList::getPropertyName(int propertyNum) const {
	if (propertyNum >= 0 && propertyNum < propertyNames.size()) {
		return propertyNames[propertyNum];
	}
	return "";
}

int AtomPropertyList::getIntPropertyValue(int propertyNum, long atomNum) const {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
	//property num is valid
		if  (propertyIsInt[propertyNum]) {
		//property is int
		int intId = propertyId[propertyNum];
			if( atomNum < integerProperties[intId]->size() ) {
				//a value for the atom exists
				return (*integerProperties[intId])[atomNum];
			}
		}
	}
	//otherwise return an arbitrary value
	return 0;
}

void AtomPropertyList::convertPropertyToFloat(int propertyNum) {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
		//property num is valid
		if  (propertyIsInt[propertyNum]) {
			//property is int -> convert
			int intId = propertyId[propertyNum];
			//change the bool flag to indicate a
			propertyIsInt[propertyNum] = false;
			//create a new float property
			floatProperties.push_back(new std::vector<double>);
			//save the corresponding id
			propertyId[propertyNum] = floatProperties.size() - 1;
			//
			floatProperties.back()->resize(integerProperties[intId]->size());
			//copy and convert all values from integer to float
			for (int i = 0; i < integerProperties[intId]->size(); i ++) {
				(*floatProperties.back())[i] =
						static_cast<double>((*integerProperties[intId])[i]);
			}
			//delete the old property
			delete integerProperties[intId];
			integerProperties[intId] = nullptr;
		}
	}
}

double AtomPropertyList::getFloatPropertyValue(int propertyNum,	long atomNum) const {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
		//property num is valid
		if  (!propertyIsInt[propertyNum]) {
		//property is float
		int floatId = propertyId[propertyNum];
			if( atomNum < floatProperties[floatId]->size() ) {
				//a value for the atom exists
				return (*floatProperties[floatId])[atomNum];
			}
		}
	}
	//otherwise return an arbitrary value
	return 0.0;
}

bool AtomPropertyList::isPropertyInt(int propertyNum) const {
	if ( propertyNum >= 0 && propertyNum < getNumProperties()) {
		return propertyIsInt[propertyNum];
	}
	return true;
}

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

#ifndef ATOMPROPERTYLIST_H_
#define ATOMPROPERTYLIST_H_

class AtomPropertyList {
public:
	AtomPropertyList();
	virtual ~AtomPropertyList();
	int addProperty(const std::string & name, long reservedSize , bool isInteger = true);
	void convertPropertyToFloat(int propertyNum);
	void addIntPropertyValue(int propertyNum, int value);
	void addFloatPropertyValue(int propertyNum, double value);
	std::string getPropertyName(int propertyNum) const;
	int getNumProperties() const;
	int getIntPropertyValue(int propertyNum, long atomNum) const;
	double getFloatPropertyValue(int propertyNum, long atomNum ) const;
	bool isPropertyInt(int propertyNum) const;
private:
	std::vector <std::string> propertyNames;
	std::vector <bool> propertyIsInt;
	std::vector <int> propertyId;
	//
	std::vector <std::vector<int> * > integerProperties;
	std::vector <std::vector<double> * > floatProperties;
};

#endif /* ATOMPROPERTYLIST_H_ */

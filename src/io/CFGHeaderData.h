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

#ifndef IO_CFGHEADERDATA_H_
#define IO_CFGHEADERDATA_H_
#include "FileEditor.h"
#include "FileImporter.h"
class CFGHeaderData {
public:
	CFGHeaderData();
	~CFGHeaderData(){};
	void parse(FileEditor * editor);
	void parse(TextReader * reader);
	long getNumLines() const {
		return numLines;
	}
	long getNumParticles() const {
		return numParticles;
	}
	int getNumAuxFields() const;
	int getEntryCount() const {
		return entryCount;
	}
	std::string getAuxField(int fieldNum) const;
	bool containsVelFields() const;
	bool isExtendedFormat() const {
		return extendedFormat;
	}
	const Matrix3 * getH0() const;
	const Matrix3 * getTransform() const;
private:
	void parseLine(std::string & line, bool & isDone);
	long numLines = 0;
	long numParticles = 0;
	double unitMultiplier = 1.;
	Matrix3 H0;
	Matrix3 transform;
	double rateScale = 1.;
	bool extendedFormat = false;
	bool containsVelocities = false;
	std::vector<std::string> auxFields;
	int entryCount = 0;
};


#endif /* IO_CFGHEADERDATA_H_ */

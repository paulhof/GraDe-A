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
//This file has been taken from OVITO.

#ifndef CFG_FILE_IMPORTER_H
#define CFG_FILE_IMPORTER_H

#include "FileImporter.h"
#include "CFGHeaderData.h"
/**
 * \brief File parser for AtomEye CFG files.
 */
class CFGImporter : public FileImporter
{
public:

	/// \brief Constructs a new instance of this class.
	CFGImporter(std::string &filename, AtomContainer * inData) : FileImporter(filename, inData, cfg){}
	/// \brief Checks if the given file has format that can be read by this importer.
	virtual bool checkFileFormat();
	virtual void parseFile();
	void readAtomPositions(int particleIndex, const char* s);
	double parseField(int particleIndex, int columnIndex, const char* token, const char* token_end);
private:
	void skipExtendedLines();
	void readAtom();
	std::vector<std::string> curAtomData;
	long curAtomNum = -1;
	Matrix3 transform;
	double translate[3];
	double * atomPositions = nullptr;
	CFGHeaderData header;
};

#endif // CFG_FILE_IMPORTER_H

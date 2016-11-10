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

#ifndef IMPORT_CFGEDITOR_H_
#define IMPORT_CFGEDITOR_H_
#include "../GradeA_Defs.h"
#include "FileEditor.h"
#include "FileImporter.h"
#include "CFGHeaderData.h"

class CFGEditor {
public:
	CFGEditor(std::string inFileName);
	virtual ~CFGEditor();
	void readNextAtom();
	void setCurAux(int auxNum, std::string value);
	std::string getCurAux(int auxNum) const;
	std::string getAuxName(int auxNum) const;
	int getNumAuxFields() const;
	int getCurNumAuxFields() const;
	long getNumAtoms() const;
	long getCurAtomNum() const;
	void flush();
	bool eof() const;
	void close();
private:
	void parseHeader();
	void skipExtendedLines();
	std::string curFieldString();
	int auxFieldNum (int auxNum) const;
	std::vector<std::string> curAtomData;
	long curAtomNum = -1;
	CFGHeaderData header;
	FileEditor editor;
	const char numPosFields = 3;
	char numVelFields = 0;
	bool headerParsed = false;
};



#endif /* IMPORT_CFGEDITOR_H_ */

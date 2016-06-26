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

#ifndef CSVTABLEREADER_H_
#define CSVTABLEREADER_H_
#include <fstream>

#include "../GradeA_Defs.h"
class CSVTableReader {
public:
	CSVTableReader();
	CSVTableReader(std::string inFilename, int inNumHeaderLines = 0, std::string inSeparator = ";");
	virtual ~CSVTableReader();
	bool parse();
	void setFileName(std::string inFilename);
	void setCustomSeparator(std::string inSeparator);
	int getNumHeaderLines() const;
	long getNumLines() const;
	long getNumColumns() const;
	std::string getColumnName(int colNum) const;
	std::string entry(long lineNum, int colNum) const;
	std::string headerEntry(int lineNum, int colNum) const;
	std::string headerLine(int lineNum) const;
private:
	void parseHeaderLine(std::string & inHeaderLine);
	void parseTitleLine(std::string& inTitleLine);
	void parseLine(std::string& inLine);
	//!\brief Adds a column to the object.
	//!\param[in] The name string of the column to add.
	bool addNewColumn(std::string columnName);
	bool addEntryToHeaderLine(std::string entry);
	bool addEntryToCurrentLine(std::string entry, std::string inColumnName);
	bool addEntryToCurrentLine(std::string entry, long columnId);
	//!\brief Adds a header line to the object.
	//!\param[in] line The line text to add.
	void addHeaderLine(std::string line);
	void addNewLine();
	bool openFile();
	bool closeFile();
	void safeAddToLine(std::string entry, long columnId);
	std::vector<std::vector<std::string >> headerLines;
	std::vector<std::string> columnNames;
	std::vector<std::vector<std::string >> lines;
	std::string separator = ";";
	std::string fileName;
	std::ifstream filestream;
	int numHeaderLines; //! number of headerlines to parse before data
};

#endif /* CSVTABLEREADER_H_ */

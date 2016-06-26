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

#ifndef CSVTABLEWRITER_H_
#define CSVTABLEWRITER_H_
#include <fstream>

#include "../GradeA_Defs.h"
/*! \brief Class which writes table data into a csv-type file.
 The class is capable of using custom separator strings and additional header lines.
*/
class CSVTableWriter {
public:
	CSVTableWriter();
	CSVTableWriter(std::string inFilename, std::string inSeparator = ";");
	//!\brief Method which writes the csv table into the file identified by fileName.
	void write();
	//!\brief Adds a header line to the object.
	//!\param[in] line The line text to add.
	void addHeaderLine(std::string line);
	//!\brief Adds a column to the object.
	//!\param[in] columnName The name string of the column to add.
	bool addNewColumn(std::string columnName);
	bool addEntryToCurrentLine(std::string entry, std::string inColumnName);
	bool addEntryToCurrentLine(std::string entry, long columnId);
	void addNewLine();
	void setFileName(std::string inFilename);
	std::string getSeparator();
	void setCustomSeparator(std::string inSeparator);
	void disableTitleLine();
	virtual ~CSVTableWriter();
private:
	void flush();
	void writeLine(std::string inLine);
	void writeText(std::string inText);
	bool openFile();
	bool closeFile();
	void safeAddToLine(std::string entry, long columnId);
	std::vector<std::string> headerLines;
	std::vector<std::string> columnNames;
	bool titleLineActive = true;
	std::vector<std::vector<std::string >> lines;
	std::string separator;
	std::string fileName;
	std::ofstream filestream;
};

#endif /* CSVTABLEWRITER_H_ */

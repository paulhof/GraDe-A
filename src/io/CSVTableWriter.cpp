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

#include "CSVTableWriter.h"

CSVTableWriter::CSVTableWriter() {
	fileName = "";
	separator = ";";
}

CSVTableWriter::CSVTableWriter(std::string inFilename, std::string inSeparator)
{
	setFileName(inFilename);
	setCustomSeparator(inSeparator);
}

CSVTableWriter::~CSVTableWriter() {

}

void CSVTableWriter::write()
{
	if(!openFile()){
		std::cerr << "Cannot open file \"" << fileName << "\" for writing." << std::endl;
		return;
	}

	//header lines
	for(int hNum = 0; hNum < headerLines.size(); hNum++){
		writeLine(headerLines[hNum]);
	}
	std::string lineString;
	if(titleLineActive){
		//title line
		for(int colNum = 0; colNum < columnNames.size(); colNum++){
			lineString += columnNames[colNum] + separator;
		}
		//delete last separator
		lineString.erase(lineString.end()-separator.length());
		writeLine(lineString);
	}
	//entries
	for(long lineNum = 0; lineNum < lines.size(); lineNum++){
		lineString = "";
		for(int colNum = 0; colNum < columnNames.size(); colNum++){
			if(colNum < lines[lineNum].size()){
				lineString += lines[lineNum][colNum] + separator;
			} else {
				//added because variable number of columns can cause problems with several parsers
				lineString += separator;
			}
		}
		//delete last separator
		lineString.erase(lineString.end()-separator.length());

		if(lineNum < lines.size()-1){
			writeLine(lineString);
		} else {
			//last line without \n
			writeText(lineString);
		}
	}
	flush();
	if(!closeFile()){
		std::cerr << "Cannot close file \"" << fileName << "\" after writing." << std::endl;
	}
}

bool CSVTableWriter::addNewColumn(std::string columnName)
{
	columnNames.push_back(columnName);
	return true;
}

bool CSVTableWriter::addEntryToCurrentLine(std::string entry, std::string inColumnName)
{
	for (int i = 0; i < columnNames.size(); i++){
		if(inColumnName == columnNames[i] ){
			safeAddToLine(entry,i);
			return true;
		}
	}
	return false;
}

bool CSVTableWriter::addEntryToCurrentLine(std::string entry, long  columnId)
{
	if(columnId >= columnNames.size()){
		return false;
	}
	safeAddToLine(entry,columnId);
	return true;
}

void CSVTableWriter::addNewLine()
{
	lines.resize(lines.size()+1);

}

void CSVTableWriter::setFileName(std::string inFilename)
{
	fileName = inFilename;
}

void CSVTableWriter::setCustomSeparator(std::string inSeparator)
{
	separator = inSeparator;
}



bool CSVTableWriter::openFile()
{
	try{
		filestream.open(fileName.c_str());
	} catch (...){
		return false;
	}
	return filestream.is_open();
}

void CSVTableWriter::writeLine(std::string inLine)
{
	filestream << inLine << "\n";
}

void CSVTableWriter::writeText(std::string inText) {
	filestream << inText;
}

bool CSVTableWriter::closeFile()
{
	try{
		filestream.close();
	} catch (...){
		return false;
	}
	return !filestream.is_open();
}



void CSVTableWriter::addHeaderLine(std::string line)
{
	headerLines.push_back(line);
}

void CSVTableWriter::flush() {
	filestream.flush();
}

void CSVTableWriter::safeAddToLine(std::string entry, long  columnId)
{
	if(lines.size() <= 0) {
		lines.resize(1);
	}
	if(lines.back().size() <= columnId){
		lines.back().resize(columnId+1, "");
	}
	lines.back()[columnId] = entry;
}

std::string CSVTableWriter::getSeparator() {
	return separator;
}

void CSVTableWriter::disableTitleLine() {
	titleLineActive = false;
}

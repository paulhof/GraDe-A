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

#include "../io/FileImporter.h"
FileImporter::FileImporter(std::string &inFilename, AtomContainer * inData, FileType inFileType){
	filename = inFilename;
	data = inData;
	fileType = inFileType;
}

FileImporter::~FileImporter() {
}

const char* TextReader::readLine()
{
	lineNumber++;
	getline();
	if( eof() ) throw EndOfFile();
	return line_str.c_str();
}

bool TextReader::eof()
{
 return fileStream->eof();
}

long TextReader::getLineNumber()
{
	return lineNumber;
}

TextReader::TextReader(std::string & inFileName)
{
	filename = inFileName;
	lineNumber = 0;
	fileStream = new std::ifstream;
	fileStream->open(filename.c_str());
	if(fileStream->fail())
			throw Exception("TextReader-FileStream failed to open \"" + filename + "\"");
}

TextReader::~TextReader()
{
	if(!fileStream->is_open()) std::cerr << "TextReader-FileStream has no file opened - will fail when closing" << std::endl;
	if(fileStream->fail())
			throw Exception("TextReader-FileStream failed before closing \"" + filename +"\"");
	fileStream->close();
	if(fileStream->fail())
		throw Exception("TextReader-FileStream failed to close \"" + filename +"\"");
	delete fileStream;
}

std::string TextReader::lineString()
{
	return line_str;
}

int TextReader::getline()
{
	std::getline(*fileStream, line_str);
	return line_str.length()+1;
}








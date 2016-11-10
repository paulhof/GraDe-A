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

#ifndef ORIENTATORFILEQUEUE_H_
#define ORIENTATORFILEQUEUE_H_
#include "GradeA_Defs.h"
#ifdef _MSC_VER
	#include <windows.h>
#else
	#include <dirent.h>
#endif

typedef struct {
	std::string prefix;
	std::string number;
	std::string postfix;
	bool isNumber;
}FileNameID;

class OrientatorFileQueue {
public:
	OrientatorFileQueue();
	virtual ~OrientatorFileQueue();
	void softFindFiles();
	void autoFindFiles();
	void sortFilesByNumber();
	std::string nextFileName();
	std::string curFileName() const;
	std::string previousFileName() const;
	std::string fileName(int fileNum);
	int numFiles() const;
	void init(std::string inSubDirectoryName, std::string inFileNamePreFix, std::string inFileNamePostFix);
	void initByWildcard(std::string inWildCard, int inStartFileNum = 0, int inEndFileNum = INT_MAX);
	void setFileFixByFileWildcard(std::string inFileWildCard);
	void setDirectory(std::string inDirName);
	std::string getFileNamePreFix() const;
	std::string getCurFileNameId() const;
	int getCurFileNameNum() const;
	int getFileNameNum(int fileNum) const;
	std::string getFileNameId(int fileNum) const;
	std::string getFileNamePostFix() const;
	std::string getFileNameEnding() const;
	std::string getDirName() const;
	//
	std::string outCsvFileName(int fileNum) const;
	std::string outCsvFileNameWildCard() const;
	std::string outCfgFileName(int fileNum) const;
	std::string outOrientationCsvFileName(int fileNum) const;
private:
	bool addFile(std::string fileName, bool softMode = false);
	FileNameID fileId(std::string fileNameIdString) const;
	void decomposePostFix(std::string inPostFix);
	void adjustSubDirectoryName();
	std::vector<std::string> filesInDirectory(std::string dir);
	std::string preSeparator = "_";
	std::string fileNamePreFix = "";
	std::string fileNamePostFix = "";
	std::string fileNameEnding = "";
	std::vector<FileNameID> fileNumbers;
	std::vector<std::string> fileNameIds;
	std::string subDirectoryName;
	int curFileNum;
	//
	int startFileNum = 0;
	int endFileNum = INT_MAX;
};

#endif /* ORIENTATORFILEQUEUE_H_ */

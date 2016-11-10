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

#include "OrientatorFileQueue.h"

OrientatorFileQueue::OrientatorFileQueue() {
	curFileNum = 0;
	fileNamePreFix = "";
	fileNamePostFix = "";
}

void OrientatorFileQueue::init(std::string inSubDirectoryName, std::string inFileNamePreFix, std::string inFileNamePostFix)
{
	subDirectoryName = inSubDirectoryName;
	adjustSubDirectoryName();
	fileNamePreFix = inFileNamePreFix;
	decomposePostFix(inFileNamePostFix);
}

void OrientatorFileQueue::initByWildcard(std::string inWildCard, int inStartFileNum, int inEndFileNum)
{
	if (inWildCard.find_first_of("'\"") == 0) {
		inWildCard = inWildCard.substr(1);
	}
	if (inWildCard.find_last_of("'\"") == inWildCard.length()-1) {
		inWildCard = inWildCard.substr(0,inWildCard.length()-1);
	}
	startFileNum = inStartFileNum;
	endFileNum = inEndFileNum;
	size_t dirPos;
	std::string dirString = {DIRCHAR};
	dirPos = inWildCard.rfind(dirString);
	std::string fileNameWildCard;
	if(dirPos == std::string::npos) {
		subDirectoryName = ".";
		subDirectoryName += DIRCHAR;
		fileNameWildCard = inWildCard;
	} else {
		subDirectoryName = inWildCard.substr(0,dirPos+1);
		fileNameWildCard = inWildCard.substr(dirPos+1,std::string::npos);
	}
	setFileFixByFileWildcard(fileNameWildCard);
}

void OrientatorFileQueue::setFileFixByFileWildcard(std::string inFileWildCard)
{
	size_t idPos;
	idPos = inFileWildCard.rfind('*');
	if(idPos == std::string::npos){
		fileNamePreFix = "";
		fileNamePostFix = inFileWildCard;
		decomposePostFix(fileNamePostFix);
		return;
	}
	fileNamePreFix = inFileWildCard.substr(0,idPos);
	fileNamePostFix = inFileWildCard.substr(idPos+1, std::string::npos);
	decomposePostFix(fileNamePostFix);
	if (fileNamePostFix.length() > 0){
		if (fileNamePostFix.back() == '_' ){
				preSeparator = "";
		}
	}
}

void OrientatorFileQueue::setDirectory(std::string inDirName){
	subDirectoryName = inDirName;
	adjustSubDirectoryName();
}

void OrientatorFileQueue::decomposePostFix(std::string inPostFix) {
	size_t dotPos = inPostFix.rfind('.');
	fileNamePostFix = inPostFix.substr(0,dotPos);
	fileNameEnding = inPostFix.substr(dotPos, std::string::npos);
}

void OrientatorFileQueue::autoFindFiles()
{
	std::vector<std::string> fileNameList = filesInDirectory(subDirectoryName);
	for(int fileNum = 0; fileNum < fileNameList.size(); fileNum++){
		if(	addFile(fileNameList[fileNum])){
			std::cout << "Found file " << fileNameList[fileNum] << std::endl;
		}
	}
	//now sort the filenames by number
	sortFilesByNumber();
}

void OrientatorFileQueue::softFindFiles() {
	std::vector<std::string> fileNameList = filesInDirectory(subDirectoryName);
	for(int fileNum = 0; fileNum < fileNameList.size(); fileNum++){
		if(	addFile(fileNameList[fileNum],true)){
			std::cout << "Found file " << fileNameList[fileNum] << std::endl;
		}
	}
	//now sort the filenames by number
	sortFilesByNumber();
}

bool sortByNumber(FileNameID id1, FileNameID id2){
	//put the numbers first and order them in ascending order
	//bug resolved: keep in mind:
	//http://en.cppreference.com/w/cpp/concept/Compare
	//if id1 ==  id2 and both isNumber is false return true
	//will lead to undefined behavior
	if((!id1.isNumber) && (!id2.isNumber)) return false;
	if((!id1.isNumber) && (id2.isNumber)) return false;
	if((id1.isNumber) && (!id2.isNumber)) return true;
	return atoi(id1.number.c_str()) < atoi(id2.number.c_str());
}

void OrientatorFileQueue::sortFilesByNumber(){
	std::sort(fileNumbers.begin(),fileNumbers.end(),sortByNumber);
	for(int i = 0; i < fileNameIds.size(); i++){
		fileNameIds[i] = fileNumbers[i].prefix +  fileNumbers[i].number +  fileNumbers[i].postfix;
	}
}

void OrientatorFileQueue::adjustSubDirectoryName(){
	if(subDirectoryName.back() != DIRCHAR){
		subDirectoryName += DIRCHAR;
	}
}

OrientatorFileQueue::~OrientatorFileQueue() {
}

std::string OrientatorFileQueue::nextFileName()
{
	return fileName(++curFileNum);
}

std::string OrientatorFileQueue::curFileName() const {
	return subDirectoryName + fileNamePreFix + fileNameIds[curFileNum] + fileNamePostFix + fileNameEnding;
}

std::string OrientatorFileQueue::previousFileName() const{
	if (curFileNum == 0 ){
		return "";
	}
	return subDirectoryName + fileNamePreFix + fileNameIds[curFileNum-1] + fileNamePostFix + fileNameEnding;
}

int OrientatorFileQueue::numFiles() const {
	return fileNameIds.size();
}

std::string OrientatorFileQueue::getFileNameEnding() const{
	return fileNameEnding;
}

bool OrientatorFileQueue::addFile(std::string fileName, bool softMode) {
	size_t pos;
	//fileName must begin with prefix
	pos = fileName.find(fileNamePreFix);
	if(pos != 0){
		return false;
	}
	//fileName must end with postfix and ending
	std::string totalPostFix = fileNamePostFix+fileNameEnding;
	pos = fileName.rfind(totalPostFix);
	if(pos != fileName.length() - totalPostFix.length()){
		return false;
	}
	//now we can extract the part in between prefix and postfix
	std::string curFileIdString = fileName.substr(fileNamePreFix.length(),pos-fileNamePreFix.length());
	FileNameID curFileId = fileId(curFileIdString);

	//allow only files with number greater or equal to startFileNum
	if((!softMode) && curFileId.isNumber){
		//if wildcard has no postfix, allow only files id-string beginning with number,
		//because else output cfg names will be found in subsequent runs
		if(fileNamePostFix.length() == 0){
			if(curFileId.prefix.length() != 0){
				return false;
			}
		}
		int curFileIdNum = atoi(curFileId.number.c_str());
		if( curFileIdNum < startFileNum || curFileIdNum > endFileNum){
			return false;
		}
	}
	//Otherwise add file to list
	fileNameIds.push_back(curFileIdString);
	fileNumbers.push_back(curFileId);
	return true;
}

FileNameID OrientatorFileQueue::fileId(std::string fileNameIdString) const{
	std::string numString = "0123456789";
	FileNameID curFileId;
	size_t startpos, endpos;
	//find the first number in the string
	startpos = fileNameIdString.find_first_of(numString);
	if(startpos == std::string::npos){
		//no number has been found
		curFileId.prefix = fileNameIdString;
		curFileId.number = "";
		curFileId.postfix = "";
		curFileId.isNumber = false;
		return curFileId;
	}
	endpos = fileNameIdString.find_first_not_of(numString,startpos);
	if(endpos != std::string::npos){
		//Id contains postfix

		curFileId.prefix = fileNameIdString.substr(0,startpos);
		curFileId.number = fileNameIdString.substr(startpos,endpos-startpos);
		curFileId.postfix = fileNameIdString.substr(endpos,fileNameIdString.length()-endpos);
		curFileId.isNumber = true;
		return curFileId;
	}
	//Id contains prefix and number
	curFileId.prefix = fileNameIdString.substr(0,startpos);
	curFileId.number = fileNameIdString.substr(startpos,fileNameIdString.length()-startpos);
	curFileId.postfix = "";
	curFileId.isNumber = true;
	return curFileId;
}

std::vector<std::string> OrientatorFileQueue::filesInDirectory(std::string dirName)
{
	std::vector<std::string> fileNameList;
#ifdef _MSC_VER
	//Windows API needs asterisk at the end of dirName
	if(dirName.size() > 0) {
		if (dirName.back() == DIRCHAR){
			dirName.push_back('*');
		}
	}
	HANDLE findHandle;
	WIN32_FIND_DATA FindFileData;
	if((findHandle = FindFirstFile(dirName.c_str(), &FindFileData)) != INVALID_HANDLE_VALUE){
	    do{
	    	fileNameList.push_back(FindFileData.cFileName);
	    } while(FindNextFile(findHandle, &FindFileData));
	    FindClose(findHandle);
	}
#else
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (dirName.c_str())) != nullptr) {
	  while ((ent = readdir (dir)) != nullptr) {
		  fileNameList.push_back(ent->d_name);
	 }
	  closedir (dir);
	}

#endif
	return fileNameList;
}

std::string OrientatorFileQueue::fileName(int fileNum){
	if (fileNum < 0 ) return "";
	if (fileNum >= fileNameIds.size()) return "";
	curFileNum = fileNum;
	return subDirectoryName + fileNamePreFix + fileNameIds[fileNum] + fileNamePostFix + fileNameEnding;
}

std::string OrientatorFileQueue::getFileNamePreFix() const {
	return fileNamePreFix;
}

std::string OrientatorFileQueue::getFileNamePostFix() const {
	return fileNamePostFix;
}

std::string OrientatorFileQueue::getCurFileNameId() const {
	return getFileNameId(curFileNum);
}

std::string OrientatorFileQueue::getFileNameId(int fileNum) const{
	if (fileNum < 0 ) {
		return "";
	}
	if (fileNum >= fileNameIds.size()){
		return "";
	}
	return fileNameIds[fileNum];
}

int OrientatorFileQueue::getFileNameNum(int fileNum) const {
	if (fileNum < 0 ) {
		return -1;
	}
	if (fileNum >= fileNameIds.size()){
		return -1;
	}
	return atoi(fileNumbers[fileNum].number.c_str());
}

std::string OrientatorFileQueue::getDirName() const{
	return subDirectoryName;
}

int OrientatorFileQueue::getCurFileNameNum() const {
	return atoi(fileNumbers[curFileNum].number.c_str());
}

std::string OrientatorFileQueue::outCsvFileName(int fileNum) const{
	std::string postSeparator = "_";
	if (getFileNameId(fileNum).front() == '_' ){
				postSeparator = "";
	}
	return getFileNamePreFix() + getFileNamePostFix() + preSeparator + "GrainData"
				+ postSeparator + getFileNameId(fileNum) + ".csv";
}

std::string OrientatorFileQueue::outCfgFileName(int fileNum)  const{
	std::string postSeparator = "_";
	if (getFileNameId(fileNum).front() == '_' ){
		postSeparator = "";
	}
	return getFileNamePreFix() + getFileNamePostFix() + preSeparator + "AtomData" +
				postSeparator +getFileNameId(fileNum) + ".cfg";
}

std::string OrientatorFileQueue::outCsvFileNameWildCard() const {
	return getFileNamePreFix() + getFileNamePostFix() + preSeparator + "GrainData_*.csv";
}

std::string OrientatorFileQueue::outOrientationCsvFileName(int fileNum) const {
	std::string postSeparator = "_";
	if (getFileNameId(fileNum).front() == '_' ){
				postSeparator = "";
	}
	return getFileNamePreFix() + getFileNamePostFix() + preSeparator + "OriData"
					+ postSeparator + getFileNameId(fileNum) + ".csv";
}

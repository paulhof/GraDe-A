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

#ifndef GRAINTIMEEVOLUTIONWRITER_H_
#define GRAINTIMEEVOLUTIONWRITER_H_
#include "../GradeA_Defs.h"
#include "CSVTableWriter.h"
#include "CSVTableReader.h"
#include "GrainCSVFileFormat.h"
#include "../OrientatorFileQueue.h"

#define NOT_INSIDE -1

class GrainTimeEvolutionWriter {
public:
	GrainTimeEvolutionWriter(std::string inCsvFileNameWildCard);
	virtual ~GrainTimeEvolutionWriter();
	//! runs the grain time evolution writing....
	void run();
private:
	void write();
	void addCurrentTimeStepDataToWriters();
	void initPropertiesFromFirstFile();
	void initWritersForEachProperty();
	int grainIdPos(std::string grainId);
	void addNewGrainId(std::string grainId);
	void resetReader();
	void resetWriters();
	std::string inputFileWildCard;
	std::vector<std::string> propertyNames;
	std::vector<std::string> grainIdOrder;
	//
	CSVTableWriter * writers = nullptr;
	int numWriters = 0;
	CSVTableReader * reader = nullptr;
	OrientatorFileQueue csvQueue;
};

#endif /* GRAINTIMEEVOLUTIONWRITER_H_ */

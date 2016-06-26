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

#ifndef IMPORT_GRAINCSVFILEFORMAT_H_
#define IMPORT_GRAINCSVFILEFORMAT_H_
#include "../GradeA_Defs.h"
#include "../AtomContainer.h"
#include "CSVTableWriter.h"
#include "../ContainerData.h"
#include "../GradeA_Version.h"
#include "CSVTableReader.h"

namespace GrainCSVFile {
	enum class Column {
			GrainID,
			NumAtoms , NumRegularAtoms , NumOrphanAtoms,
			PosX, PosY , PosZ,
			phi_1 , PHI, phi_2 ,
			OrientationSpread,
			quaternion_0, quaternion_1, quaternion_2, quaternion_3
	};
	enum class OptionalColumn {
		Volume,
		misOrientation,
		cubMisOrientation,
		travelledDistance
	};
	class Format {
	public:
		Format();
		virtual ~Format();
		bool isRightFormat(CSVTableReader * csvTable) const;
		bool init(CSVTableReader * csvTable, ContainerData * data) const;
		void prepareTable (CSVTableWriter * writer) const;
		void fillTable(CSVTableWriter * csvTable, const ContainerData * container, const CubicLattice & material) const;
		void fillTable(CSVTableWriter * csvTable, const AtomContainer * container, const CubicLattice & material) const;
		int getNumColumns() const;

		const std::string getColName(int colNum) const {
			if(colNum < 0 || colNum >= getNumColumns()) {
				return "";
			}
			return colNames[colNum];
		}

		const int getContainerSizeColNum() const {
			return containerSizeColNum;
		}

		const int getContainerSizeLineNum() const {
			return containerSizeLineNum;
		}

		const int getNumHeaderLines() const {
			return numHeaderLines;
		}

		const int getPeriodicColNum() const {
			return periodicColNum;
		}

		const std::string& getPeriodicFlag() const {
			return periodicFlag;
		}

		const int getPeriodicLineNum() const {
			return periodicLineNum;
		}
		const int getColNum(Column col) const;
	private:
		//Header line setup
		const int numHeaderLines = 2;
		const int containerSizeLineNum = 1;
		const int periodicLineNum = 1;
		const int containerOriginLineNum = 1;
		const int containerSizeColNum = 0;	//size 0-2
		const int periodicColNum = 3;		//periodic 3
		const int containerOriginColNum = 4;//origin 4-6
		const std::string periodicFlag = "p";
		//Column setup
		const std::vector<std::string> colNames {
			"Grain ID",
			"NumAtoms", "NumRegularAtoms", "NumOrphanAtoms",
			"PosX","PosY","PosZ",
			"phi_1", "PHI", "phi_2",
			"OriSpread",
			"q0","q1","q2","q3"
		};
		const std::vector<std::string> optionalColNames {
			"Volume",
			"misOri","cubMisOri",
			"travelDistance"
		};
	};
}

GrainCSVFile::Format const csvFormat;
#endif /* IMPORT_GRAINCSVFILEFORMAT_H_ */

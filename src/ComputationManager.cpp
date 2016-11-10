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

#include "ComputationManager.h"
#include "CubicLattices.h"

#if !defined(WINDOWS) || defined(CYGWIN)
#include <sys/stat.h>
#define mkdir(x) mkdir(x,0777);
#else
#include <direct.h>
#define mkdir(x) _mkdir(x)
#endif

bool isNumeric(const std::string& input) {
    return true;
}

ComputationManager::ComputationManager(bool inPeriodic, double latticeParameter, double inAngularThreshold, std::string chemElementName, bool inPrintOrientations) {
	material = new FccLattice(latticeParameter, VOLUMEUNIT,chemElementName);
	periodic = inPeriodic;
	printOrientations = inPrintOrientations;
	NN_searchRadiusSqrMin = SQR(0.9 * HALFSQRT2 *  latticeParameter);
	NN_searchRadiusSqrMax = SQR(1.1 * HALFSQRT2 * latticeParameter);
	grainAngularThreshold = inAngularThreshold;
	boxSize = 1.1 * sqrt(NN_searchRadiusSqrMax);
}

void ComputationManager::run(std::string fileNameWildCard, std::string inInitGrainFileName, int startFileNum, int endFileNum) {
#ifdef _MSC_VER
	std::cout << "This is the MSVC version" << std::endl;
#endif
	initGrainFileName = inInitGrainFileName;
	mkdir(TIMEEVOSUBDIR);
	queue.initByWildcard(fileNameWildCard, startFileNum, endFileNum);
	queue.autoFindFiles();
	std::cout <<"Found " << queue.numFiles() << " files to compute." << std::endl;
	parallelRun();
	std::cout << LINE << "\n"
	<<"-- Finished whole Calculation --\n"
	<< LINE <<std::endl;
}

void ComputationManager::parallelRun() {
	OrientatorFileQueue privateQueue = queue;
	//Measure computation time
	StopWatch parallelWatch;
	parallelWatch.trigger();
	//Even though this is a thread-parallelization, the threads do not share any data.
	//The amount of memory increases linearly with number of threads.
#pragma omp parallel shared(std::cout) firstprivate(privateQueue) default(none)
{
#pragma omp single
{
	std::cout << "Running computation in parallel with " << omp_get_num_threads() << " threads" << std::endl;
}
	//DO NOT WRITE TO ANY MEMBER OF THIS CLASS OBJECT INSIDE THE PARALLEL REGION (shared memory) !!!!
	//ESPECIALLY DONT USE the queue object, use privateQueue instead!

	//Therefor construct a manager object for each thread
	ComputationManager threadManager (periodic, material->getLatticeParameter(), grainAngularThreshold, material->getName(), printOrientations);
#pragma omp for
	for(int iF = 0; iF < privateQueue.numFiles(); iF++){
#pragma omp critical
	{
		std::cout << "Thread " << omp_get_thread_num() << ": Running file " << iF ;
		std::cout << " with filename \"" << privateQueue.fileName(iF) << "\" " << std::endl;
	}
		threadManager.runFile(iF, privateQueue);
	}
}
	parallelWatch.trigger();
	std::cout << LINE << "\n"
		<<"-- Main Computation (Parallel part) Done --\n"
		<< LINE <<std::endl;
	StopWatch serialWatch;
	serialWatch.trigger();

	//Run a grain tracking in order to be able to plot e.g. mass over time
	GrainTracker tracker(&queue,material,initGrainFileName);
	tracker.run();

	//write time evolution files
	std::cout << "Running TimeEvo-Writer with \"" << queue.outCsvFileNameWildCard() << "\"" << std::endl;
	GrainTimeEvolutionWriter timeEvo(queue.outCsvFileNameWildCard());
	timeEvo.run();

	serialWatch.trigger();
	std::cout << "\nParallel time = " << parallelWatch.getString() << std::endl;
	std::cout << "Serial time = " << serialWatch.getString() << std::endl;
	std::cout << "Serial overhead = " << serialWatch.getDuration()/(serialWatch.getDuration()+parallelWatch.getDuration())*100.0 << " %" <<  std::endl;
}

void ComputationManager::runFile(int fileNum, OrientatorFileQueue & inQueue) {
	queue = inQueue;
	if (fileNum >= queue.numFiles() ){
		return;
	}
	runSingleFile(fileNum);
}

void ComputationManager::runSingleFile(int fileNum) {
	//Init a new container object in order to store atom position data
	initContainer();
	std::string inputFileName;
	inputFileName = queue.fileName(fileNum);

	//import-object utilized to read data into the container
	CFGImporter import(inputFileName, container);

	//check whether file has cfg format, exit if not
	bool isCfg;

	try {
		isCfg = import.checkFileFormat();
	} catch (...) {
		std::cerr << "File \"" << queue.curFileName()
				<< "\" skipped - could not be read." << std::endl;
		delete container;
		container = nullptr;
		return;
	}

	if (!isCfg) {
		std::cerr << "File \"" << queue.curFileName()
				<< "\" skipped - has wrong format." << std::endl;
		delete container;
		container = nullptr;
		return;
	}

	//now since file has cfg format, parse it
	try {
		import.parseFile();
	} catch (...) {
		std::cerr << "File \"" << queue.curFileName()
				<< "\" skipped - could not be parsed." << std::endl;
		delete container;
		container = nullptr;
		return;
	}
	//---------------------------------------------------------------------
	//MAIN EXECUTION:

#pragma omp critical
{
	std::cout << "Thread " << omp_get_thread_num() << ": Orientation Calculation started. " << std::endl;
}
	//analyze the dataset
	//atom-wise orientation calculation
	container->calculateAtomOrientations(NN_searchRadiusSqrMin, NN_searchRadiusSqrMax);
#pragma omp critical
{
	std::cout << LINE << "\n"
	<<"Thread " << omp_get_thread_num() <<": Orientation Calculation Done (1/3)\n"
	<< LINE <<std::endl;
}
	//print orientations
	if( printOrientations){
		std::cout << "Writing orientation csv file: " << queue.outOrientationCsvFileName(fileNum) << " for " << container->getNumOrientations() << " orientations" << std::endl;
		OrientatorPrinter * oriPrinter = new OrientatorPrinter(container->getOrientator(), queue.outOrientationCsvFileName(fileNum));
		oriPrinter->print();
		delete oriPrinter;
	}
	//grain identification
#pragma omp critical
{
	std::cout << "Thread " << omp_get_thread_num() << ": Grain Identification started." << std::endl;
}
	container->identifyGrains(grainAngularThreshold, NN_searchRadiusSqrMin, NN_searchRadiusSqrMax);
#pragma omp critical
{
	std::cout << LINE << "\n"
	<<"Thread " << omp_get_thread_num() <<": Grain Identification Done (2/3)\n"
	<< LINE <<std::endl;
}
	//---------------------------------------------------------------------
	//OUTPUT:
	//write the atom data into a cfg file
	std::string outputCfgFileName = queue.outCfgFileName(fileNum);

	std::cout << "Writing cfg file: " << outputCfgFileName << std::endl;
	writeCfgFile(outputCfgFileName);

	//write the output grain data into a csv file
	std::string outputCsvFileName = queue.outCsvFileName(fileNum);
	std::cout << "Writing csv file: " << outputCsvFileName << std::endl;
	writeCsvTableFile(outputCsvFileName);

	delete container;
	container = nullptr;

	std::cout << "Finished calculation of file " << fileNum+1 << " with filename " << queue.curFileName() << std::endl << std::endl;
}

void ComputationManager::initContainer() {
	if (periodic) {
		container = new PeriodicAtomContainer(boxSize);
	} else {
		container = new AtomContainer(boxSize);
	}
}

ComputationManager::~ComputationManager() {
	delete material;
	material = nullptr;
}

void ComputationManager::writeCfgFile(std::string fileName) {
	AtomIO output;
	if(!output.openCfgFile(fileName)){
		//opening failed
		return;
	}
	std::vector<std::string> propertyNames;
	container->getAtomPropertyNames(propertyNames);
	std::vector<std::string> grainAvgPropertyNames;
	container->getGrainAvgPropertyNames(propertyNames);
	//append the average-propertynames to propertynames
	propertyNames.insert(propertyNames.end(),grainAvgPropertyNames.begin(),grainAvgPropertyNames.end());
	long nAtoms = container->getNumAtoms();
	output.writeCfgFileHeader(container->getSize(), nAtoms, propertyNames.size(), propertyNames.data(), material->getName());
	std::vector<std::string> atomsProperties;

	double pos[DIM];
	for(long iA = 0; iA < nAtoms; iA++){
		container->getAtomsPosition(iA, pos);
		container->getAtomsProperties(iA, atomsProperties);
		output.appendAtomToCfgFile(pos, atomsProperties);
	}
	output.closeCfgFile();
}

void ComputationManager::writeCsvTableFile(std::string fileName) {
	csvTable = new CSVTableWriter(fileName);
	csvFormat.fillTable(csvTable,container,*material);
	csvTable->write();
	delete csvTable;
	csvTable = nullptr;
}

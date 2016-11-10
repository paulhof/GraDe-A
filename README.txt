#GraDe-A README FILE
#GraDe-A is released under GPLv3, for license details see LICENSE.txt

GraDe-A: Grain Detection Algorithm
GraDe-A is a software tool, which is able to detect 3D-grain-entities from MD-simulation data.
The software currently supports fcc-strucutured materials.
The code has been developed by Paul Hoffrogge during a Masterthesis-project at IMM, RWTH-Aachen.

-------------------------------------------------------------------
Contact:
Dr.-Ing. Luis Antonio Barrales-Mora
@ Institut für Metallkunde und Metallphysik, RWTH-Aachen
barrales@imm.rwth-aachen.de
+49 (0) 241 80 28 298
-------------------------------------------------------------------

Install as described in INSTALL.txt.

Set number of threads with OMP_NUM_THREADS environment variable.
For Windows copy dependent .dll files in the same folder as the executable.

Run the grade-A executable with proper parameters.
E.g. for an Al-dataset of the form
inputfile_1.cfg
inputfile_2.cfg
...
inputfile_x.cfg
do:

Windows:
set OMP_NUM_THREADS=8
grade-A.exe "inputfile_*.cfg" p 4.05 1.0

Linux:
OMP_NUM_THREADS=8
./grade-A "inputfile_*.cfg" p 4.05 1.0

Please write the glob-pattern with "", the corresponding files are found by the software itself.
The program generates a folder "./TimeEvo/" and writes output-files for each inputfile.

Armadillo is open source software released under MPL2, for license details see either ./armadillo/LICENSE.txt or http://mozilla.org/MPL/2.0/
Eigen is open source software released under MPL2, for license details see http://mozilla.org/MPL/2.0/

Cfg-file-reading routines have been taken from Ovito and edited.
See http://ovito.org/ for full source.
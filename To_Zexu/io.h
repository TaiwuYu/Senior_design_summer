#ifndef IO_H
#define IO_H

void write_vtk(double *rfield, char filename[], char fieldname[] );
void write_vtk2(double *rfield, char filename[], char fieldname[] );
// to write vtk file

void OutputMed(int timestep);//to print the configuration in file
void init_from_file(double *tpr, char *filename);//to initialize the field "tpr" with "filenamae"
#endif

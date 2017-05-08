
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


/* Functions written by Cristoph Bach, 2010-07-27.
 * project:  Simulation and parameter estimation using ETAS and shakemaps
 */


#include "read_matrix.h"

#include <stdio.h>

#include "../defines.h"
#include "../util/files.h"

int read_matrix_N(int NL, char *infile, int columns, int headerlines, double **data, long *rows) {
/* Reads an ascii file into a double array.
 *
 * Input:
 *  N: no. of lines to read. set to -1 to read entire file.
 *  infile: input file name
 *  columns: no. of columns
 *  headerlines: no. of header lines, to be skipped.
 *
 * Output:
 *  data: contains the values from the file. Memory must be allocated previously.
 *  rows: no. of rows.
 */


	FILE *fin;
	double Zd, dum;
	int s, z, i, ans, dumerror;
	long ZZ, N;
	int c;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen("**Error: unable to open input file %s (read_matrix).**\n", infile);
		print_logfile("**Error: unable to open input file %s (read_matrix).**\n", infile);
		return (1);
	}

	for (i = 1; i <= headerlines; i++){
		//Read line of any length:
		c=1;
		while (c!='\n' && c!=EOF) c=fgetc(fin);
	}

	ans = 1;
	N = 0;

	while (ans != EOF) {
		ans = fscanf(fin, "%lf", &dum);
		if (ans==0){
			print_screen("**Error: invalid format at line %d of file %s. Only numeric characters allowed. (read_matrix).**\n", (int) ceil(N/(1.*columns)) + headerlines, infile);
			print_logfile("**Error: invalid format at line %d of file %s. Only numeric characters allowed. (read_matrix).**\n", (int) ceil(N/(1.*columns)) + headerlines, infile);
			return(1);
		}
		if (ans != EOF) N++;
	}

	fclose(fin);
	Zd = N * 1.0 / (1. * columns);
	if (rows) *rows = (long) Zd;	//if null pointer is passed, ignore.
	ZZ = (long) Zd;

	if (Zd - 1.0 * ZZ != 0.0) {
		print_screen("Error: Mismatch in the number of columns!");
		print_logfile("Error: Mismatch in the number of columns!");
		return(1);
	}

	if (NL>-1) ZZ=NL;

	fin = fopen(infile, "r");
	for (i = 1; i <= headerlines; i++){
		c=1;
		while (c!='\n' && c!=EOF) c=fgetc(fin);
	}
	for (z = 1; z <= ZZ; z++)
		for (s = 1; s <= columns; s++)
			dumerror = fscanf(fin, "%lf", &data[s][z]);
	fclose(fin);

	return(0);

}

int countline(char *filename){
/* Returns the number of lines in file "filename"
 */

	FILE *fin;
	int dum=0;
	int counter =0;

	if((fin = fopen(filename, "r"))==NULL) {
		print_screen(" **Error: unable to open input file %s. (countline.c)**\n", filename);
		print_logfile(" **Error: unable to open input file %s. (countline.c) **\n", filename);
		return -1;
	}

	while(dum!=EOF) {
	      dum=fgetc(fin);
	      if(dum=='\n')  counter++;
	    }

	fclose(fin);
	return counter;
}

int countcol_header(char *filename, int headerlines){
	/* Returns the number of columns in file "filename", using the first line after headerlines.
	 */

	char title[3000];
	FILE *fin;
	int dum=0;
	int counter =0;

	if((fin = fopen(filename, "r"))==NULL) {
		print_screen(" **Error: unable to open input file %s.**\n", filename);
		print_logfile(" **Error: unable to open input file %s.**\n", filename);
		return -1;
	}

	if (countline(filename)<=headerlines){
		print_screen(" **Error: headerlines should be less than total no. of lines in file %s.**\n", filename);
		print_logfile(" **Error: headerlines should be less than total no. of lines in file %s.**\n", filename);
		return -1;
	}

	for (int i = 1; i <= headerlines; i++){
			fgetline(fin, title, 3000);
	}
	dum=(int) '\t';
	while(dum=='\t' || dum==' ' || dum=='\n') dum=fgetc(fin);

	while (dum!='\n'){
		while(dum!='\t' && dum!=' ' && dum!='\n') {
			dum=fgetc(fin);
		}
		if(dum!='\n')  {
			counter++;
			while(dum=='\t' || dum==' ') dum=fgetc(fin);
			if (dum=='\n') counter-=1;
		}
	}
	fclose(fin);

	return counter+1;
}













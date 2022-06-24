
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


#include "read_csv.h"
#include <stdio.h>
#include <string.h>

#include "../defines.h"
#include "../util/files.h"

#define MAXCHAR 1000

int read_csv_N(int NL, char *infile, int columns, int headerlines, double **data, long *rows) {
/* Reads an ascii file into a double array.
 *
 * Input:
 *  N: no. of lines to read. set to -1 to read entire file.
 *  infile: input file name
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
 	char row[MAXCHAR];
	char *token;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen("**Error: unable to open input file %s (read_csv).**\n", infile);
		print_logfile("**Error: unable to open input file %s (read_csv).**\n", infile);
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
		ans = fscanf(fin, "%s", row);
		token = strtok(row,",");
                if (ans == EOF) break;
                while(token != NULL){
                    N++;
                    token = strtok(NULL, ",");
                }
	}

	fclose(fin);
	Zd = N * 1.0 / (1. * columns);
	if (rows) *rows = (long) Zd;	//if null pointer is passed, ignore.
	ZZ = (long) Zd;

	if (N % columns) {
		print_screen("Error: Mismatch in the number of columns. N=%d, columns=%d\n", N, columns);
		print_logfile("Error: Mismatch in the number of columns. N=%d, columns=%d\n", N, columns);
		return(1);
	}
	print_screen("File file contains %d rows, %d columns\n [read_csv]", ZZ, columns);
	print_logfile("File file contains %d rows, %d columns\n [read_csv]", ZZ, columns);


	if (NL>-1) ZZ=NL;

	fin = fopen(infile, "r");
	for (i = 1; i <= headerlines; i++){
		c=1;
		while (c!='\n' && c!=EOF) c=fgetc(fin);
	}
	for (z = 1; z <= ZZ; z++){

		fscanf(fin, "%s", row);
	        token = strtok(row,",");
		for (s = 1; s <= columns; s++){
			dumerror =sscanf(token, "%lf", &data[s][z]);
			token = strtok(NULL, ",");
		}
	}
	fclose(fin);

	return(0);

}

int countcol_csv(char *infile) {
/* Count columns in CSV file.
 */

	FILE *fin;
	int columns;
 	char row[MAXCHAR];
	char *token;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen("**Error: unable to open input file %s (read_csv).**\n", infile);
		print_logfile("**Error: unable to open input file %s (read_csv).**\n", infile);
		return (1);
	}

	fscanf(fin, "%s,", row);
	token = strtok(row,",");
	columns = 0;
        while(token != NULL){
	    columns++;
            token = strtok(NULL, ",");
        }

	fclose(fin);
	return columns;
}

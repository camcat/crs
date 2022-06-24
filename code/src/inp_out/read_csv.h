
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

#ifndef READ_CSV_H
#define READ_CSV_H

#include <stdio.h>

#include "../defines.h"
#include "../util/files.h"
#include "read_matrix.h"


#define read_csv(...) read_csv_N(-1, __VA_ARGS__)

//-----------------------------------------------------------------------------
// INPUT:       infile     :  Name der einzulesenden Datei
//              columns    :  Anzahl der Spalten in dieser Datei
//              headerlines:    "     "  Titelzeilen
// OUTPUT:      data       :  Matrix (S,Z) mit S  Spalten Z  Zeilen
//              rows       :  Zeilenzahl
//-----------------------------------------------------------------------------
int read_csv_N(int N, char *infile, int columns, int headerlines, double **data, long *rows);
int countcol_csv(char *infile);

#endif // READ_CSV_H

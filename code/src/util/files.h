
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


/* Functions from code written by Cristoph Bach, 2010
 * project:  Simulation and parameter estimation using ETAS and shakemaps
 */

int fgetword(FILE *fp, char *word, int count);
/* Liest ein von Trennzeichen begrenztes Wort ein */

char *sgetint(char *, int *);
/* liest aus String s eine von Blanks begrenzte int-Zahl ein */

int fgetsequ(FILE *, char *, int);
/*  liest eine Sequenz bestimmter Laenge ein */

int fgetline(FILE *, char *, int);
/* liest eine ganze Zeile ein */

int isalnum_0(int);
/* prueft, ob Zeichen alphanumerisch. Null ist Begrenzer! */

int isalnumde(int);
/* prueft, ob Zeichen alphanumerisch im erweiterten Zeichensatz */

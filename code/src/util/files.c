
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "files.h"

/* fgetword: ein Wort der Laenge count aus FILE fp einlesen,
   wenn count == 0, das ganze Wort einlesen. */

int fgetword(FILE *fp, char *word, int count) 
{
	int c, i;
	char *w;

	w = word;
	i = 0;
	while (isalnumde(c = getc(fp))) { /* Codierung angeben!  */
	  *w++ = c;
	  ++i;
	  if (count != 0 && i == count)
	    break;
	}
	*w = '\0';				
	return c;
}

/* liest aus String s ab Stelle sp int-Zahl num ein, gibt Zeiger auf naechste
Zahl zurueck. */

char *sgetint(char *sp, int *num)
{
  char s[10] = "", *p = s;
  
  for ( ; isdigit(*sp); *p++ = *sp++)
    ;
  *num = atoi(s);
  return s[0] == '\0' ? NULL : ++sp;
}


/* fgetsequ: eine Sequenz der Laenge count aus FILE fp einlesen,
 Leerzeichen == 0. */

int fgetsequ(FILE *fp, char *word, int count) 
{
	int c, i;
	char *w;

	w = word;
	for (i = 1; i < count + 1; i++) {
	  if ((c = getc(fp)) == EOF) 
	    return EOF;
	  *w++ = c;
	}
	*w = '\0';
	return c;
}

/* Liest aus file fp eine Zeile der Laenge max in String line ein */

int fgetline(FILE *fp, char *line, int max)
{
  if (fgets(line, max, fp) == NULL)
    return 0;
  else
    return strlen(line);
}

/* Prueft, ob c alphanum oder 0 */

int isalnum_0(int c)
{
  return (c == '0') ? 0 : isalnumde(c);
}

/* Prueft, ob c erweitert alphanum */

int isalnumde(int c)
{
  return (c > 127) ? -1 : isalnum(c);
}


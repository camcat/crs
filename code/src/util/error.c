
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


#include "error.h"
#include "../util/moreutil.h"
#include "../defines.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef _CRS_MPI
	#include "mpi.h"
#endif


void error_quit_fun(const char *fun, int quit, const char * format, ...){
    char buffer[3000];
    va_list args;
    va_start (args, format);

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		 vsprintf (buffer,format, args);

		 //print error to log file:
		 if (flog){
			  if (extra_verbose) {
				  fprintf(flog,"[%s]: %s",fun,buffer);
			  }
			  else{
				  fprintf(flog,"%s",buffer);
			  }
			  fflush(flog);
		 }

		 //print error to screen:
		  if (extra_verbose) {
			  fprintf(stderr,"[%s]: %s",fun,buffer);
		  }
		  else{
			  fprintf(stderr,"%s",buffer);
		  }
		  fflush(stderr);
		  va_end (args);
	}

	if (quit) exit(EXIT_FAILURE);
}

void print_logfile_fun(const char *fun, const char * format, ...){
/*
 *  Print message to log file.
 *  If flag: extra_verbose=1, also print out name of the calling function.
 */

  char buffer[3000];
  int procId = 0;
  va_list args;
  va_start (args, format);

  	#ifdef _CRS_MPI
  		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  	#endif

  	if (!flog) return;

  	if(procId == 0) {
		  vsprintf (buffer,format, args);
		  if (extra_verbose) {
			  fprintf(flog,"[%s]: %s",fun,buffer);
		  }
		  else{
			  fprintf(flog,"%s",buffer);
		  }

		  fflush(flog);
  	}
  	va_end (args);
}

void print_screen_fun(const char *fun, const char * format, ...){
/*
 *  Print message to screen.
 *  If flag: extra_verbose=1, also print out name of the calling function.
 *  print to stdout or stderr depending on whether message contains the word "error".
 */
  char buffer[3000];
  int procId = 0;
  FILE *out_channel;	//stderr or stdout, depending on string content.
  va_list args;
  va_start (args, format);

  	#ifdef _CRS_MPI
  		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  	#endif

  	if(procId == 0) {
  		if (!quiet){
		  vsprintf (buffer,format, args);
		  out_channel = (search_string(buffer, "Error") || search_string(buffer, "error"))? stderr : stdout;
		  if (extra_verbose) {
			  fprintf(out_channel, "[%s]: %s",fun,buffer);
		  }
		  else{
			  fprintf(out_channel, "%s",buffer);
		  }
  		  fflush(out_channel);
  		}
  	}
  	va_end (args);
}

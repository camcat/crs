
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


//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>

#include <stdio.h>

#include "../defines.h"

// Define macros to print out forecast (no. of earthquakes, binned into magnitudes) and cmbmap (stress field, not binned into magnitudes).
#define csep_forecast(...) csep_forecast_general(__VA_ARGS__,1)
#define csep_cmbmap(...) csep_forecast_general(__VA_ARGS__,0)

void csep_forecast_general(char *filename, struct crust crst, double *rates, int original_resolution,  int use_mags);
void write_csep_forecast(char *filename, double *lats, double *lons, double *deps, double dlat, double dlon, double ddep,
		double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag);
void write_csep_forecast_points(char *filename, double *lats, double *lons, double *deps, 
                double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag);

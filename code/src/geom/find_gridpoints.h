
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


#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../util/moreutil.h"

#include "../util/util1.h"

int find_closest_point(double *ys, double *xs, double *depths, int N, double y, double x, double Depth);
int find_gridpoints(double *ys, double *xs, double *dAs, double *depths, int N, double y, double x, double SD, double Depth, double SDd, int cut_sd, int *ngridj, int **ngridpointj, double **weightsj, int inside, int d3);
int find_gridpoints_d(double *ys, double *xs, double *depths, int *already_selected, int Nsel0, int N, double y_eq, double x_eq, double Depth, double m, double dDCFS, int *ngridj, int **ngridpointj);
int find_gridpoints_exact(double *ys, double *xs, double *depths, double dx, double dy, double dz, int N, int Nselmax, double y, double x, double SD, double Depth, double SDd,
		int cut_sd, int *ngridj, int **ngridpointj, double **weightsj, int inside, int d3);
double exact_prob_1d(double r, double dr, double sd);
double exact_prob(double rx, double ry, double rz, double dx, double dy, double dz, double sdx, double sdy, double sdz, int d3);
int all_nearestneighbours(double *x, double *y, int N, int **pts, double **dist);
int all_2ndnearestneighbours(double *x, double *y, int N, int **pts, double **dist);

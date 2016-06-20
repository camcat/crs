
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


#ifndef MEM_MGMT_H_
#define MEM_MGMT_H_


#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "../defines.h"

#include "../util/util1.h"

void reduce_eqkfm_memory(struct eqkfm *eqkfm0, int NF);
void check_empty_eqkfm(struct eqkfm eqkfm0, double toll, int *is_str, int *is_dip, int *is_open);
void shift_cat(struct catalog *cat, int N);
void init_crst(struct crust *crst);
void init_cat1(struct catalog *cat, int Zsel);
struct set_of_models *set_of_models_array(long n1, long n2);
struct eqkfm *eqkfm_array(long n1, long n2);
struct pscmp *pscmp_array(long n1, long n2);
struct pscmp *pscmp_arrayinit(struct crust v0, long n1, long n2);
void freefull_pscmparray(struct pscmp *v, long n1, long n2);
void free_cat(struct catalog cat);

#endif /* MEM_MGMT_H_ */


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


#ifndef EQKFM_COPY_H_
#define EQKFM_COPY_H_

#include "../defines.h"

void empty_eqkfm(struct eqkfm *eqkfm0);
void copy_eqkfm_nolocation_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_attributes(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_focmec(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_all(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_noslipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);

#endif /* EQKFM_COPY_H_ */

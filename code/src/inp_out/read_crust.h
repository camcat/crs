
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


#ifndef READ_CRUST_H_
#define READ_CRUST_H_

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../okada/prestress.h"
#include "../util/error.h"

#include "../util/util1.h"
#include "read_csep_template.h"
#include "read_matrix.h"

int read_crust(char *fnametemplate, char *focmecgridfile, struct crust *crst, double resxy, double resz, int multiple_focmecfiles);
int read_focmecgridfile(char *fname, struct crust *crst);

#endif /* READ_CRUST_H_ */


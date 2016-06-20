
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


#ifndef SETUP_EQKFM_H_
#define SETUP_EQKFM_H_


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"
#include "../inp_out/read_eqkfm.h"
#include "../inp_out/read_focmec.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../okada/okadaDCFS.h"
#include "../seis/soumod1.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "../util/splines_eqkfm.h"
#include "eqkfm_copy.h"
#include "find_timesteps.h"
#include "mem_mgmt.h"
#include "struct_conversions.h"

void set_current_slip_model(struct eqkfm *eqkfm0, int slipmodel_index);
int setup_catalogetc(char *catname, char **focmeccat, int nofmcat, struct tm reftime, double dDCFS, double Mag_source, double Mag_main, struct crust crst,
		struct catalog *cat, struct eqkfm **eqkfm1, double ***focmec, int **firstelements, struct flags flag, int *NFM, int *Ntot,
		double dt, double dM, double xytoll, double ztoll, double dR, double tw, double tstart, double tend);
int setup_aseismic_element(struct eqkfm *eqkfm0res, char **slipmodels, char *cmb_format, int no_snap,
						double mu, double disc, double tmain, double *tsnap, int nsel, int *sel_pts, int cuts_surf, double lat0, double lon0);
int setup_aseismic_eqkfm(struct slipmodels_list slipmodels, struct crust crst, struct eqkfm **eqkfm0res);
int setup_eqkfm_element(struct eqkfm *eqkfm0res, char **slipmodel, char *cmb_format, int no_slipmodels, double mu,
		double tmain, int nsel, int *sel_pts, double *mmain, int cuts_surf, int *NF0, double lat0, double lon0, int same_geom);
int load_newdata(double *t0, double t1, struct set_of_models *allmodels, int Nmain, int *NFaults, char **slipmodels, char **multimodels, int *no_slipmodels, int *Nmain_now);
int setup_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct Coeff_LinkList **Coefficients_aseismic, struct pscmp **DCFS_out, struct crust crst, struct eqkfm *eqkfm0,
		int Nm, int *Nfaults, struct eqkfm *eqkfm_aft, int no_afterslip, int *Nfaults_aft);
int update_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults);
#endif /* SETUP_EQKFM_H_ */

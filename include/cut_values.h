#ifndef CUT_VALUES_H_
#define CUT_VALUES_H_

namespace cutVals{

	const double Z_min = 0.3;
	const double Z_max = 1.0;
	const int bins_Z = 14;

	const double Q2_min = 2.0;
	const double Q2_max = 8.0;
	const int bins_Q2 = 12;

	const double xB_min = .1;
	const double xB_max = .66;
	const int bins_xB = 14;

	
	const double e_PCAL_W_min=19.0;
	const double e_PCAL_V_min=19.0;
	const double e_E_PCAL_min=0.07;
	const double PCAL_ECIN_SF_min=0.2;
	
	const double y_max=0.75;
	const double W_min=2.5;
	const double W_max=100;
	const double theta_min=5.; //For both e and pi
	const double theta_max=35.;
	const double P_pi_min=1.25;
	const double P_pi_max=5.0;
	const double P_e_min=3.0;
	const double P_e_max=10.6;

	const double Mx_min = 1.7;
	const double Mx_max = 5.0;
	
	const double Vz_e_min_inbending = -5.0;	
	const double Vz_e_max_inbending = 1.0;

	//Legacy cuts ("loose")
	const double SamplingFraction_min=0.17;
	const double Vz_e_min_inbending_loose = -13.0;	
	const double Vz_e_max_inbending_loose = 12.0;
	const double Vz_e_min_outbending_loose = -18.0;
	const double Vz_e_max_outbending_loose = 10.0;
	const double deltaVz_loose	=	20.0;
	const double Mx_min_loose = 1.5;

	const double pips_parameters[6][2] = {{5.83055147, 1.16171268},
                       {5.78469108, 0.80349681},
                       {5.844136,   0.53329847},
                       {5.81600646, 0.62782215},
                       {5.75400247, 0.88392704},
                       {5.4041441,  2.12325929}};

	const double pims_parameters[6][2] = {{ 6.72342364, 14.60912754},
                       { 6.85593523, 14.2618346 },
                       { 6.66919406, 14.64725313},
                       { 6.77248325, 14.33494692},
                       { 6.63926263, 14.4344934 },
                       { 6.85481318, 14.3687773 }};

	const double phi_theta_bowl_theta_min[6][4][2] = {{{8.02,15.42},  {7.33,13.70},  {7.14,12.08},  {7.04,10.39}},
                                {{7.80, 15.39}, {7.24, 13.69}, {7.09,  12.10},{6.87,10.28}},
                                {{7.66,15.49},  {7.14,13.73},  {6.93,12.08},  {6.96, 10.37}},
                                {{7.58,15.45},  {7.09,13.65},  {6.89,12.05},  {6.86,10.23}},
                                {{7.65,15.30},  {7.15,13.59},  {6.94,12.00},  {6.84,10.22}},
                                {{7.77,15.38},  {7.04,13.68},  {6.90,12.05},  {6.88,10.34}}};

	const double phi_theta_bowl_phi0[6][4][2] = {{{-18.02,18.11},   {-13.71,12.08},   {-10.36,10.19},   {-6.81,7.35}},
                                {{41.9, 78.25},    {46.0, 73.75},    {49.56, 70.3},    {53.2, 67.2}},
                                {{101.62,138.25},  {105.81,133.38},  {109.44,130.00},  {113.19,126.94}},
                                {{161.88,197.62},  {165.88,193.12},  {169.50,189.62},  {173.38,186.62}},
                                {{-137.75,-101.88},{-133.75,-106.44},{-130.25,-110.25},{-126.44,-112.81}},
                                {{-77.31,-41.91},  {-72.88,-46.41},  {-69.31,-49.88},  {-65.31,-52.38}}};

	const double theta_bowl_width = 0.9*550.;

	const int bins_p = 4;
	const double p_bin_edges[5] = { 1.25, 2.25, 2.50, 3.5, 5.00};
	const double p_bin_edges_3d[5] = { 1.25, 2.00, 2.50, 3.5, 5.00};
	
}

#endif

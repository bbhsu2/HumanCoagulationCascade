/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C++ Implementation by Bernard Hsu and anonymous
	3/27/2014
*/

#include <iostream>
#include <iomanip>
#include <cstring>
#include "Coagulation.h"
#include "RungeKutta.h"

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

// module level constants
const double Coagulation::v1 = 50000.0, Coagulation::k1 = 0.000001;
const double Coagulation::v2 = 50.0, Coagulation::k2 = 1.0;
const double Coagulation::v3 = 7.0, Coagulation::k3 = 10.0;
const double Coagulation::v4 = 7.0, Coagulation::k4 = 1.0;
const double Coagulation::v5 = 10.0, Coagulation::k5 = 10.0;
const double Coagulation::v6 = 0.1, Coagulation::k6 = 10.0;
const double Coagulation::v7 = 0.02, Coagulation::k7 = 10.0;
const double Coagulation::v8 = 2.0, Coagulation::k8 = 0.10;
const double Coagulation::v9 = 1.0e-9, Coagulation::k9 = 10.0;
const double Coagulation::v10 = 5.0e4, Coagulation::k10 = 10.0;
const double Coagulation::v11 = 50.0, Coagulation::k11 = 1.0;
const double Coagulation::v12 = 100.0, Coagulation::k12 = 10.0;
const double Coagulation::v13 = 9.0, Coagulation::k13 = 500.0;
const double Coagulation::v14 = 20000.0, Coagulation::k14 = 0.5;
const double Coagulation::v15 = 500.0, Coagulation::k15 = 500.0;
const double Coagulation::v16 = 7.0, Coagulation::k16 = 10.0;
const double Coagulation::v17 = 7.0, Coagulation::k17 = 10.0;
const double Coagulation::v18 = 7.0, Coagulation::k18 = 100.0;
const double Coagulation::v19 = 1.0, Coagulation::k19 = 1.0;
const double Coagulation::v20 = 7.0, Coagulation::k20 = 1.0;
const double Coagulation::v21 = 7.0, Coagulation::k21 = 5000.0;
const double Coagulation::v22 = 5.0, Coagulation::k22 = 10000.0;
const double Coagulation::v23 = 2.0, Coagulation::k23 = 1.0;
const double Coagulation::v24 = 7.0, Coagulation::k24 = 1.0;
const double Coagulation::v25 = 2.0, Coagulation::k25 = 1.0;
const double Coagulation::c26 = 0.01;
const double Coagulation::c27 = 0.5;
const double Coagulation::c28 = 0.5;
const double Coagulation::c29 = 0.5;
const double Coagulation::c30 = 0.1;
const double Coagulation::c31 = 0.5;
const double Coagulation::c32 = 0.5;
const double Coagulation::v33 = 70.0, Coagulation::k33 = 1.0;
const double Coagulation::v34 = 900.0, Coagulation::k34 = 200.0;
const double Coagulation::v35 = 70.0, Coagulation::k35 = 1.0;
const double Coagulation::v36 = 1000.0, Coagulation::k36 = 1.0;
const double Coagulation::c37 = 0.5;
const double Coagulation::v38 = 1.0, Coagulation::k38 = 10.0;
const double Coagulation::v39 = 1.0, Coagulation::k39 = 10.0;
const double Coagulation::v40 = 0.2, Coagulation::k40 = 10.0;
const double Coagulation::v41 = 7.0, Coagulation::k41 = 1.0;
const double Coagulation::v42 = 70.0, Coagulation::k42 = 1.0;
const double Coagulation::v43 = 7.0, Coagulation::k43 = 1.0;
const double Coagulation::c44 = 0.85 * 1.0 / 7.1; //this is only for Unfractionated Heparin//
const double Coagulation::c45 = 0.85;
const double Coagulation::c46 = 0.85 * 1.0; ////Need correction factor for unfractionated heparin
const double Coagulation::Imax = 1.0, Coagulation::IC50 = 0.34;
const double Coagulation::p_XII = 375.0, Coagulation::d_XII = 0.012;
const double Coagulation::d_XIIa = 20.0;
const double Coagulation::p_VIII = 0.70, Coagulation::d_VIII = 0.058;
const double Coagulation::d_VIIIa = 20.0;
const double Coagulation::p_IX = 89.60, Coagulation::d_IX = 0.029;
const double Coagulation::p_IXa = 0.0, Coagulation::d_IXa = 20.0;
const double Coagulation::p_XI = 30.6, Coagulation::d_XI = 0.1;
const double Coagulation::d_XIa = 20.0;
const double Coagulation::p_VII = 10.0, Coagulation::d_VII = 0.12;
const double Coagulation::d_VIIa = 20.0;
const double Coagulation::p_X = 174.3, Coagulation::d_X = 0.018;
const double Coagulation::d_Xa = 20.0;
const double Coagulation::p_V = 26.7, Coagulation::d_V = 0.043;
const double Coagulation::d_Va = 20.0;
const double Coagulation::d_Xa_Va = 20.0;
const double Coagulation::p_II = 1394.4, Coagulation::d_II = 0.010;
const double Coagulation::d_IIa = 67.4;
const double Coagulation::d_TAT = 0.2;
const double Coagulation::p_Fg = 8945.5, Coagulation::d_Fg = 0.032;
const double Coagulation::d_F = 0.050;
const double Coagulation::d_XF = 0.050;
const double Coagulation::d_FDP = 3.5;
const double Coagulation::d_D = 0.1;
const double Coagulation::p_XIII = 70.3, Coagulation::d_XIII = 0.0036;
const double Coagulation::d_XIIIa = 0.69;
const double Coagulation::p_Pg = 2154.30, Coagulation::d_Pg = 0.05;
const double Coagulation::d_P = 20.0;
const double Coagulation::p_PC = 60.0, Coagulation::d_PC = 0.050;
const double Coagulation::d_APC = 20.4;
const double Coagulation::p_Tmod = 50.0, Coagulation::d_Tmod = 0.050;
const double Coagulation::d_IIa_Tmod = 20.0;
const double Coagulation::d_IXa_VIIIa = 20.0;
const double Coagulation::d_TF = 0.05;
const double Coagulation::d_VII_TF = 0.7;
const double Coagulation::d_VIIa_TF = 20.0;
const double Coagulation::p_TFPI = 2.5, Coagulation::d_TFPI = 20.0;
const double Coagulation::d_Xa_TFPI = 20.0;
const double Coagulation::d_VIIa_TF_Xa_TFPI = 20.0;
const double Coagulation::p_PS = 300.0, Coagulation::d_PS = 0.0165;
const double Coagulation::d_APC_PS = 20.0;
const double Coagulation::p_Pk = 450.0, Coagulation::d_Pk = 0.05;
const double Coagulation::d_K = 20.0;
const double Coagulation::ka = 1.0;
const double Coagulation::Vd = 10.0;
const double Coagulation::ke = 0.693; //1 hour halflife assumed
//=========================Vit K Dependent factor initialization====================================
const double Coagulation::init_II = 1394.0;
const double Coagulation::init_VKH2 = 0.1;
const double Coagulation::init_VII = 10.0;
const double Coagulation::init_IX = 89.6;
const double Coagulation::init_X = 174.3;
const double Coagulation::init_PC = 60.0;
const double Coagulation::init_PS = 300.0;
const double Coagulation::d_VK = 0.2052;
const double Coagulation::d_VKH2 = 0.228;
const double Coagulation::VitK_k12 = 0.0587;
const double Coagulation::VitK_k21 = 0.0122;
const double Coagulation::Vc = 24.0;
const double Coagulation::init_VK = 1.0;
const double Coagulation::init_VKO = 0.1;
const double Coagulation::d_VKO = 0.228; 
double Coagulation::compartment = 1.0;

void Coagulation::Initialize(double (&vars)[num_factors])
{
	vars[Fg] = 8945.5; //Fibrinogen, I
	vars[F] = 0.0; //Fibrin = Ia
	vars[XF] = 0.0; //Crosslinked Fibrin
	vars[II] = 1394.4; //Prothrombin
	vars[IIa] = 0.0; //Thrombin
	vars[V] = 26.7;
	vars[Va] = 0.0;
	vars[VII] = 10.0; //Coag cascade entry point
	vars[VIIa] = 0.0;
	vars[VIII] = 0.7;
	vars[VIIIa] = 0.0;
	vars[IX] = 89.6;
	vars[IXa] = 0.0;
	vars[X] = 174.3;
	vars[Xa] = 0.0;
	vars[XI] = 30.6;
	vars[XIa] = 0.0;
	vars[XII] = 375.0;
	vars[XIIa] = 0.0;
	vars[XIII] = 70.3;
	vars[XIIIa] = 0.0;
	vars[Pg] = 2154.3;
	vars[P] = 0.0;
	vars[DegProd] = 0.0; //written as DP in the diagram
	vars[FDP] = 0.0;

	//==================Activators.  Set one of these to 0==================
	vars[CA] = 0.0; // Contact Activator for XII
	vars[TF] = 300.0;  //Tissue Factor

	//=========================Proteins=====================================
	vars[PC] = 60.0; //Protein C
	vars[APC] = 0.0; //Activted Protein C
	vars[PS] = 300.0; //Protein S
	vars[Pk] = 450.0; //Prekallikrien
	vars[K] = 0.0;  //Kallikrien
	vars[Tmod] = 50.0; //Thrombomodulin Activates II
	vars[TFPI] = 2.5; //Tissue Factor Pathway Inhibitor
	vars[D] = 0.0; //D-dimer

	//=========================VitK=========================================
	vars[VKH2] = 0.1;
	vars[VK] = 1.0;
	vars[VKO] = 0.1;  //exists in multiple isoforms
	vars[VK_p] = 0.0;

	//=========================Complexes====================================
	vars[Xa_TFPI] = 0.0;
	vars[VIIa_TF_Xa_TFPI] = 0.0;
	vars[VII_TF] = 0.0;
	vars[VIIa_TF] = 0.0;
	vars[Xa_Va] = 0.0;
	vars[IXa_VIIIa] = 0.0;
	vars[IIa_Tmod] = 0.0;
	vars[APC_PS] = 0.0;
	vars[TAT] = 0.0; //Thrombin-Antithrombin complex
	vars[ATIII_Heparin] = 0.0;
	vars[Xa_ATIII_Heparin] = 0.0;
	vars[IXa_ATIII_Heparin] = 0.0;
	vars[IIa_ATIII_Heparin] = 0.0;
	vars[Cwarf] = 0.0; //10.0
	vars[Awarf] = 0.0; //10.0 amount of warfarin
}

double Coagulation::Hyperbolic_rate_law(double v, double substrate, double enzyme, double k)
{
	return compartment * (v * enzyme * substrate) / (k + enzyme);
}

double Coagulation::VKH2mediated_factor_production(double vars[], double d_factor, double factor_initial)
{
	return compartment * d_factor * factor_initial * vars[VKH2] * (1.0 / init_VKH2);
}

double Coagulation::Irreversible_association(double s1, double s2, double cp)
{
	return compartment * s1 * s2 / cp;
}

double Coagulation::Factor_production(double initial, double degradation)
{
	return compartment * (initial * degradation);
}

double Coagulation::Factor_degradation(double degradation, double factor)
{
	return compartment * (degradation * factor);
}

double Coagulation::Warfarin_inhibited_first_order_kinetics(double Imax, double Cwarf, double IC50, double substrate, double degradation)
{
	return compartment * (degradation * substrate * (1.0 - Imax * Cwarf / (IC50 + Cwarf)));
}

double Coagulation::VK_Transport(double VK, double VK_p)
{
	return VitK_k12 * VK - VitK_k21 * (1.0 / Vc) * VK_p;
}

void Coagulation::coag_ode(double (&finish)[num_factors], double (&rate_basis)[num_factors])
{
	double r[49] = {};

	r[1] = Hyperbolic_rate_law(v1, rate_basis[VIII], rate_basis[IIa], k1);
	r[2] = Hyperbolic_rate_law(v2, rate_basis[VIIIa], rate_basis[APC_PS], k2);
	r[3] = Hyperbolic_rate_law(v3, rate_basis[IX], rate_basis[XIa], k3);
	r[4] = Hyperbolic_rate_law(v4, rate_basis[XI], rate_basis[XIIa], k4);
	r[5] = Hyperbolic_rate_law(v5, rate_basis[XI], rate_basis[IIa], k5);
	r[6] = Hyperbolic_rate_law(v6, rate_basis[VII], rate_basis[IIa], k6);
	r[7] = Hyperbolic_rate_law(v7, rate_basis[X], rate_basis[IXa], k7);
	r[8] = Hyperbolic_rate_law(v8, rate_basis[X], rate_basis[IXa_VIIIa], k8);
	r[9] = Hyperbolic_rate_law(v9, rate_basis[X], rate_basis[VIIa], k9);
	r[10] = Hyperbolic_rate_law(v10, rate_basis[V], rate_basis[IIa], k10);
	r[11] = Hyperbolic_rate_law(v11, rate_basis[Va], rate_basis[APC_PS], k11);
	r[12] = Hyperbolic_rate_law(v12, rate_basis[II], rate_basis[Xa_Va], k12);
	r[13] = Hyperbolic_rate_law(v13, rate_basis[II], rate_basis[Xa], k13);
	r[14] = Hyperbolic_rate_law(v14, rate_basis[Fg], rate_basis[IIa], k14);
	r[15] = Hyperbolic_rate_law(v15, rate_basis[Fg], rate_basis[P], k15);
	r[16] = Hyperbolic_rate_law(v16, rate_basis[F], rate_basis[XIIIa], k16);
	r[17] = Hyperbolic_rate_law(v17, rate_basis[F], rate_basis[P], k17);
	r[18] = Hyperbolic_rate_law(v18, rate_basis[XF], rate_basis[P], k18);
	r[19] = Hyperbolic_rate_law(v19, rate_basis[XF], rate_basis[APC_PS], k19);
	r[20] = Hyperbolic_rate_law(v20, rate_basis[XIII], rate_basis[IIa], k20);
	r[21] = Hyperbolic_rate_law(v21, rate_basis[Pg], rate_basis[IIa], k21);
	r[22] = Hyperbolic_rate_law(v22, rate_basis[Pg], rate_basis[F], k22);
	r[23] = Hyperbolic_rate_law(v23, rate_basis[Pg], rate_basis[APC_PS], k23);
	r[24] = Hyperbolic_rate_law(v24, rate_basis[PC], rate_basis[IIa_Tmod], k24);
	r[25] = Hyperbolic_rate_law(v25, rate_basis[Xa_Va], rate_basis[APC_PS], k25);
	r[26] = Irreversible_association(rate_basis[IXa], rate_basis[VIIIa], c26);
	r[27] = Irreversible_association(rate_basis[Xa], rate_basis[Va], c27);
	r[28] = Irreversible_association(rate_basis[IIa], rate_basis[Tmod], c28);
	r[29] = Irreversible_association(rate_basis[TF], rate_basis[VIIa], c29);
	r[30] = Irreversible_association(rate_basis[TF], rate_basis[VII], c30);
	r[31] = Irreversible_association(rate_basis[VIIa_TF], rate_basis[Xa_TFPI], c31);
	r[32] = Irreversible_association(rate_basis[Xa], rate_basis[TFPI], c32);
	r[33] = Hyperbolic_rate_law(v33, rate_basis[VII_TF], rate_basis[Xa], k33);
	r[34] = Hyperbolic_rate_law(v34, rate_basis[X], rate_basis[VIIa_TF], k34);
	r[35] = Hyperbolic_rate_law(v35, rate_basis[IX], rate_basis[VIIa_TF], k35);
	r[36] = Hyperbolic_rate_law(v36, rate_basis[VII_TF], rate_basis[TF], k36);
	r[37] = Irreversible_association(rate_basis[APC], rate_basis[PS], c37);
	r[38] = Hyperbolic_rate_law(v38, rate_basis[VII], rate_basis[Xa], k38);
	r[39] = Hyperbolic_rate_law(v39, rate_basis[VII], rate_basis[VIIa_TF], k39);
	r[40] = Hyperbolic_rate_law(v40, rate_basis[VII], rate_basis[IXa], k40);
	r[41] = Hyperbolic_rate_law(v41, rate_basis[XII], rate_basis[CA], k41);
	r[42] = Hyperbolic_rate_law(v42, rate_basis[XII], rate_basis[K], k42);
	r[43] = Hyperbolic_rate_law(v43, rate_basis[Pk], rate_basis[XIIa], k43);
	r[44] = Irreversible_association(rate_basis[IIa], rate_basis[ATIII_Heparin], c44);
	r[45] = Irreversible_association(rate_basis[Xa], rate_basis[ATIII_Heparin], c45);
	r[46] = Irreversible_association(rate_basis[IXa], rate_basis[ATIII_Heparin], c46);
	r[47] = Warfarin_inhibited_first_order_kinetics(Imax, rate_basis[Cwarf], IC50, rate_basis[VK], d_VK);
	r[48] = Warfarin_inhibited_first_order_kinetics(Imax, rate_basis[Cwarf], IC50, rate_basis[VKO], d_VKO);

	finish[XII] = Factor_production(p_XII, d_XII) - r[41] - r[42] - Factor_degradation(d_XII, rate_basis[XII]);
	finish[XIIa] = r[41] + r[42] - Factor_degradation(d_XIIa, rate_basis[XIIa]);
	finish[VIII] = Factor_production(p_VIII, d_VIII) - r[1] - Factor_degradation(d_VIII, rate_basis[VIII]);
	finish[VIIIa] = r[1] - r[2] - r[26] - Factor_degradation(d_VIIIa, rate_basis[VIIIa]);
	finish[IX] = VKH2mediated_factor_production(rate_basis, d_IX, init_IX) - r[35] - r[3] - Factor_degradation(d_IX, rate_basis[IX]);
	finish[IXa] = r[35] + r[3] - r[26] - r[46] - Factor_degradation(d_IXa, rate_basis[IXa]);
	finish[XI] = Factor_production(p_XI, d_XI) - r[5] - r[4] - Factor_degradation(d_XI, rate_basis[XI]);
	finish[XIa] = r[5] + r[4] - Factor_degradation(d_XIa, rate_basis[XIa]);
	finish[VII] = VKH2mediated_factor_production(rate_basis, d_VII, init_VII) - r[30] - r[6] - r[38] - r[39] - r[40] - Factor_degradation(d_VII, rate_basis[VII]);
	finish[VIIa] = -r[29] + r[6] + r[38] + r[39] + r[40] - Factor_degradation(d_VIIa, rate_basis[VIIa]);
	finish[X] = VKH2mediated_factor_production(rate_basis, d_X, init_X) - r[9] - r[34] - r[7] - r[8] - d_X * rate_basis[X];
	finish[Xa] = r[9] + r[34] + r[7] + r[8] - r[27] - r[32] - r[45] - Factor_degradation(d_Xa, rate_basis[Xa]);
	finish[V] = Factor_production(p_V, d_V) - r[10] - Factor_degradation(d_V, rate_basis[V]);
	finish[Va] = r[10] - r[11] - r[27] - Factor_degradation(d_Va, rate_basis[Va]);
	finish[Xa_Va] = r[27] - r[25] - Factor_degradation(d_Xa_Va, rate_basis[Xa_Va]);
	finish[II] = VKH2mediated_factor_production(rate_basis, d_II, init_II) - r[12] - r[13] - Factor_degradation(d_II, rate_basis[II]);
	finish[IIa] = r[12] + r[13] - r[28] - r[44] - Factor_degradation(d_IIa, rate_basis[IIa]);
	finish[TAT] = Factor_production(d_IIa, rate_basis[IIa]) - Factor_degradation(d_TAT, rate_basis[TAT]);
	finish[Fg] = Factor_production(p_Fg, d_Fg) - r[14] - r[15] - Factor_degradation(d_Fg, rate_basis[Fg]);
	finish[F] = r[14] - r[16] - r[17] - Factor_degradation(d_F, rate_basis[F]);
	finish[XF] = r[16] - r[18] - r[19] - Factor_degradation(d_XF, rate_basis[XF]);
	finish[FDP] = r[15] + r[17] + Factor_degradation(d_Fg, rate_basis[Fg]) + Factor_degradation(d_F, rate_basis[F]) - Factor_degradation(d_FDP, rate_basis[FDP]);
	finish[D] = r[18] + r[19] + Factor_degradation(d_XF, rate_basis[XF]) - Factor_degradation(d_D, rate_basis[D]);
	finish[XIII] = Factor_production(p_XIII, d_XIII) - r[20] - Factor_degradation(d_XIII, rate_basis[XIII]);
	finish[XIIIa] = r[20] - Factor_degradation(d_XIIIa, rate_basis[XIIIa]);
	finish[Pg] = Factor_production(p_Pg, d_Pg) - r[21] - r[22] - r[23] - Factor_degradation(d_Pg, rate_basis[Pg]);
	finish[P] = r[21] + r[23] + r[22] - Factor_degradation(d_P, rate_basis[P]);
	finish[PC] = VKH2mediated_factor_production(rate_basis, d_PC, init_PC) - r[24] - Factor_degradation(d_PC, rate_basis[PC]);
	finish[APC] = r[24] - r[37] - Factor_degradation(d_APC, rate_basis[APC]);
	finish[Tmod] = Factor_production(p_Tmod, d_Tmod) - r[28] - Factor_degradation(d_Tmod, rate_basis[Tmod]);
	finish[IIa_Tmod] = r[28] - Factor_degradation(d_IIa_Tmod, rate_basis[IIa_Tmod]);
	finish[IXa_VIIIa] = r[26] - Factor_degradation(d_IXa_VIIIa, rate_basis[IXa_VIIIa]);
	finish[TF] = -r[29] - r[30] - Factor_degradation(d_TF, rate_basis[TF]);
	finish[VII_TF] = r[30] - r[36] - r[33] - Factor_degradation(d_VII_TF, rate_basis[VII_TF]);
	finish[VIIa_TF] = r[29] + r[36] + r[33] - r[31] - Factor_degradation(d_VIIa_TF, rate_basis[VIIa_TF]);
	finish[TFPI] = Factor_production(p_TFPI, d_TFPI) - r[32] - Factor_degradation(d_TFPI, rate_basis[TFPI]);
	finish[Xa_TFPI] = r[32] - r[31] - Factor_degradation(d_Xa_TFPI, rate_basis[Xa_TFPI]);
	finish[VIIa_TF_Xa_TFPI] = r[31] - Factor_degradation(d_VIIa_TF_Xa_TFPI, rate_basis[VIIa_TF_Xa_TFPI]);
	finish[PS] = VKH2mediated_factor_production(rate_basis, d_PS, init_PS) - r[37] - Factor_degradation(d_PS, rate_basis[PS]);
	finish[APC_PS] = r[37] - Factor_degradation(d_APC_PS, rate_basis[APC_PS]);
	finish[Pk] = Factor_production(p_Pk, d_Pk) - r[43] - Factor_degradation(d_Pk, rate_basis[Pk]);
	finish[K] = r[43] - Factor_degradation(d_K, rate_basis[K]);
	finish[VK] = -r[47] + r[48] - VK_Transport(rate_basis[VK], rate_basis[VK_p]) + Factor_production(init_VK, d_VK) - Factor_degradation(d_VK, rate_basis[VK]);
	finish[VKH2] = r[47] - Factor_degradation(d_VKH2, rate_basis[VKH2]);
	finish[VKO] = -r[48] + Factor_degradation(d_VKH2, rate_basis[VKH2]);
	finish[VK_p] = VK_Transport(rate_basis[VK], rate_basis[VK_p]);
	finish[Awarf] = -ka * rate_basis[Awarf];
	finish[Cwarf] = ((ka * rate_basis[Awarf]) / Vd) - (ke * rate_basis[Cwarf]);
	finish[ATIII_Heparin] = -Irreversible_association(rate_basis[IIa], rate_basis[ATIII_Heparin], c44) - Irreversible_association(rate_basis[Xa], rate_basis[ATIII_Heparin], c45) - Irreversible_association(rate_basis[IXa], rate_basis[ATIII_Heparin], c46) - ke * rate_basis[ATIII_Heparin];
	finish[CA] = 0.0; //finish[CA] = dCA_dt(rate_basis)
	finish[DegProd] = finish[FDP] + finish[D];
}

void Coagulation::create_output(double (&vars)[num_factors], double t)
{

}

void Coagulation::coag_simulate()
{
	double coag_cur_state[num_factors] = {};
	double coag_next_state[num_factors] = {};
	double dt[] = { 1.0e-10, 1.0e-8, 1.0e-5 };
	double t_max[] = { 1.0e-8, 1.0e-6, 1.0e0 };
	int incr[] = { 10, 10, 100 };
	double t = 0.0;
	int c = 0;

	Initialize(coag_cur_state);

	//sw.WriteLine("time" + ',' + "XII" + ',' + " XIIa" + ',' + " VIII" + ',' + " VIIIa" + ',' + " IX" + ',' + " IXa" + ',' + " XI" + ',' + " XIa" + ',' + " VII" + ',' + " VIIa" + ',' + " X" + ',' + " Xa" + ',' + "V" + ',' + " Va" + ',' + " Xa_Va" + ',' + " II" + ',' + " IIa" + ',' + " TAT" + ',' + " Fg" + ',' + " F" + ',' + " XF" + ',' + " FDP" + ',' + " D" + ',' + " XIII" + ',' + "XIIIa" + ',' + "Pg" + ',' + "P" + ',' + "PC" + ',' + " APC" + ',' + " Tmod" + ',' + " IIa_Tmod" + ',' + " IXa_VIIIa" + ',' + " TF" + ',' + " VII_TF" + ',' + " VIIa_TF" + ',' + "TFPI" + ',' + "Xa_TFPI" + ',' + " VIIa_TF_Xa_TFPI" + ',' + " PS" + ',' + " APC_PS" + ',' + "Pk" + ',' + "K" + ',' + "VK" + ',' + "VKH2" + ',' + "VKO" + ',' + "VK_p" + ',' + "Awarf" + ',' + "Cwarf" + ',' + "ATIII_Heparin" + ',' + "DP");
	for (size_t i = 0; i < ARRAY_SIZE(dt); i++)
	{
		while (t < t_max[i])
		{
			rk2(coag_next_state, coag_cur_state, dt[i], coag_ode);
			memcpy(coag_cur_state, coag_next_state, sizeof(coag_cur_state));
			t += dt[i];
			c++;

			if (c == incr[i])
			{
				create_output(coag_cur_state, t);
				c = 0;
			}
		}
	}

	std::cout << std::setprecision(15);
	std::cout <<"Fg: " << coag_cur_state[Fg]<< std::endl;
	std::cout <<"F: " << coag_cur_state[F]<< std::endl;
	std::cout <<"XF: " << coag_cur_state[XF]<< std::endl;
	std::cout <<"DegProd: " << coag_cur_state[DegProd]<< std::endl;
	std::cout <<"II: " << coag_cur_state[II]<< std::endl;
	std::cout <<"IIa: " << coag_cur_state[IIa]<< std::endl;
	std::cout <<"IIa_Tmod: " << coag_cur_state[IIa_Tmod] << std::endl;
	std::cout <<"PC: " << coag_cur_state[PC]<< std::endl;
	std::cout <<"APC: " << coag_cur_state[APC]<< std::endl;
	std::cout <<"PS: " << coag_cur_state[PS]<< std::endl;
	std::cout <<"APC_PS: " << coag_cur_state[APC_PS]<< std::endl;
	std::cout <<"TAT: " << coag_cur_state[TAT]<< std::endl;
	std::cout <<"Pk: " << coag_cur_state[Pk]<< std::endl;
	std::cout <<"K: " << coag_cur_state[K]<< std::endl;
	std::cout <<"VII: " << coag_cur_state[VII]<< std::endl;
	std::cout <<"VII_TF: " << coag_cur_state[VII_TF]<< std::endl;
	std::cout <<"VIIa: " << coag_cur_state[VIIa]<< std::endl;
	std::cout <<"VIIa_TF: " << coag_cur_state[VIIa_TF]<< std::endl;
	std::cout <<"VIII: " << coag_cur_state[VIII]<< std::endl;
	std::cout <<"VIIIa: " << coag_cur_state[VIIIa]<< std::endl;
	std::cout <<"IX: " << coag_cur_state[IX]<< std::endl;
	std::cout <<"IXa: " << coag_cur_state[IXa]<< std::endl;
	std::cout <<"IXa_VIIIa: " << coag_cur_state[IXa_VIIIa]<< std::endl;
	std::cout <<"X: " << coag_cur_state[X]<< std::endl;
	std::cout <<"Xa: " << coag_cur_state[Xa]<< std::endl;
	std::cout <<"Xa_Va: " << coag_cur_state[Xa_Va]<< std::endl;
	std::cout <<"Xa_TFPI: " << coag_cur_state[Xa_TFPI]<< std::endl;
	std::cout <<"VIIa_TF_Xa_TFPI: " << coag_cur_state[VIIa_TF_Xa_TFPI]<< std::endl;
	std::cout <<"XI: " << coag_cur_state[XI]<< std::endl;
	std::cout <<"XIa: " << coag_cur_state[XIa]<< std::endl;
	std::cout <<"XII: " << coag_cur_state[XII]<< std::endl;
	std::cout <<"XIIa: " << coag_cur_state[XIIa]<< std::endl;
	std::cout << std::endl;
	std::cout <<"TF: " << coag_cur_state[TF]<< std::endl;
	std::cout <<"CA: " << coag_cur_state[CA]<< std::endl;
}

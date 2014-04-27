/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C++ Implementation by Bernard Hsu and anonymous
	3/27/2014
*/

#ifndef Coagulation_h
#define Coagulation_h

class Coagulation
{
	enum ECoagFactors
	{
		F,
		XF,
		II,
		IIa,
		V,
		Va,
		VII,
		VIIa,
		VIII,
		VIIIa,
		IX,
		IXa,
		X,
		Xa,
		XI,
		XIa,
		XII,
		XIIa,
		XIII,
		XIIIa,
		Pg,
		P,
		DegProd,
		FDP,
		//=========================Proteins========================================
		PC,
		APC,
		PS,
		Pk,
		K,
		CA,
		Tmod,
		TFPI,
		TF,
		D,
		//=========================VitK============================================
		VKH2,
		VK,
		VKO,
		VK_p,
		//=========================Complexes=======================================
		Xa_TFPI,
		VIIa_TF_Xa_TFPI,
		VII_TF,
		VIIa_TF,
		Xa_Va,
		IXa_VIIIa,
		IIa_Tmod,
		APC_PS,
		TAT,
		ATIII_Heparin,
		Xa_ATIII_Heparin,
		IXa_ATIII_Heparin,
		IIa_ATIII_Heparin,
		Cwarf,
		Awarf,
		Fg,
		num_factors
	};

public:
	// class level constants
	static const double v1, k1;
	static const double v2, k2;
	static const double v3, k3;
	static const double v4, k4;
	static const double v5, k5;
	static const double v6, k6;
	static const double v7, k7;
	static const double v8, k8;
	static const double v9, k9;
	static const double v10, k10;
	static const double v11, k11;
	static const double v12, k12;
	static const double v13, k13;
	static const double v14, k14;
	static const double v15, k15;
	static const double v16, k16;
	static const double v17, k17;
	static const double v18, k18;
	static const double v19, k19;
	static const double v20, k20;
	static const double v21, k21;
	static const double v22, k22;
	static const double v23, k23;
	static const double v24, k24;
	static const double v25, k25;
	static const double c26;
	static const double c27;
	static const double c28;
	static const double c29;
	static const double c30;
	static const double c31;
	static const double c32;
	static const double v33, k33;
	static const double v34, k34;
	static const double v35, k35;
	static const double v36, k36;
	static const double c37;
	static const double v38, k38;
	static const double v39, k39;
	static const double v40, k40;
	static const double v41, k41;
	static const double v42, k42;
	static const double v43, k43;
	static const double c44;
	static const double c45;
	static const double c46;
	static const double Imax, IC50;
	static const double p_XII, d_XII;
	static const double d_XIIa;
	static const double p_VIII, d_VIII;
	static const double d_VIIIa;
	static const double p_IX, d_IX;
	static const double p_IXa, d_IXa;
	static const double p_XI, d_XI;
	static const double d_XIa;
	static const double p_VII, d_VII;
	static const double d_VIIa;
	static const double p_X, d_X;
	static const double d_Xa;
	static const double p_V, d_V;
	static const double d_Va;
	static const double d_Xa_Va;
	static const double p_II, d_II;
	static const double d_IIa;
	static const double d_TAT;
	static const double p_Fg, d_Fg;
	static const double d_F;
	static const double d_XF;
	static const double d_FDP;
	static const double d_D;
	static const double p_XIII, d_XIII;
	static const double d_XIIIa;
	static const double p_Pg, d_Pg;
	static const double d_P;
	static const double p_PC, d_PC;
	static const double d_APC;
	static const double p_Tmod, d_Tmod;
	static const double d_IIa_Tmod;
	static const double d_IXa_VIIIa;
	static const double d_TF;
	static const double d_VII_TF;
	static const double d_VIIa_TF;
	static const double p_TFPI, d_TFPI;
	static const double d_Xa_TFPI;
	static const double d_VIIa_TF_Xa_TFPI;
	static const double p_PS, d_PS;
	static const double d_APC_PS;
	static const double p_Pk, d_Pk;
	static const double d_K;
	static const double ka;
	static const double Vd;
	static const double ke; //ke is 1 hour t1/2 assumed
	//=========================Vit K Dependent factor initialization====================================
	static const double init_II;
	static const double init_VKH2;
	static const double init_VII;
	static const double init_IX;
	static const double init_X;
	static const double init_PC;
	static const double init_PS;
	static const double d_VK;
	static const double d_VKH2;
	static const double VitK_k12;
	static const double VitK_k21;
	static const double Vc;
	static const double init_VK;
	static const double init_VKO;
	static const double d_VKO;
	static double compartment;

	static double Hyperbolic_rate_law(double v, double substrate, double enzyme, double k);
	static double VKH2mediated_factor_production(double vars[], double d_factor, double factor_initial);
	static double Irreversible_association(double s1, double s2, double cp);
	static double Factor_production(double initial, double degradation);
	static double Factor_degradation(double degradation, double factor);
	static double Warfarin_inhibited_first_order_kinetics(double Imax, double Cwarf, double IC50, double substrate, double degradation);
	static double VK_Transport(double VK, double VK_p);
	static void coag_ode(double (&out)[num_factors], double (&rate_basis)[num_factors]);
	static void create_output(double (&vars)[num_factors], double t);
	static void coag_simulate();

private:
	static void Initialize(double (&vars)[num_factors]);
};

#endif

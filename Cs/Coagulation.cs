/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C# Implementation by Bernard Hsu
	3/23/2014
*/

using System;
using System.IO;

namespace Coag
{
    public class Coagulation
    {
        static string filename = "coag.csv";
        public StreamWriter sw = new StreamWriter(filename);
        //indexers
        private int Fg = 54 - 1;
        private int F = 1 - 1;
        private int XF = 2 - 1;
        private int II = 3 - 1;
        private int IIa = 4 - 1;
        private int V = 5 - 1;
        private int Va = 6 - 1;
        private int VII = 7 - 1;
        private int VIIa = 8 - 1;
        private int VIII = 9 - 1;
        private int VIIIa = 10 - 1;
        private int IX = 11 - 1;
        private int IXa = 12 - 1;
        private int X = 13 - 1;
        private int Xa = 14 - 1;
        private int XI = 15 - 1;
        private int XIa = 16 - 1;
        private int XII = 17 - 1;
        private int XIIa = 18 - 1;
        private int XIII = 19 - 1;
        private int XIIIa = 20 - 1;
        private int Pg = 21 - 1;
        private int P = 22 - 1;
        private int DegProd = 23 - 1;
        private int FDP = 24 - 1;
        //=========================Proteins========================================
        private int PC = 25 - 1;
        private int APC = 26 - 1;
        private int PS = 27 - 1;
        private int Pk = 28 - 1;
        private int K = 29 - 1;
        private int CA = 30 - 1;
        private int Tmod = 31 - 1;
        private int TFPI = 32 - 1;
        private int TF = 33 - 1;
        private int D = 34 - 1;
        //=========================VitK============================================
        private int VKH2 = 35 - 1;
        private int VK = 36 - 1;
        private int VKO = 37 - 1;
        private int VK_p = 38 - 1;
        //=========================Complexes=======================================
        private int Xa_TFPI = 39 - 1;
        private int VIIa_TF_Xa_TFPI = 40 - 1;
        private int VII_TF = 41 - 1;
        private int VIIa_TF = 42 - 1;
        private int Xa_Va = 43 - 1;
        private int IXa_VIIIa = 44 - 1;
        private int IIa_Tmod = 45 - 1;
        private int APC_PS = 46 - 1;
        private int TAT = 47 - 1;
        private int ATIII_Heparin = 48 - 1;
        private int Xa_ATIII_Heparin = 49 - 1;
        private int IXa_ATIII_Heparin = 50 - 1;
        private int IIa_ATIII_Heparin = 51 - 1;
        private int Cwarf = 52 - 1;
        private int Awarf = 53 - 1;
        private int num_factors = 54;
        
        //class level complexes
        const double v1 = 50000.0, k1 = 0.000001;
        const double v2 = 50.0, k2 = 1.0;
        const double v3 = 7.0, k3 = 10.0;
        const double v4 = 7.0, k4 = 1.0;
        const double v5 = 10.0, k5 = 10.0;
        const double v6 = 0.1, k6 = 10.0;
        const double v7 = 0.02, k7 = 10.0;
        const double v8 = 2.0, k8 = 0.10;
        const double v9 = 1.0e-9, k9 = 10.0;
        const double v10 = 5.0e4, k10 = 10.0;
        const double v11 = 50.0, k11 = 1.0;
        const double v12 = 100.0, k12 = 10.0;
        const double v13 = 9.0, k13 = 500.0;
        const double v14 = 20000.0, k14 = 0.5;
        const double v15 = 500.0, k15 = 500.0;
        const double v16 = 7.0, k16 = 10.0;
        const double v17 = 7.0, k17 = 10.0;
        const double v18 = 7.0, k18 = 100.0;
        const double v19 = 1.0, k19 = 1.0;
        const double v20 = 7.0, k20 = 1.0;
        const double v21 = 7.0, k21 = 5000.0;
        const double v22 = 5.0, k22 = 10000.0;
        const double v23 = 2.0, k23 = 1.0;
        const double v24 = 7.0, k24 = 1.0;
        const double v25 = 2.0, k25 = 1.0;
        const double c26 = 0.01;
        const double c27 = 0.5;
        const double c28 = 0.5;
        const double c29 = 0.5;
        const double c30 = 0.1;
        const double c31 = 0.5;
        const double c32 = 0.5;
        const double v33 = 70.0, k33 = 1.0;
        const double v34 = 900.0, k34 = 200.0;
        const double v35 = 70.0, k35 = 1.0;
        const double v36 = 1000.0, k36 = 1.0;
        const double c37 = 0.5;
        const double v38 = 1.0, k38 = 10.0;
        const double v39 = 1.0, k39 = 10.0;
        const double v40 = 0.2, k40 = 10.0;
        const double v41 = 7.0, k41 = 1.0;
        const double v42 = 70.0, k42 = 1.0;
        const double v43 = 7.0, k43 = 1.0;
        const double c44 = 0.85 * 1.0 / 7.1; 
        const double c45 = 0.85;
        const double c46 = 0.85 * 1.0; 
        const double Imax = 1.0, IC50 = 0.34;
        const double p_XII = 375.0, d_XII = 0.012;
        const double d_XIIa = 20.0;
        const double p_VIII = 0.70, d_VIII = 0.058;
        const double d_VIIIa = 20.0;
        const double p_IX = 89.60, d_IX = 0.029;
        const double p_IXa = 0.0, d_IXa = 20.0;
        const double p_XI = 30.6, d_XI = 0.1;
        const double d_XIa = 20.0;
        const double p_VII = 10.0, d_VII = 0.12;
        const double d_VIIa = 20.0;
        const double p_X = 174.3, d_X = 0.018;
        const double d_Xa = 20.0;
        const double p_V = 26.7, d_V = 0.043;
        const double d_Va = 20.0;
        const double d_Xa_Va = 20.0;
        const double p_II = 1394.4, d_II = 0.010;
        const double d_IIa = 67.4;
        const double d_TAT = 0.2;
        const double p_Fg = 8945.5, d_Fg = 0.032;
        const double d_F = 0.050;
        const double d_XF = 0.050;
        const double d_FDP = 3.5;
        const double d_D = 0.1;
        const double p_XIII = 70.3, d_XIII = 0.0036;
        const double d_XIIIa = 0.69;
        const double p_Pg = 2154.30, d_Pg = 0.05;
        const double d_P = 20.0;
        const double p_PC = 60.0, d_PC = 0.050;
        const double d_APC = 20.4;
        const double p_Tmod = 50.0, d_Tmod = 0.050;
        const double d_IIa_Tmod = 20.0;
        const double d_IXa_VIIIa = 20.0;
        const double d_TF = 0.05;
        const double d_VII_TF = 0.7;
        const double d_VIIa_TF = 20.0;
        const double p_TFPI = 2.5, d_TFPI = 20.0;
        const double d_Xa_TFPI = 20.0;
        const double d_VIIa_TF_Xa_TFPI = 20.0;
        const double p_PS = 300.0, d_PS = 0.0165;
        const double d_APC_PS = 20.0;
        const double p_Pk = 450.0, d_Pk = 0.05;
        const double d_K = 20.0;
        const double ka = 1.0;
        const double Vd = 10.0;
        const double ke = 0.693; //1 hour half life

        //Vit K Dependent factor initialization
        const double init_II = 1394.0;
        const double init_VKH2 = 0.1;
        const double init_VII = 10.0;
        const double init_IX = 89.6;
        const double init_X = 174.3;
        const double init_PC = 60.0;
        const double init_PS = 300.0;
        const double d_VK = 0.2052;
        const double d_VKH2 = 0.228;
        const double VitK_k12 = 0.0587;
        const double VitK_k21 = 0.0122;
        const double Vc = 24.0;
        const double init_VK = 1.0;
        const double init_VKO = 0.1;
        const double d_VKO = 0.228;

        double compartment = 1.0; 

        private void Initialize(double[] vars)
        {
            vars[Fg] = 8945.5; 
            vars[F] = 0.0; 
            vars[XF] = 0.0;
            vars[II] = 1394.4; 
            vars[IIa] = 0.0; 
            vars[V] = 26.7;
            vars[Va] = 0.0;
            vars[VII] = 10.0; 
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
            vars[DegProd] = 0.0;
            vars[FDP] = 0.0;

            //Cascade entrypoint. Only set one of these to 0
            vars[CA] = 300.0; // Contact Activator for XII
            vars[TF] = 0.0;  //Tissue Factor

            //Proteins
            vars[PC] = 60.0; 
            vars[APC] = 0.0; 
            vars[PS] = 300.0; 
            vars[Pk] = 450.0; 
            vars[K] = 0.0;  
            vars[Tmod] = 50.0;
            vars[TFPI] = 2.5; 
            vars[D] = 0.0; 

            //VitK
            vars[VKH2] = 0.1;
            vars[VK] = 1.0;
            vars[VKO] = 0.1; 
            vars[VK_p] = 0.0;

            //Complexes
            vars[Xa_TFPI] = 0.0;
            vars[VIIa_TF_Xa_TFPI] = 0.0;
            vars[VII_TF] = 0.0;
            vars[VIIa_TF] = 0.0;
            vars[Xa_Va] = 0.0;
            vars[IXa_VIIIa] = 0.0;
            vars[IIa_Tmod] = 0.0;
            vars[APC_PS] = 0.0;
            vars[TAT] = 0.0; 
            vars[ATIII_Heparin] = 0.0;
            vars[Xa_ATIII_Heparin] = 0.0;
            vars[IXa_ATIII_Heparin] = 0.0;
            vars[IIa_ATIII_Heparin] = 0.0;
            vars[Cwarf] = 0.0; 
            vars[Awarf] = 0.0; 
        }

        public double Hyperbolic_rate_law(double v, double substrate, double enzyme, double k)
        {
            return compartment * (v * enzyme * substrate) / (k + enzyme);
        }
        
        public double VKH2mediated_factor_production(double[] vars, double d_factor, double factor_initial)
        {
            return compartment * d_factor * factor_initial * vars[VKH2] * (1.0 / init_VKH2);
        }
        
        public double Irreversible_association(double s1, double s2, double cp)
        {
            return compartment * s1 * s2 / cp;
        }
        
        public double Factor_production(double initial, double degradation)
        {
            return compartment * (initial * degradation);
        }
        
        public double Factor_degradation(double degradation, double factor)
        {
            return compartment * (degradation * factor);
        }
        
        public double Warfarin_inhibited_first_order_kinetics(double Imax, double Cwarf, double IC50, double substrate, double degradation)
        {
            return compartment * (degradation * substrate * (1 - Imax * Cwarf / (IC50 + Cwarf)));
        }
        
        public double VK_Transport(double VK, double VK_p)
        {
            return VitK_k12 * VK - VitK_k21 * (1.0 / Vc) * VK_p;
        }
        
        public void create_output(double[] vars, double t)
        {
            sw.WriteLine(t + "," + vars[XII] + "," + vars[XIIa] + "," + + vars[VIII] + ", " + vars[VIIIa] + ", " + vars[IX] + ", " + vars[IXa] + ", " + vars[XI] + ", " + vars[XIa] + "," + vars[VII] + "," + vars[VIIa] + "," + vars[X] + "," + vars[Xa] + "," + vars[V] + "," + vars[Va] + "," + vars[Xa_Va] + "," + vars[II] + "," + vars[IIa] + "," + vars[TAT] + "," + vars[Fg] + "," + vars[F] + "," + vars[XF] + "," + vars[FDP] + "," + vars[D] + "," + vars[XIII] + "," + vars[XIIIa] + "," + vars[Pg] + "," + vars[P] + "," + vars[PC] + "," + vars[APC] + "," + vars[Tmod] + "," + vars[IIa_Tmod] + "," + vars[IXa_VIIIa] + "," + vars[TF] + "," + vars[VII_TF] + "," + vars[VIIa_TF] + "," + vars[TFPI] + "," + vars[Xa_TFPI] + "," + vars[VIIa_TF_Xa_TFPI] + "," + vars[PS] + "," + vars[APC_PS] + "," + vars[Pk] + "," + vars[K] + "," + vars[VK] + "," + vars[VKH2] + "," + vars[VKO] + "," + vars[VK_p] + "," + vars[Awarf] + "," + vars[Cwarf] + "," + vars[ATIII_Heparin] + "," + vars[DegProd]);
        }


        public double[] coag_ode(double[] rate_basis)
        {
            double[] r = new double[49];
            double[] finish = new double[rate_basis.Length];

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
            //finish[CA] = dCA_dt(rate_basis)
            finish[DegProd] = finish[FDP] + finish[D];

            return finish;
        }

        public void coag_simulate()
        {
            double[] coag_cur_state = new double[num_factors];
            double[] dt = new double[] { 1.0e-10, 1.0e-8, 1.0e-5 };
            double[] t_max = new double[] { 1.0e-8, 1.0e-6, 1.0e0 };
            int[] incr = new int[] { 10, 10, 100 };
            double t = 0.0;
            int c = 0;

            this.Initialize(coag_cur_state);

            sw.WriteLine("time" + ',' + "XII" + ',' + " XIIa" + ',' + " VIII" + ',' + " VIIIa" + ',' + " IX" + ',' + " IXa" + ',' + " XI" + ',' + " XIa" + ',' + " VII" + ',' + " VIIa" + ',' + " X" + ',' + " Xa" + ',' + "V" + ',' + " Va" + ',' + " Xa_Va" + ',' + " II" + ',' + " IIa" + ',' + " TAT" + ',' + " Fg" + ',' + " F" + ',' + " XF" + ',' + " FDP" + ',' + " D" + ',' + " XIII" + ',' + "XIIIa" + ',' + "Pg" + ',' + "P" + ',' + "PC" + ',' + " APC" + ',' + " Tmod" + ',' + " IIa_Tmod" + ',' + " IXa_VIIIa" + ',' + " TF" + ',' + " VII_TF" + ',' + " VIIa_TF" + ',' + "TFPI" + ',' + "Xa_TFPI" + ',' + " VIIa_TF_Xa_TFPI" + ',' + " PS" + ',' + " APC_PS" + ',' + "Pk" + ',' + "K" + ',' + "VK" + ',' + "VKH2" + ',' + "VKO" + ',' + "VK_p" + ',' + "Awarf" + ',' + "Cwarf" + ',' + "ATIII_Heparin" + ',' + "DP");

            for (int i = 0; i < dt.Length; i++)
            {
                while (t < t_max[i])
                {
                    coag_cur_state = RungeKutta.rk2(coag_cur_state, dt[i], coag_ode);
                    t += dt[i];
                    c++;

                    if (c == incr[i])
                    {
                        create_output(coag_cur_state, t);
                        c = 0;
                    }
                }
            }

            Console.WriteLine("Fg:" + coag_cur_state[Fg]);
            Console.WriteLine("F:" + coag_cur_state[F]);
            Console.WriteLine("XF:" + coag_cur_state[XF]);
            Console.WriteLine("DegProd:" + coag_cur_state[DegProd]);
            Console.WriteLine("II:" + coag_cur_state[II]);
            Console.WriteLine("IIa:" + coag_cur_state[IIa]);
            Console.WriteLine("IIa_Tmod:" + coag_cur_state[IIa_Tmod]);
            Console.WriteLine("PC:" + coag_cur_state[PC]);
            Console.WriteLine("APC:" + coag_cur_state[APC]);
            Console.WriteLine("PS:" + coag_cur_state[PS]);
            Console.WriteLine("APC_PS:" + coag_cur_state[APC_PS]);
            Console.WriteLine("TAT:" + coag_cur_state[TAT]);
            Console.WriteLine("Pk:" + coag_cur_state[Pk]);
            Console.WriteLine("K:" + coag_cur_state[K]);
            Console.WriteLine("VII:" + coag_cur_state[VII]);
            Console.WriteLine("VII_TF:" + coag_cur_state[VII_TF]);
            Console.WriteLine("VIIa:" + coag_cur_state[VIIa]);
            Console.WriteLine("VIIa_TF:" + coag_cur_state[VIIa_TF]);
            Console.WriteLine("VIII:" + coag_cur_state[VIII]);
            Console.WriteLine("VIIIa:" + coag_cur_state[VIIIa]);
            Console.WriteLine("IX:" + coag_cur_state[IX]);
            Console.WriteLine("IXa:" + coag_cur_state[IXa]);
            Console.WriteLine("IXa_VIIIa:" + coag_cur_state[IXa_VIIIa]);
            Console.WriteLine("X:" + coag_cur_state[X]);
            Console.WriteLine("Xa:" + coag_cur_state[Xa]);
            Console.WriteLine("Xa_Va:" + coag_cur_state[Xa_Va]);
            Console.WriteLine("Xa_TFPI:" + coag_cur_state[Xa_TFPI]);
            Console.WriteLine("VIIa_TF_Xa_TFPI:" + coag_cur_state[VIIa_TF_Xa_TFPI]);
            Console.WriteLine("XI:" + coag_cur_state[XI]);
            Console.WriteLine("XIa:" + coag_cur_state[XIa]);
            Console.WriteLine("XII:" + coag_cur_state[XII]);
            Console.WriteLine("XIIa:" + coag_cur_state[XIIa]);
            Console.WriteLine();
            Console.WriteLine("TF:" + coag_cur_state[TF]);
            Console.WriteLine("CA:" + coag_cur_state[CA]);
            sw.Close();
        }
    }
}


/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C# Implementation by Bernard Hsu
	3/23/2014
*/

using System;

namespace Coag
{
    public static class RungeKutta
    {
        public delegate double[] RkDelegate(double[] x);

        static double sixth = 1.0 / 6.0;

        public static double[] rk2(double[] x, double dt, RkDelegate f)
        {
            double[] r = new double[x.Length];
            double[] k1 = new double[x.Length];
            double[] k2 = new double[x.Length];
            double[] fullstep = new double[x.Length];
            double[] halfdt = new double[x.Length];
            for (int i = 0; i < halfdt.Length; i++)
            {
                fullstep[i] = dt;
                halfdt[i] = 0.5 * dt;
            }

            k1 = f(x);
            k2 = f(ArrayCombiner(x, k1, halfdt));
            r = ArrayCombiner(x, k2, fullstep);
            return r;
        }

        private static double[] ArrayCombiner(double[] x, double[] k, double[] dt)
        {
            for (int i = 0; i < k.Length; i++)
            {
                k[i] *= dt[i];
                k[i] += x[i];
            }
            return k;
        }
    }
}


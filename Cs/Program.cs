/*
	Human Coagulation Cascade
	http://github.com/bbhsu2/HumanCoagulationCascade
	C# Implementation by Bernard Hsu
	3/23/2014
*/

using System;

namespace Coag
{
    class Program
    {
        static void Main(string[] args)
        {
            Coagulation coag = new Coagulation();
            coag.coag_simulate();
            Console.ReadLine();
        }
    }
}

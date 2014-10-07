using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab7
{
    class DifEquation
    {
        public delegate double func2D(double x, double y);
        public delegate double func3D(double x);

        private double x0    { get; set; }
        private double y0    { get; set; }
        private double h     { get; set; }
        private func3D deriv { get; set; }

        DifEquation(double x0, double y0, func3D deriv, double h)
        {
            this.x0 = x0;
            this.y0 = y0;
            this.h = h;
            this.deriv = deriv;
        }
    }
}
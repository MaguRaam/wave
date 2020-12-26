#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <blitz/array.h>

using namespace blitz;

std::string int_to_string(unsigned int value, const unsigned int digits)
{
    std::string lc_string = std::to_string(value);

    if (lc_string.size() < digits)
    {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                  ? 1
                                                  : 0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }

    return lc_string;
}

//write data tecplot:
void write(const Array<double,1>& u,const Array<double,1>& x, double t, int nt)
{
    std::ofstream tpl;
    const std::string filename = "../plot/plot_" + int_to_string(nt, 3) + ".dat";
    tpl.open(filename);
    tpl.flags(std::ios::dec | std::ios::scientific);  
    tpl.precision(6);
    tpl << "TITLE = \"Wave Equation 1D\" " << std::endl
    << "VARIABLES = \"x\", \"u\" " << std::endl;
    tpl << "Zone I = " << x.numElements() << std::endl;
    tpl << "SOLUTIONTIME = " << t << std::endl;

    for (int i = 0; i < x.numElements(); i++)
    tpl << x(i) << "\t" << u(i)<< std::endl;
}
 

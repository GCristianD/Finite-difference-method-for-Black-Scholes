//
//  Option.cpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#include "Option.hpp"
#include <cmath>

// The formulas for the payoff of a european put option.
// See https://en.wikipedia.org/wiki/Put_option
double PutOption::Pay(double z)
{
    if (z>K) {return 0.0;}
    else {return K-z;}
}

double PutOption::Upper(BlackScholes *PModel, double t)
{   return 0.0;
}

double PutOption::Lower(BlackScholes *PModel, double t)
{
    return exp(-(T-t)*PModel->r )*K;
}
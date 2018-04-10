//
//  ParabolicPDE.cpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#include "ParabolicPDE.hpp"
#include <iostream>
#include <cmath>

BlackScholesEq::BlackScholesEq(BlackScholes* PModel_, Option* POption_ )
{
    PModel=PModel_;
    POption=POption_;
    T= POption->T;
    xlow= POption->zlow;
    xup= POption->zup;
    
}

// Defines the formulas of the coefficients of the Black-Scholes PDE,
// see https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation

double BlackScholesEq::a1(double t, double x)
{    return -pow(x* PModel->sigma,2.0)/2; }

double BlackScholesEq::a2(double t, double x)
{   return -x* PModel->r; }

double BlackScholesEq::a3(double t, double x)
{   return PModel->r; }

double BlackScholesEq::a4(double t, double x)
{   return 0.0; }

double BlackScholesEq::f(double x)
{   return POption->Pay(x); }

double BlackScholesEq::fup(double t)
{   return POption->Upper(PModel, t); }

double BlackScholesEq::flow(double t)
{   return POption->Lower(PModel, t); }


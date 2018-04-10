//
//  ParabolicPDE.hpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#ifndef ParabolicPDE_hpp
#define ParabolicPDE_hpp

#include <stdio.h>
#include "Option.hpp"

// This class includes virtual functions for the coefficients and boundary conditions of a parabolic partial differential equation in variables (t,x). It also specifies the domain [0,T]*[xlow,xup] through the variables T,xlow,xup.
class ParabEq
{
public:
    double xlow, xup,T;
    
    virtual double a1(double t, double x)=0;
    virtual double a2(double t, double x)=0;
    virtual double a3(double t, double x)=0;
    virtual double a4(double t, double x)=0;
    
    virtual double f(double x)=0;
    virtual double fup(double t)=0;
    virtual double flow(double t)=0;
};

// The Black Scholes equation is a particular case of a parabolic PDE. PModel and POption are used to pass the parameters of the model. https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation
class BlackScholesEq: public ParabEq
{
public:
    BlackScholes* PModel;
    Option* POption;
    
    BlackScholesEq(BlackScholes* PModel_, Option* POption_ );
    
    double a1(double t, double x);
    double a2(double t, double x);
    double a3(double t, double x);
    double a4(double t, double x);
    
    double f(double x);
    double fup(double t);
    double flow(double t);    
};



#endif /* ParabolicPDE_hpp */

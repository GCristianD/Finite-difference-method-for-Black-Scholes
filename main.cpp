//
//  main.cpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//


#include <iostream>
#include "Option.hpp"
#include "ParabolicPDE.hpp"
#include "FDMethod.hpp"
using namespace std;

//  We implement the pricing of a European Put option (https://en.wikipedia.org/wiki/Put_option ) based on numerically solving the Black-Scholes Partial Differential Equation https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_equation by the explicit finite difference method.
//  T = expiry date,  K = strike price, S = price of the stock, r = interest rate, sigma=volatility

int main()
{
    
    double S=90.0, sigma=0.2, r=0.04;
    double T=0.08, K=90.0, zlow=0.0, zup=2.0*S;
    int Tpoints=3000, Xpoints=1000;           // Numbers of points in the discretizations
    
    BlackScholes Model(S ,sigma, r);
    PutOption Put(K,T,zlow,zup);
    BlackScholesEq BS(&Model,&Put);
    
    ExpMet Met(&BS,Tpoints,Xpoints);
    Met.Solve();
    
    cout << "Cost of option: " << Met.v(0.0,S) << endl;
    
    return 0;
}

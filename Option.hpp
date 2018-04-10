//
//  Option.hpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#ifndef Option_hpp
#define Option_hpp

#include <iostream>
#include <stdio.h>

// This class stores the parameters S0, r, sigma of the Black-Scholes model.
class BlackScholes
{
public:
    double S0, sigma, r;
    BlackScholes(double S_, double sigma_, double r_) {S0=S_; sigma=sigma_; r=r_;   }
};

// This class contains virtual functions for the payoff and boundary conditions (Upper and Lower)
class Option
{
public:
    
    double zlow,zup,T;
    
    virtual double Pay(double z)=0;
    virtual double Upper(BlackScholes *PModel, double t)=0;
    virtual double Lower(BlackScholes *PModel, double t)=0;
};


// This class implements the payoff and the boundary conditions of the Put Options. For the explanation of the formulas see https://en.wikipedia.org/wiki/Put_option 
class PutOption: public Option
{
public:
    double K;
    PutOption(double K_, double T_, double zlow_, double zup_)
        {K=K_; T=T_; zlow=zlow_; zup=zup_; }
    double Pay(double z);
    double Upper(BlackScholes* PModel, double t);
    double Lower(BlackScholes* PModel, double t);
};





#endif /* Option_hpp */

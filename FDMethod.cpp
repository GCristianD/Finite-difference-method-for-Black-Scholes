//
//  FDMethod.cpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#include "FDMethod.hpp"



FDMethod::FDMethod(ParabEq* PPDE_, int Tpoints, int Xpoints)
{
    PPDE=PPDE_;
    max_i=Tpoints;
    max_j=Xpoints;
    // Computes the step-size given the size of the domain and the number of steps
    dx=(PPDE->xup-PPDE->xlow)/max_j;
    dt=PPDE->T/max_i;
    W.resize(max_i+1);
    for (int i=0; i<=max_i; i++) W[i].resize(max_j+1);
}

// Approximates the solutions at points (t,x) which are not on the grid by taking a weighted average of the nearby 4 points on the grid. Assumes W is already defined (i.e. Solve was called).
double FDMethod::v(double t, double x)
{
    int i=(int)(t/dt);
    int j=(int)((x-PPDE->xlow)/dx);
    
    
    // Computes the weights q1,q0,p1,p0
    double q1=(t-FDMethod::t(i))/dt;
    double q0=1.0-q1;
    double p1=(x-FDMethod::x(j))/dx;
    double p0=1.0-p1;
    return q1*p1*W[i+1][j+1]+q1*p0*W[i+1][j]+q0*p1*W[i][j+1] + q0*p0*W[i][j]; // The weighted average
}


// Solves the Partial Differential Equation by using the explicit discretization method.
void ExpMet::Solve()
{
    for (int j=0; j<=max_j; j++) W[max_i][j]=f(j);  // This is computed from the given boundary condition
    for (int i=max_i; i>0; i--)
    {
        // Boundary conditions
        W[i-1][0]=flow(i-1);
        W[i-1][max_j]=fup(i-1);
        
        // We compute the interior values recursively
        for (int j=1; j<max_j; j++)
        {
            W[i-1][j]=C1(i,j)*W[i][j-1]+C2(i,j)*W[i][j]+C3(i,j)*W[i][j+1]+C4(i,j);
        }
    }
}
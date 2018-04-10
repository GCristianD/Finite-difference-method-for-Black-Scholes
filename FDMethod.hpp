//
//  FDMethod.hpp
//  FiniteDiff
//
//  Created by Cristian Gavrus on 1/17/18.
//  Copyright © 2018 Cristian Gavrus. All rights reserved.
//

#ifndef FDMethod_hpp
#define FDMethod_hpp

#include <stdio.h>
#include <vector>
#include "ParabolicPDE.hpp"
using namespace std;

// This class holds the functions that define the coefficients a1-a4 of the PDE and
// the boundary conditions f,fup,flow. 
class FDMethod
{
public:
    int max_i, max_j;          // Numbers of points on the t and x axes in our discretization
    double dx,dt;              // The step sizes
    vector<vector<double>> W;  // W is a matrix storing the values of the solution on the grid points.
    ParabEq* PPDE;             // PPDE is used for passing functions from the ParabPDE class
    
    FDMethod(ParabEq* PPDE_, int Tpoints, int Xpoints);
    
    // Computes the i-th (resp. j-th) point on that axis (on the grid)
    double t(double i) { return dt*i;}
    double x(int j) { return PPDE->xlow+dx*j;}
    
    // Returns the the PDE coefficients and boundary values at the points on the grid. 
    double a1(double i,int j) { return PPDE->a1(t(i),x(j));  }
    double a2(double i,int j) { return PPDE->a2(t(i),x(j));  }
    double a3(double i,int j) { return PPDE->a3(t(i),x(j));  }
    double a4(double i,int j) { return PPDE->a4(t(i),x(j));  }
    double f(int j) { return PPDE->f(x(j));  }
    double fup(int i) { return PPDE->fup(t(i));  }
    double flow(int i) { return PPDE->flow(t(i));  }
    
    // Approximates the solutions at points not on the grid by taking a weighted average of the
    // nearby 4 points on the grid
    double v(double t, double x);
    
};

// This subclass of FDMethod is used for implementing the explicit method to solve the PDE.
class ExpMet: public FDMethod
{
public:
    ExpMet(ParabEq* PPDE_, int Tpoints, int Xpoints): FDMethod(PPDE_, Tpoints, Xpoints) {}
    
    
    // C1-C4 are the coefficients of the difference equations used to recursively approximate the
    // values of the solutions on the grid.
    double C1(int i, int j) { return (a2(i,j)/2.0 - a1(i,j)/dx)*dt/dx;      }
    double C2(int i, int j) { return 1.0-a3(i,j)*dt+a1(i,j)/(dx*dx)*2.0*dt; }
    double C3(int i, int j) { return -(a2(i,j)/2.0 + a1(i,j)/dx)*dt/dx;    }
    double C4(int i, int j) { return - a4(i,j)*dt; }
    
    // Solves the Partial Differential Equation by using the explicit discretization method.
    void Solve();
};


#endif /* FDMethod_hpp */

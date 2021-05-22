/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 12-05-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Gauss elimination method
extern	void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

// Back-substitution function
extern	void	backsub(Matrix _U, Matrix _d, Matrix _x);

// LU decomposition method
extern	void	LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

// Function that solves LU decomposition method
extern	void	solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix _x);

// Forward-substitution function
extern	void	fwdsub(Matrix _U, Matrix _d, Matrix _x);

// Inverse matrix function using LU decomposition method
extern	void	inv(Matrix _A, Matrix _Ainv);

// Matrix multiplication
extern	Matrix	multiMat(Matrix _Ainv, Matrix _b);

// Gauss Jordan method
extern	void	gaussJordan(Matrix _A, Matrix _b, Matrix _D, Matrix _d);

// QR decomposition method
extern	void	QRdecomp(Matrix _A, Matrix _R, Matrix _Q, Matrix _H);

// Norm function
extern	double	norm(Matrix _A);

// Matrix transpose function
extern	Matrix	trans(Matrix _A);

// Eig function
extern	void	eig(Matrix _A, Matrix _R, Matrix _Q, Matrix _H, Matrix _eVal);

// Newton's method
void fzeroNewton(Matrix _x0, double _tol, Matrix _xn);

// Non-linear function 1
double func1(double _x, double _y);

// Non-linear function 2
double func2(double _x, double _y);

// Derived function of non-linear function 1 by x
double dxfunc1(double _x);

// Derived function of non-linear function 1 by y
double dyfunc1(double _y);

// Derived function of non-linear function 2 by x
double dxfunc2(double _x);

// Derived function of non-linear function 1 by y
double dyfunc2(double _y);

//// Bisection method and Newton-Raphson method ////

// Defined non-linear function
extern double func(double _x);

// The derivative function 
extern double dfunc(double _x);

// The Bisection function
extern double bisectionNL(double _a0, double _b0, double _tol);

// The Newton-Raphson function
extern double newtonRaphson(double _x0, double _tol);

// Defined non-linear function of advanced problem
extern double efunc(double _x);

// The derivative function of advanced problem
extern double defunc(double _x);

// The Bisection function for advanced problem
extern double bisectionNL2(double _a0, double _b0);

// The Hybrid function ( Bisection Method + Newton-Raphson Method )
extern double hybrid(double _x0, double _a0, double _b0, double _tol);

// Condition number function
extern double cond(Matrix _A);

// Linear least square regression function
extern Matrix linearFit(Matrix _x, Matrix _y);

// Linear spline interpolation function
extern Matrix linearInterp(Matrix _x, Matrix _y, Matrix _xq);

// Return the dy/dx results for the input data. (truncation error: O(h^2))
extern Matrix	gradient(Matrix _x, Matrix _y);

// Define a derivative function 
extern double mydFunc(const double _x);

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
extern Matrix	gradientFunc(double func(const double _x), Matrix _xin);

// Return the dy/dx results for the input data. (truncation error: O(h^2))
extern void	gradient1D(double _x[], double _y[], double _dydx[], int _m);

// Modified Newton-Raphson function 
extern double newtonRaphsonFunc(double func(const double _x), double dfunc(const double _x), double _x0, double _tol);

// Rectangle Method function
extern double rect(Matrix _x, Matrix _y);

// Trapezoidal Method function
extern double trapz(double x[], double y[], int m);

// Composite Rectangular Method function
extern double comrect(Matrix _x, Matrix _y);

// Composite Mid-Point Method function
extern double commid(Matrix _x, Matrix _y);

// Composite Trapezoidal Method function
extern double comtrape(Matrix _x, Matrix _y);

// Simpson 1/3 Method
extern double comsimpson13(Matrix _x, Matrix _y);

// Simpson 3/8 Method
extern double comsimpson38(Matrix _x, Matrix _y);

// Integration using rectangular method for discrete data inputs
extern double IntegrateRect(double _x[], double _y[], int _m);

// Integration using Simpson 13 method
extern double integral(double func(const double x), double a, double b, int n);

// Mid-Point Method function
extern double integralMid(double x[], double y[], int m);

// Integration using Simpson 38 method
extern double integral38(double func(const double x), double a, double b, int n);

#endif
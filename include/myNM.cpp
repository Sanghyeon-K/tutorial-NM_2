/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 12-05-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"

void tempFunction(int m) {
	printf("Hedllo");
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

// Gauss elimination method
// input : _A(n x n), _b(n x 1)
// output : _U(n x n), _d(n x 1)
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	double divide = 0;

	if (_A.rows != _A.cols)
	{
		printf("\n*************************************************");
		printf("\n     ERROR!: Matrix A is not a square matrix     ");
		printf("\n*************************************************\n");
		exit(0);
	}

	if (_A.rows != _b.rows)
	{
		printf("\n**************************************************************");
		printf("\n     ERROR!: Matrix A and Vector b cannot be multiplied     ");
		printf("\n**************************************************************\n");
		exit(0);
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_U.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 0; i < _b.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			_d.at[i][j] = _b.at[i][j];
		}
	}

	for (int k = 0; k < _U.rows - 1; k++)
	{
		for (int i = k + 1; i < _U.rows; i++)
		{
			divide = _U.at[i][k] / _U.at[k][k];
			for (int j = k; j < _U.cols; j++)
			{
				_U.at[i][j] = _U.at[i][j] - divide * _U.at[k][j];
			}
			_d.at[i][0] = _d.at[i][0] - divide * _d.at[k][0];
		}
	}
}

// Back-substitution function
// input : _U(n x n), _d(n x 1)
// output : _x(n x 1)
void backsub(Matrix _U, Matrix _d, Matrix _x)
{
	double sum = 0;

	for (int i = 0; i < _d.rows; i++)
	{
		for (int j = 0; j < _d.cols; j++)
		{
			_x.at[i][j] = _d.at[i][j];
		}
	}

	for (int i = _U.cols - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < _U.cols; j++)
		{
			sum += (_U.at[i][j] * _x.at[j][0]);
		}
		_x.at[i][0] = (_x.at[i][0] - sum) / _U.at[i][i];
		sum = 0;
	}
}

// LU decomposition method
// input : _A(n x n)
// output : _L(n x n), _U(n x n), _P(n x n)
void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
{
	Matrix I = eye(_U.rows, _U.cols);
	double maxV;
	double SP;
	double maxSP;
	double divide = 0;
	double temp = 0;
	int maxrow;
	int count;

	//printf("==========================================================================\n");
	//printf("                         LU Decomposition Procedure                       \n");
	//printf("==========================================================================\n\n");

	if (_A.rows != _A.cols)
	{
		printf("\n*************************************************");
		printf("\n     ERROR!: Matrix A is not a square matrix     ");
		printf("\n*************************************************\n");
		exit(0);
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_U.at[i][j] = _A.at[i][j];
		}
	}

	for (int k = 0; k < _U.rows; k++)
	{
		//printf(" At k = %d\n\n", k);

		maxV = 0;
		SP = 0;
		maxSP = 0;
		maxrow = 0;
		count = 0;

		// Scaled partial pivoting
		for (int i = k; i < _U.rows; i++)
		{
			for (int j = k; j < _U.cols; j++)
			{
				if (fabs(_U.at[i][j]) > maxV && _U.at[i][j] != 0)
				{
					maxV = fabs(_U.at[i][j]);
				}
				count++;
			}
			SP = fabs(_U.at[i][k]) / maxV;
			if (SP > maxSP)
			{
				maxSP = SP;
				maxrow = i;
			}
			if (count == _U.cols - k)
			{
				count = 0;
				maxV = 0;
			}
		}

		if (maxrow == k) // No row exchange
		{
			// Apply row reduction
			if (k < _U.rows - 1) 
			{
				for (int i = k + 1; i < _U.rows; i++)
				{
					divide = _U.at[i][k] / _U.at[k][k];
					_L.at[i][k] = divide;
					for (int j = k; j < _U.cols; j++)
					{
						_U.at[i][j] = _U.at[i][j] - divide * _U.at[k][j];
					}
				}
			}
			if (k == _U.rows - 1)
			{
				for (int i = 0; i < _L.rows; i++)
				{
					for (int j = 0; j < _L.cols; j++)
					{
						_L.at[i][j] = _L.at[i][j] + I.at[i][j];
					}
				}
			}
		}
		else // Apply row exchange
		{
			// Row exchange at Matrix U, L, P
			for (int j = 0; j < _U.cols; j++)
			{
				temp = _U.at[k][j];
				_U.at[k][j] = _U.at[maxrow][j];
				_U.at[maxrow][j] = temp;
			}
			for (int j = 0; j < _L.cols; j++)
			{
				temp = _L.at[k][j];
				_L.at[k][j] = _L.at[maxrow][j];
				_L.at[maxrow][j] = temp;
			}
			for (int j = 0; j < _P.cols; j++)
			{
				temp = _P.at[k][j];
				_P.at[k][j] = _P.at[maxrow][j];
				_P.at[maxrow][j] = temp;
			}

			// Apply row reduction
			if (k < _U.rows - 1)
			{
				for (int i = k + 1; i < _U.rows; i++)
				{
					divide = _U.at[i][k] / _U.at[k][k];
					_L.at[i][k] = divide;
					for (int j = k; j < _U.cols; j++)
					{
						_U.at[i][j] = _U.at[i][j] - divide * _U.at[k][j];
					}
				}
			}
			if (k == _U.rows - 1)
			{
				for (int i = 0; i < _L.rows; i++)
				{
					for (int j = 0; j < _L.cols; j++)
					{
						_L.at[i][j] = _L.at[i][j] + I.at[i][j];
					}
				}

			}
		}
		//printMat(_P, " [ matrix P ]");
		//printMat(_L, " [ matrix L ]");
		//printMat(_U, " [ matrix U ]");
	}
}

// Function that solves LU decomposition method
// input : _L(n x n), _U(n x n), _P(n x n), _b(n x 1)
// output : _x(n x 1)
void solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix _x)
{
	Matrix y = zeros(_L.rows, _b.cols);
	Matrix Pb = zeros(_P.rows, _b.cols);

	if (_P.cols != _b.rows)
	{
		printf("\n**************************************************************");
		printf("\n      ERROR!: Matrix P and Vector b cannot be multiplied      ");
		printf("\n**************************************************************\n");
		exit(0);
	}

	for (int i = 0; i < _P.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			Pb.at[i][j] = 0;
			for (int k = 0; k < _P.cols; k++)
			{
				Pb.at[i][j] += _P.at[i][k] * _b.at[k][j];
			}
		}
	}
	fwdsub(_L, Pb, y);
	backsub(_U, y, _x);
}

// Forward-substitution function
// input : _L(n x n), _b(n x 1)
// output : _y(n x 1)
void fwdsub(Matrix _L, Matrix _b, Matrix _y)
{
	double sum = 0;

	for (int i = 0; i < _b.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			_y.at[i][j] = _b.at[i][j];
		}
	}

	for (int i = 0; i < _L.cols; i++)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			sum += (_L.at[i][j] * _y.at[j][0]);
		}
		_y.at[i][0] = (_y.at[i][0] - sum) / _L.at[i][i];
		sum = 0;
	}
}
 
// Inverse matrix function using LU decomposition method
// input : _A(n x n)
// output : _Ainv(n x n)
void inv(Matrix _A, Matrix _Ainv)
{
	Matrix invL = zeros(_A.rows, _A.cols);
	Matrix invU = zeros(_A.rows, _A.cols);
	Matrix invP = eye(_A.rows, _A.cols);
	Matrix invb = eye(_A.rows, _A.cols);
	Matrix invy = zeros(_A.rows, 1);
	Matrix invx = zeros(_A.rows, 1);
	double fr = 1;

	//printf("==========================================================================\n");
	//printf("            Inverse Matrix Procedure [ Using LU Decomposition ]           \n");
	//printf("==========================================================================\n");

	if (_A.rows != _A.cols)
	{
		printf("\n*************************************************");
		printf("\n     ERROR!: Matrix A is not a square matrix     ");
		printf("\n*************************************************\n");
		exit(0);
	}

	LUdecomp(_A, invL, invU, invP);

	for (int i = 0; i < invU.rows; i++)
	{
		for (int j = 0; j < invU.cols; j++)
		{
			if (i == j)
			{
				fr *= invU.at[i][j];
			}
		}
	}

	if (fr == 0)
	{
		printf("\n*************************************************");
		printf("\n   ERROR!: Matrix A is not a full rank matrix    ");
		printf("\n*************************************************\n");
		exit(0);
	}

	for (int j = 0; j < _A.cols; j++)
	{
		for (int i = 0; i < _A.rows; i++)
		{
			invy.at[i][0] = invb.at[i][j];
		}
		solveLU(invL, invU, invP, invy, invx);
		for (int i = 0; i < _A.rows; i++)
		{
			_Ainv.at[i][j] = invx.at[i][0];
		}
	}
}

// Matrix multiplication
// input : _Ainv(a x b), _b(b x c)
// output : _x(a x c)
Matrix multiMat(Matrix _Ainv, Matrix _b)
{
	Matrix x = zeros(_Ainv.rows, _b.cols);

	for (int i = 0; i < _Ainv.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			for (int k = 0; k < _Ainv.cols; k++)
			{
				x.at[i][j] += _Ainv.at[i][k] * _b.at[k][j];
			}
		}
	}

	return x;
}

// Gauss Jordan method
// input : _A(n x n), _b(n x 1)
// output : _D(n x n), _d(n x 1)
void gaussJordan(Matrix _A, Matrix _b, Matrix _D, Matrix _d)
{
	double maxV;
	double SP;
	double maxSP;
	double temp = 0;
	double divide1 = 0;
	double divide2 = 0;
	int maxrow;
	int count;

	if (_A.rows != _A.cols)
	{
		printf("\n*************************************************");
		printf("\n     ERROR!: Matrix A is not a square matrix     ");
		printf("\n*************************************************\n");
		exit(0);
	}

	if (_A.rows != _b.rows)
	{
		printf("\n**************************************************************");
		printf("\n     ERROR!: Matrix A and Vector b cannot be multiplied     ");
		printf("\n**************************************************************\n");
		exit(0);
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_D.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 0; i < _b.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			_d.at[i][j] = _b.at[i][j];
		}
	}

	for (int k = 0; k < _D.rows; k++)
	{
		maxV = 0;
		SP = 0;
		maxSP = 0;
		maxrow = 0;
		count = 0;

		// Scaled partial pivoting
		for (int i = k; i < _D.rows; i++)
		{
			for (int j = k; j < _D.cols; j++)
			{
				if (fabs(_D.at[i][j]) > maxV && _D.at[i][j] != 0)
				{
					maxV = fabs(_D.at[i][j]);
				}
				count++;
			}
			SP = fabs(_D.at[i][k]) / maxV;
			if (SP > maxSP)
			{
				maxSP = SP;
				maxrow = i;
			}
			if (count == _D.cols - k)
			{
				count = 0;
				maxV = 0;
			}
		}

		if (maxrow == k)
		{
			divide1 = _D.at[k][k];

			for (int j = k; j < _D.cols; j++)
			{
				_D.at[k][j] = _D.at[k][j] / divide1;
			}
			_d.at[k][0] = _d.at[k][0] / divide1;

			for (int i = 0; i < _D.rows; i++)
			{
				if (i != k)
				{
					divide2 = _D.at[i][k];

					for (int j = k; j < _D.cols; j++)
					{
						_D.at[i][j] = _D.at[i][j] - divide2 * _D.at[k][j];
					}
					_d.at[i][0] = _d.at[i][0] - divide2 * _d.at[k][0];
				}
			}
		}
		else
		{
			for (int j = 0; j < _D.cols; j++)
			{
				temp = _D.at[k][j];
				_D.at[k][j] = _D.at[maxrow][j];
				_D.at[maxrow][j] = temp;
			}

			for (int j = 0; j < _d.cols; j++)
			{
				temp = _d.at[k][j];
				_d.at[k][j] = _d.at[maxrow][j];
				_d.at[maxrow][j] = temp;
			}

			divide1 = _D.at[k][k];

			for (int j = k; j < _D.cols; j++)
			{
				_D.at[k][j] = _D.at[k][j] / divide1;
			}
			_d.at[k][0] = _d.at[k][0] / divide1;

			for (int i = 0; i < _D.rows; i++)
			{
				if (i != k)
				{
					divide2 = _D.at[i][k];

					for (int j = k; j < _D.cols; j++)
					{
						_D.at[i][j] = _D.at[i][j] - divide2 * _D.at[k][j];
					}
					_d.at[i][0] = _d.at[i][0] - divide2 * _d.at[k][0];
				}
			}
		}
	}
}

// QR decomposition method
// input : _A(n x n)
// output : _R(n x n), _Q(n x n), _H(n x n)
void QRdecomp(Matrix _A, Matrix _R, Matrix _Q, Matrix _H)
{
	Matrix c = zeros(_A.rows, 1);
	Matrix e = zeros(_A.rows, 1);
	Matrix v = zeros(_A.rows, 1);
	Matrix vt = zeros(1, _A.rows);
	Matrix vvt = zeros(_A.rows, _A.rows);
	Matrix vtv = zeros(1, 1);
	Matrix I = eye(_A.rows, _A.rows);
	Matrix QH = zeros(_A.rows, _A.rows);
	Matrix HR = zeros(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_R.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 0; i < _Q.rows; i++)
	{
		for (int j = 0; j < _Q.cols; j++)
		{
			if (j == i)
			{
				_Q.at[i][j] = 1;
			}
		}
	}

	for (int k = 0; k < _R.rows - 1; k++)
	{
		for (int i = 0; i < _R.rows; i++)
		{
			c.at[i][0] = _R.at[i][k];
		}
		
		for (int i = 0; i < k; i++)
		{
			c.at[i][0] = 0;
		}

		for (int i = 0; i < _R.rows; i++)
		{
			if (i == k)
			{
				if (_R.at[i][i] >= 0)
				{
					e.at[i][0] = 1;
				}
				else
				{
					e.at[i][0] = -1;
				}
			}
			else
			{
				e.at[i][0] = 0;
			}
		}

		for (int i = 0; i < _A.rows; i++)
		{
			v.at[i][0] = c.at[i][0] + norm(c) * e.at[i][0];
		}

		vt = trans(v);

		vvt = multiMat(v, vt);

		vtv = multiMat(vt, v);

		for (int i = 0; i < _H.rows; i++)
		{
			for (int j = 0; j < _H.cols; j++)
			{
				_H.at[i][j] = I.at[i][j] - ( 2 * vvt.at[i][j]) / vtv.at[0][0];
			}
		}
		
		QH = multiMat(_Q, _H);

		for (int i = 0; i < _Q.rows; i++)
		{
			for (int j = 0; j < _Q.cols; j++)
			{
				_Q.at[i][j] = QH.at[i][j];
			}
		}

		HR = multiMat(_H, _R);

		for (int i = 0; i < _R.rows; i++)
		{
			for (int j = 0; j < _R.cols; j++)
			{
				_R.at[i][j] = HR.at[i][j];
			}
		}
	}
}

// Norm function
// input : _A(n x 1)
// output : sqrt(sum)
double norm(Matrix _A)
{
	double sum = 0;
	for (int i = 0; i < _A.rows; i++)
	{
		sum += pow(_A.at[i][0], 2);
	}

	return sqrt(sum);
}

// Matrix transpose function
// input : _A(n x 1)
// output : _TM(1 x n)
Matrix trans(Matrix _A)
{
	double temp = 0;
	Matrix TM = zeros(_A.cols, _A.rows);

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			TM.at[j][i] = _A.at[i][j];
		}
	}

	return(TM);
}

// Eig function
// input : _A(n x n), _R(n x n), _Q(n x n), _H(n x n)
// output : _eVal(n x 1)
void eig(Matrix _A, Matrix _R, Matrix _Q, Matrix _H, Matrix _eVal)
{
	int Nmax = 100;

	Matrix A = zeros(_A.rows, _A.cols);

	if (_A.rows != _A.cols)
	{
		printf("\n*************************************************");
		printf("\n     ERROR!: Matrix A is not a square matrix     ");
		printf("\n*************************************************\n");
		exit(0);
	}
	
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			A.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 1; i <= Nmax; i++)
	{
		QRdecomp(A, _R, _Q, _H);

		A = multiMat(_R, _Q);

		initMat(_R, 0);
		initMat(_Q, 0);
		initMat(_H, 0);
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			if (i == j)
			{
				_eVal.at[i][0] = A.at[i][j];
			}
		}
	}
}

double func1(double _x, double _y)
{
	return (_y - (exp(_x / 2) + exp(-_x / 2)) / 2);
}

double func2(double _x, double _y)
{
	return (9 * pow(_x, 2) + 25 * pow(_y, 2) - 225);
}

double dxfunc1(double _x)
{
	return (-(exp(_x / 2) - exp(-_x / 2)) / 4);
}

double dyfunc1(double _y)
{
	return 1;
}

double dxfunc2(double _x)
{
	return 18 * _x;
}

double dyfunc2(double _y)
{
	return 50 * _y;
}

// Newton's method
// input : _x0(n x 1), double tol
// output : _xn(n x 1)
void fzeroNewton(Matrix _x0, double _tol, Matrix _xn)
{
	Matrix F = zeros(_x0.rows, 1);
	Matrix J = zeros(_x0.rows, _x0.rows);
	Matrix D = zeros(_x0.rows, _x0.rows);
	Matrix x = zeros(_x0.rows, 1);
	double ep1 = 1;
	double ep2 = 1;
	int Nmax = 100;
	int k = 0;

	for (int i = 0; i < _x0.rows; i++)
	{
		for (int j = 0; j < _x0.cols; j++)
		{
			_xn.at[i][j] = _x0.at[i][j];
		}
	}

	do
	{
		printf(" [ k = %d ] \n\n", k);

		F.at[0][0] = (-1) * func1(_xn.at[0][0], _xn.at[1][0]);
		F.at[1][0] = (-1) * func2(_xn.at[0][0], _xn.at[1][0]);

		printMat(_xn, " xn ");
		printMat(F, " F ");

		J.at[0][0] = dxfunc1(_xn.at[0][0]);
		J.at[0][1] = dyfunc1(_xn.at[1][0]);
		J.at[1][0] = dxfunc2(_xn.at[0][0]);
		J.at[1][1] = dyfunc2(_xn.at[1][0]);

		printMat(J, " J ");

		gaussJordan(J, F, D, x);

		printMat(x, " x ");
		_xn.at[0][0] += x.at[0][0];
		_xn.at[1][0] += x.at[1][0];
		printMat(_xn, " xn ");

		ep1 = fabs(x.at[0][0]);
		ep2 = fabs(x.at[1][0]);

		initMat(D, 0);

		k++;
	} while (k < Nmax && ep1 > _tol && ep2 > _tol);
}

//// Bisection method and Newton-Raphson method ////

// Defined non-linear function
double func(double _T)
{
	double a = 2.624 * pow(10, 4);
	double b = -0.03734;
	double R = 800;
	double y;

	y = a * exp(b * _T) - R;
	
	return y;
}

// The derivative function
double dfunc(double _T)
{
	double a = 2.624 * pow(10, 4);
	double b = -0.03734;
	double dy;

	dy = a * b * exp(b * _T);

	return dy;
}

// The Bisection function
double bisectionNL(double _a0, double _b0, double _tol)
{
	int k = 0;
	int Nmax = 1000;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;

	// Error checking : if func(a0) * func(b0) > 0
	if (func(a) * func(b) > 0)
	{
		xn = 0;
		printf("\n******* Error : func(%f) and func(%f) have same sign *******\n\n", a, b);
		return xn;
	}

	// Checking : if func(a0) = 0 or func(b0) = 0
	if (func(a) == 0 || func(b) == 0)
	{
		if (func(a) == 0)
		{
			xn = a;
			printf("\n******* Solution of this non-linear function is %f *******\n\n", a);
			return xn;
		}
		else if (func(b) == 0)
		{
			xn = b;
			printf("\n******* Solution of this non-linear function is %f *******\n\n", b);
			return xn;
		}
	}

	do {
		xn = (a + b) / 2;
		// Error checking : if dfunc(xn)=0
		if (dfunc(xn) == 0)
		{
			printf("\n******* Error : Function is not continous at X(n) = %f  *******\n\n", xn);
			xn = 0;
			break;
		}
		ep = fabs(func(xn));
		printf("Iteration : %d \t", k);
		printf("X(n) : %f \t", xn);
		printf("Tolerance : %.10f\n", ep);
		if (func(a) * func(xn) < 0)
			b = xn;
		else
			a = xn;
		k++;
	} while (k<Nmax && ep>_tol);

	return xn;
}


// The Newton-Raphson function
double newtonRaphson(double _T0, double _tol)
{
	int k = 0;
	int Nmax = 1000;
	double T0 = _T0;
	double T = 0;
	double ep = 1000;
	double hk = 0;

	do {
		hk = -func(T0) / dfunc(T0);
		T = T0 + hk;
		// Error checking : if dfunc(xn)=0
		if (dfunc(T) == 0)
		{
			printf("\n******* Error : Function is not continous at T(K) = %f  *******\n\n", T);
			T = 0;
			break;
		}
		ep = fabs(T - T0);
		printf("Iteration : %d \t", k);
		printf("T(k) : %f \t", T);
		printf("Tolerance : %.10f\n", ep);
		T0 = T;
		k++;
	} while (k<Nmax && ep>_tol);

	return T;
}

// Defined non-linear function of advanced problem
double efunc(double _x)
{
	double y;

	y = pow(_x, -1) - 2;

	return y;
}

// The derivative function of advanced problem
double defunc(double _x)
{
	double dy;

	dy = -pow(_x, -2);

	return dy;
}

// The Bisection function for advanced problem
double bisectionNL2(double _a0, double _b0)
{
	double a = _a0;
	double b = _b0;
	double xn = 0;

	xn = (a + b) / 2;
	if (efunc(a) * efunc(xn) < 0)
		b = xn;
	else
		a = xn;

	return xn;
}

// The Hybrid function ( Bisection Method + Newton-Raphson Method )
double hybrid(double _x0, double _a0, double _b0, double _tol)
{
	int k = 0;
	int Nmax = 1000;
	double x0 = _x0;
	double x = 0;
	double ep = 1000;
	double hk = 0;

	do {
		hk = -efunc(x0) / defunc(x0);
		x = x0 + hk;
		// Error checking : if defunc(xn)=0
		if (defunc(x) == 0)
		{
			printf("\n******* Error : Function is not continous at X(n) = %f  *******\n\n", x);
			x = 0;
			break;
		}
		// Using Bisection method if Newton-Raphson method gives the solution out of bound [ a0, b0 ]
		// Otherwise use Newton-Raphson method
		if (x < _a0 || x > _b0)
		{
			if (x < _a0)
			{
				printf("\n* Newton-Raphson method gives X(n) = %lf which is less than %f *\n", x, _a0);
				printf("[ Using Bisection method ]\n");
				x = bisectionNL2(x, _b0);
			}
			else if (x > _b0)
			{
				printf("\n* Newton-Raphson method gives X(n) = %lf which is less than %f *\n", x, _b0);
				printf("[ Using Bisection method ]\n");
				x = bisectionNL2(_a0, x);
			}
		}
		else
		{
			printf("\n[ Using Newton-Raphson method ]\n");
		}
		ep = fabs(x - x0);
		printf("Iteration : %d \t", k);
		printf("X(n) : %f \t", x);
		printf("Tolerance : %.10f\n", ep);
		x0 = x;
		k++;
	} while (k<Nmax && ep>_tol);

	return x;
}

// Condition number function
// input : _A(n x n)
// output : sqrt(maxEV / minEV)
double cond(Matrix _A)
{
	Matrix R = zeros(_A.cols, _A.cols);
	Matrix Q = zeros(_A.cols, _A.cols);
	Matrix H = zeros(_A.cols, _A.cols);
	Matrix eVal = zeros(_A.cols, 1);
	Matrix MatT = zeros(_A.cols, _A.rows);
	Matrix MatM = zeros(_A.cols, _A.cols);
	double maxEV = 0;
	double minEV = 10000000000000000;

	MatT = trans(_A);

	MatM = multiMat(MatT, _A);

	eig(MatM, R, Q, H, eVal);
	
	for (int i = 0; i < eVal.rows; i++)
	{
		if (fabs((eVal.at[i][0])) > maxEV)
		{
			maxEV = fabs((eVal.at[i][0]));
		}
		if (fabs(eVal.at[i][0]) < minEV && eVal.at[i][0] != 0)
		{
			minEV = fabs(eVal.at[i][0]);
		}
	}

	return sqrt(maxEV / minEV);
}

// Linear least square regression function
// input : _x(n x 1), _y(n x 1)
// output : z(2 x 1)
Matrix linearFit(Matrix _x, Matrix _y)
{
	int nx = 0;
	int ny = 0;
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;
	double a0 = 0;
	double a1 = 0;

	nx = _x.rows;
	ny = _y.rows;

	if ((nx != ny) || (nx < 2))
	{
		printf(" ERROR: Length of x oand y must be equal and should be at least 2. \n\n");
		exit(0);
	}

	for (int i = 0; i < _x.rows; i++)
	{
		Sx += _x.at[i][0];
		Sy += _y.at[i][0];
		Sxy += _x.at[i][0] * _y.at[i][0];
		Sxx += pow(_x.at[i][0], 2);
	}

	a0 = (Sxx * Sy - Sxy * Sx) / (nx * Sxx - pow(Sx, 2));
	a1 = (nx * Sxy - Sx * Sy) / (nx * Sxx - pow(Sx, 2));

	double z_array[] = { a0, a1 };
	return arr2Mat(z_array, 2, 1);
}

// Linear spline interpolation function
// input : _x(n x 1), _y(n x 1), _xq(m x 1)
// output : yq(m x 1)
Matrix linearInterp(Matrix _x, Matrix _y, Matrix _xq)
{
	int k = 0;
	int nx = 0;
	int ny = 0;
	Matrix y = zeros(_xq.rows, 1);

	nx = _x.rows;
	ny = _y.rows;

	if ((nx != ny) || (nx < 2))
	{
		printf(" ERROR: Length of x oand y must be equal and should be at least 2. \n\n");
		exit(0);
	}

	for (int i = 0; i < _xq.rows; i++)
	{
		if (i == 0)
		{
			y.at[i][0] = _y.at[i][0];
		}
		else
		{
			for (int j = 1; j < _x.rows; j++)
			{
				if (_xq.at[i][0] <= _x.at[j][0])
				{
					k = j;
					break;
				}
			}
			y.at[i][0] = _y.at[k - 1][0] * ((_xq.at[i][0] - _x.at[k][0]) / (_x.at[k - 1][0] - _x.at[k][0])) + _y.at[k][0] * ((_xq.at[i][0] - _x.at[k - 1][0]) / (_x.at[k][0] - _x.at[k - 1][0]));
		}
	}

	return y;
}

// Return the dy/dx results for the input data. (truncation error: O(h^2))
// input : _x(n x 1), _y(n x 1)
// output : dy(n x 1)
Matrix	gradient(Matrix _x, Matrix _y) {

	int nx = _x.rows;
	int ny = _y.rows;

	Matrix dy = zeros(nx, 1);

	if ((nx != ny) || (nx < 2))
	{
		printf(" ERROR: Length of x oand y must be equal and should be at least 2. \n\n");
		exit(0);
	}

	for (int i = 0; i < nx; i++)
	{
		if (i == 0)
		{
			if (nx == 2)
			{
				dy.at[i][0] = (_y.at[i + 1][0] - _y.at[i][0]) / (_x.at[i + 1][0] - _x.at[i][0]);
			}
			else
			{
				dy.at[i][0] = (-3 * _y.at[i][0] + 4 * _y.at[i + 1][0] - _y.at[i + 2][0]) / (2 * (_x.at[i + 1][0] - _x.at[i][0]));
			}
		}
		else if ((i > 0) && (i < nx - 1))
		{
			dy.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * (_x.at[i + 1][0] - _x.at[i][0]));
		}
		else
		{
			if (nx == 2)
			{
				dy.at[i][0] = (_y.at[i][0] - _y.at[i-1][0]) / (_x.at[i][0] - _x.at[i-1][0]);
			}
			else
			{
				dy.at[i][0] = (_y.at[i - 2][0] - 4 * _y.at[i - 1][0] + 3 * _y.at[i][0]) / (2 * (_x.at[i][0] - _x.at[i - 1][0]));
			}
		}
	}

	return dy;
}

// Return the dy/dx results for the input data. (truncation error: O(h^2))
// input : _x(m), _y(m), m(int)
// output : _dydx(m)
void gradient1D(double _x[], double _y[], double _dydx[], int _m)
{
	if (_m < 2)
	{
		printf(" ERROR: Length of x oand y should be at least 2. \n\n");
		exit(0);
	}

	for (int i = 0; i < _m; i++)
	{
		if (i == 0)
		{
			if (_m == 2)
			{
				_dydx[i] = (_y[i + 1] - _y[i]) / (_x[i + 1] - _x[i]);
			}
			else
			{
				_dydx[i] = (-3 * _y[i] + 4 * _y[i + 1] - _y[i + 2]) / (2 * (_x[i + 1] - _x[i]));
			}
		}
		else if ((i > 0) && (i < _m - 1))
		{
			_dydx[i] = (_y[i + 1] - _y[i - 1]) / (2 * (_x[i + 1] - _x[i]));
		}
		else
		{
			if (_m == 2)
			{
				_dydx[i] = (_y[i] - _y[i - 1]) / (_x[i] - _x[i - 1]);
			}
			else
			{
				_dydx[i] = (_y[i - 2] - 4 * _y[i - 1] + 3 * _y[i]) / (2 * (_x[i] - _x[i - 1]));
			}
		}
	}
}

// Define a function that defines the target equation.
// input : _x(double)
// output : _x^3(double)
double myFunc(const double _x) {
	return  _x * _x * _x;
}

// Define a derivative function
// input : _x(double)
// output : 3*_x^2(double)
double mydFunc(const double _x) {
	return 3 * _x * _x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// input : func(double), _xin(n x 1)
// output : dy(n x 1)
Matrix	gradientFunc(double func(const double _x), Matrix _xin) {

	int n = _xin.rows;

	if (n < 2)
	{
		printf(" ERROR: Length of xin should be at least 2. \n\n");
		exit(0);
	}

	Matrix y = zeros(n, 1);

	for (int i = 0; i < n; i++)
	{
		y.at[i][0] = func(_xin.at[i][0]);
	}

	Matrix dy = gradient(_xin, y);

	return dy;
}

// Modified Newton-Raphson function
// input : func(double), dfunc(double), _x0(double), _tol(double)
// output : x(double)
double newtonRaphsonFunc(double func(const double _x), double dfunc(const double _x), double _x0, double _tol)
{
	int k = 0;
	int Nmax = 1000;
	double x0 = _x0;
	double x = 0;
	double ep = 1000;
	double hk = 0;

	do {
		hk = -func(x0) / dfunc(x0);
		x = x0 + hk;
		// Error checking : if dfunc(xn)=0
		if (dfunc(x) == 0)
		{
			printf("\n******* Error : Function is not continous at X(K) = %f  *******\n\n", x);
			x = 0;
			break;
		}
		ep = fabs(x - x0);
		printf("Iteration : %d \t", k);
		printf("X(k) : %f \t", x);
		printf("Tolerance : %.10f\n", ep);
		x0 = x;
		k++;
	} while (k<Nmax && ep>_tol);

	return x;
}

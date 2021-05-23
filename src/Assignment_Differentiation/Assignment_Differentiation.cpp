/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 12-05-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	6		// enter your assignment number
#define eval		0		// set 0

//#include "myNM.h"

#include "../../include/myNM.h"

int main(int argc, char* argv[])
{
	// MODIFIED
	// HELLO ±³¼ö´Ô

	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/

	double t_a[21] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0 };
	double pos_a[21] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
	double vel_a[21] = { 0 };
	double acc_a[21] = { 0 };
	int m = 21;

	Matrix t = arr2Mat(t_a, m, 1);
	Matrix pos = arr2Mat(pos_a, m, 1);
	Matrix xin = arr2Mat(t_a, m, 1);

	double x0 = 2.0;
	double tol = 0.0000001;

	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);

	gradient1D(t_a, pos_a, vel_a, m);
	gradient1D(t_a, vel_a, acc_a, m);

	Matrix dydx = gradientFunc(myFunc, xin);

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	
	// PART 1
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	printMat(t, " t ");
	printMat(pos, " pos ");
	printMat(vel, " vel ");
	printMat(acc, " acc ");

	// PART 1
	printf("\n**************************************************");
	printf("\n|           PART 1. (using 1D arrays)            |");
	printf("\n**************************************************\n");

	printf(" vel-1D = \n\n");
	for (int i = 0; i < m; i++)
	{
		printf("%15.6f\n", vel_a[i]);
	}

	printf("\n\n acc-1D = \n\n");
	for (int i = 0; i < m; i++)
	{
		printf("%15.6f\n", acc_a[i]);
	}

	// PART 2
	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	printMat(xin, " xin ");
	printMat(dydx, " dydx ");

	// PART 2
	printf("\n**************************************************");
	printf("\n|        PART 2. (Newton-Raphson function)       |");
	printf("\n**************************************************\n");

	double NR_result = newtonRaphsonFunc(myFunc, mydFunc, x0, tol);
	printf("Final solution of Newton-Raphson method : %lf\n\n", NR_result);

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	freeMat(t);
	freeMat(pos);
	freeMat(vel);
	freeMat(acc);
	freeMat(xin);
	freeMat(dydx);

	system("pause");
	return 0;
}
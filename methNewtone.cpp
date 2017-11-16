// methNewtone.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "iostream"
#include "fstream"
#include <iomanip>

using namespace std;

void findAndSepRoots(int n, double* koefs, double eps);
double* solveSqrEq(int n, double* koefs);
bool NewtConditions(int n, double* koefs, double* segment);
bool intersects(double x, double* segment);
double func(int n, double *koefs, double x);
void showFunc(int n, double* koefs);
void methNewtone(int n, double* koefs, double* segment, double eps);


int _tmain(int argc, _TCHAR* argv[])
{
	setlocale(LC_ALL, "Russian");
	ifstream fin("input.txt");

	//порядок уравнения
	int n;
	fin >> n;

	//коэффиценты уравнения в порядке убывния степени, включая нули и перед нулевой степенью
	double *koefs = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
		fin >> koefs[i];

    fin.close();

	double eps;

	showFunc(n, koefs);
	do{
		cout << "Введите точность: ";
		cin >> eps;
	} while (eps <= 0);

	findAndSepRoots(n, koefs, eps);
	return 0;
}

double* diff(int n, double* koefs) {
	int k;
	if (n == 0)
		k = n;
	else
		k = n - 1;

    double* koefDiff = new double[k+1];

	if (n > 0) {
		for (int i = 0; i < n; i++)
			koefDiff[i] = koefs[i] * (n - i);
    }
	else{
		koefDiff[0] = 0;
	}

	return koefDiff;
}

void findAndSepRoots(int n, double* koefs, double eps) {
	double** segments = NULL;
    double* koefDiff = diff(n, koefs);

	double A = abs(koefs[1]);
	for (int i = 2; i < n + 1; i++){
		if (A < abs(koefs[i]))
			A = abs(koefs[i]);
	}

	double B = abs(koefs[0]);
	for (int i = 1; i < n; i++){
		if (B < abs(koefs[i]))
			B = abs(koefs[i]);
	}

	//тогда модули всех корней строго больше r и нестрого меньше R
	double r = 1 / (1 + B / abs(koefs[n]));
	double R = 1 + A / abs(koefs[0]);

	cout << "Модули всех корней лежат между " << r << " и " << R << endl;
	int k = 2;
	segments = new double*[k];
	for (int i = 0; i < n; i++)
		segments[i] = new double[k];

	segments[0][0] = (-1) * (R + eps);
	segments[0][1] = (-1) * r;

	segments[1][0] = r;
	segments[1][1] = R + eps;

	for (int i = 0; i < k; i++) {
		if (func(n, koefs, segments[i][0]) * func(n, koefs, segments[i][1]) <= 0){
			cout << "[" << segments[i][0] << " ; " << segments[i][1] << "], достаточное условие";
			if (NewtConditions(n, koefs, segments[i])) {
				cout << " - выполнено\n";
			}
			else{
				cout << " - не выполнено\n";
			}

			methNewtone(n, koefs, segments[i], eps);
		}
	}
}

double* solveSqrEq(int n, double* koefs) {
	if (n == 2) {
		double a = 1;
		double b = koefs[1] / koefs[0];
		double c = koefs[2] / koefs[0];

		double disc = b*b - 4 * a*c;
		double* roots = new double[n];

		if (disc >= 0) {
			roots[0] = ((-1)* b - pow(disc, 0.5)) / (2 * a);
			roots[1] = ((-1)* b + pow(disc, 0.5)) / (2 * a);
			return roots;
		}
        else {
			return NULL;
		}
	}
	else{
		cout << "Данное уравнение не является квадратным\n";
		return NULL;
	}
}

bool NewtConditions(int n, double* koefs, double* segment) {
	if (n == 3) {
		double* koefDiff = diff(n, koefs);
	    int k = n - 1;

		double* rootsDiff = solveSqrEq(k, koefDiff);
		if (rootsDiff != NULL) {
			//если на отрезке есть нули производной
			if (intersects(rootsDiff[0], segment) || intersects(rootsDiff[1], segment))
				return false;
		}

		//вторая производная у кубической функции - прямая
		double* koefSecDiff = diff(k, koefDiff);
		double rootSecDiff;

		if (koefSecDiff[1] == 0)
			rootSecDiff = koefSecDiff[1];
		else
			rootSecDiff = (-1)*(koefSecDiff[1] / koefSecDiff[0]);

		if (intersects(rootSecDiff, segment))
			return false;

		return true;
	}
 
	cout << "\n Данное уравнение не является кубическим, поэтому непроверяемы условия сходимости/n";
	return false;
}

bool intersects(double x, double* segment) {
	if (x < segment[0] || x > segment[1])
		return false;
	else
		return true;
}

double func(int n, double *koefs, double x) {
	double result = 0;
	for (int i = n; i >= 0; i--) {
		result += koefs[n - i] * pow(x, i);
	}

	return result;
}

void showFunc(int n, double* koefs) {
	for (int i = 0; i < n+1; i++) {
		if (koefs[i] != 0) {
			if (koefs[i] > 0 && i > 0)
				cout << "+";

			if ((abs(koefs[i]) != 1 && i < n) || (i==n))
				cout << koefs[i];
			
			if (koefs[i] == -1 && i < n)
				cout << "-";

			if (i < n-1)
			  cout << "x^" << (n - i);

			if (i == n-1)
				cout << "x";
        }
	}
	cout  << " = 0 " << endl;
}

void methNewtone(int n, double* koefs, double* segment, double eps) {
	double* koefDiff = diff(n, koefs);
	double* koefSecDiff = diff(n-1, koefDiff);

	double xCur;
	double xNext;
	//т.к. первая и вторая производная на отрезке сохраняют знак, то проверяем в точке из отрезка
	if (func(n - 1, koefDiff, segment[0]) * func(n-2, koefSecDiff, segment[0]) > 0) {
		xNext = segment[1];
	}
	else{
		xNext = segment[0];
	}

    double length;

	do{
		xCur = xNext;
		xNext = xCur - func(n, koefs, xCur)/func(n - 1, koefDiff, xCur);
		length = abs(xNext - xCur);

		cout << "xK=" << setprecision(8) << xCur;
		cout << " x(K+1)=" << setprecision(8) << xNext;
		cout << " dist(xK, x(K+1))=" << setprecision(8) << length << endl;
	} while (length >= eps);

	cout << "Корень лежит рядом с " << setprecision(8) << xNext << endl;
}

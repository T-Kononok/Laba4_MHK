#include "pch.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

using namespace std;

int readf(ifstream &, string &, string &);

bool GausSh(int *, double **, int, int);

bool Gaus(int *, double *, double **, int);

int Max(int *, double **, int, int);

void swap(int *, int, int);

int main()
{
	setlocale(LC_ALL, "Russian");

	ifstream file("file.txt");

	if (!file)
	{
		cout << "Иди создавай файл!" << endl;
		return 0;
	}
//1
	string sx, sy;
	int n;
	n = readf(file, sx, sy);
	cout << n << "\n";

	double *x = new double[n];
	double *y = new double[n];

	for (int i = 0; i < n; i++)
	{
		istringstream strx(sx);
		istringstream stry(sy);
		strx >> x[i];
		stry >> y[i];
		sx.erase(0, sx.find(' ', 0) + 1);
		sy.erase(0, sy.find(' ', 0) + 1);
	}
//вывод 
	for (int i = 0; i < n; i++)
		cout << "  x"<<i<<"= " << left << setw(5) << x[i];
	cout << endl;
	for (int i = 0; i < n; i++)
		cout << "  y" << i << "= " << left << setw(5) << y[i];
//2
	int m = 0;
	cout << "\nm (m < n) = ";
	cin >> m;

//--------
	double Ex = 0, Ey = 0, Ex2 = 0, Exy = 0, aa, bb;
	for (int i = 0; i < n; i++)
	{
		Ex += x[i];
		Ey += y[i];
		Ex2 += x[i] * x[i];
		Exy += x[i] * y[i];
	}
	aa = (n * Exy - Ex * Ey) / (n * Ex2 - Ex * Ex);
	bb = (Ey - aa * Ex) / n;
	cout << setprecision(20) << fixed << aa << " " << bb << endl;
	cout << setprecision(5) << fixed;

//3
	double *powerx = new double[2 * m + 1];
	double promx = 0;
	for (int k = 0; k < 2 * m + 1; k++)
	{
		promx = 0;
		for (int i = 0; i < n; i++)
			promx += pow(x[i], k);
		powerx[k] = promx;
	}
//4
	double **sumx = new double *[m + 1];
	for (int i = 0; i < m + 1; i++)
		sumx[i] = new double[m + 2];

	for (int i = 0; i < m + 1; i++)
		for (int j = 0; j < m + 1; j++)
			sumx[i][j] = powerx[i + j];
	sumx[0][0] = n;
//5
	for (int k = 0; k < m + 1; k++)
	{
		promx = 0;
		for (int i = 0; i < n; i++)
			promx += pow(x[i], k) * y[i];
		sumx[k][m + 1] = promx;
	}
//вывод
	cout << endl;
	for (int i = 0; i < m+1; i++)
	{
		for (int j = 0; j < m + 1; j++)
			cout << left << setw(10) << sumx[i][j];
		cout << " = " << sumx[i][m+1] << endl << endl;
	}
	cout << endl;
//6
	int *p = new int[m + 1];
	for (int i = 0; i < m+1; i++)
		p[i] = i;

	double *a = new double[m + 1];

	if (!Gaus(p, a, sumx, m + 1))
		return 0;
//7
	double s = 1/(n - m), us = 0, ks = 0;
	for (int i = 0; i < n; i++)
	{
		ks = y[i];
		for (int j = 0; j < m + 1; j++)
			ks -= a[j] * pow(x[i],j);
		ks *= ks;
		us += ks;
	}

	s *= us;
	s *= s;
	s = pow(s, 0.5);
	cout << setprecision(50) << fixed << "q = " << s << endl;

//удаление
	delete[] p;
	delete[] a;
	delete[] x;
	delete[] y;
	delete[] powerx;
	for (int i = 0; i < m + 1; i++)
		delete[] sumx[i];
	delete[] sumx;
}

int readf(ifstream &file, string &sx, string &sy)
{
	int k = 0;
	getline(file, sx);
	k = 0;
	while (k < sx.length()-1)
	{
		if (sx.at(k) == ' ' && sx.at(k + 1) == ' ')
		{
			sx.erase(k, 1);
			k--;
		}
		k++;
	}
	if (sx.at(0) == ' ')
		sx.erase(0, 1);
	if (sx.at(sx.length()-1) == ' ')
		sx.erase(sx.length()-1, 1);

	int n = 0;
	for (int i = 0; i < sx.length(); i++)
		if (sx.at(i) == ' ')
			n++;

	getline(file, sy);
	for (int i = 0; i < n; i++)
		sy += " 0";
	k = 0;
	while (k < sy.length() - 1)
	{
		if (sy.at(k) == ' ' && sy.at(k + 1) == ' ')
		{
			sy.erase(k, 1);
			k--;
		}
		k++;
	}
	if (sy.at(0) == ' ')
		sy.erase(0, 1);

	return n+1;
}


bool Gaus(int *p, double *x, double **a, int n)
{
	int i, j;
	for (int k = 0; k < n; k++)
	{
		if (GausSh(p, a, n, k))
		{
		}
		else
			return 0;
	}

	double sumx;
	for (i = n - 1; i >= 0; i--)
	{
		sumx = 0;
		for (j = i + 1; j < n; j++)
			sumx += a[p[i]][j] * x[j];
		x[i] = a[p[i]][n] - sumx;
	}

	for (i = 0; i < n; i++)
		cout << "a" << i << " = " << setprecision(20) << fixed << x[i] << endl;

	return 1;
}

bool GausSh(int *p, double **a, int n, int k)
{
	int max = Max(p, a, n, k);
	if (a[p[max]][k] == 0)
	{
		cout << "Вы ввели неадекватную матрицу, мда." << endl;
		return 0;
	}
	swap(p, p[k], max);

	double del = a[max][k];
	for (int j = k; j < n + 1; j++)
		a[max][j] /= del;

	for (int i = k + 1; i < n; i++)
	{
		del = a[p[i]][k];
		for (int j = k; j < n + 1; j++)
		{
			a[p[i]][j] -= (a[p[k]][j] * del);
		}
	}

	return 1;
}

void swap(int *p, int k1, int k2)
{
	int s1 = 0;
	while (p[s1] != k1)
		s1++;
	int s2 = 0;
	while (p[s2] != k2)
		s2++;
	p[s1] = k2;
	p[s2] = k1;
}

int Max(int *p, double **a, int n, int k)
{
	double max = 0;
	int maxs = p[k];

	for (int i = k; i < n; i++)
	{
		if (abs(a[p[i]][k]) > max)
		{
			max = abs(a[p[i]][k]);
			maxs = p[i];
		}
	}

	return maxs;
}


#define _CRT_SECURE_NO_WARNINGS
#pragma hdrstop
#pragma argsused
using namespace std;
#ifdef _WIN32
#include <tchar.h>
#else
typedef char _TCHAR;
#define _tmain main
#endif
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <string>
#include <filesystem>
#include <vector>
//размеры сетки
#define nd 7
#define M1 1000
#define M2 100
#define b 1
#ifdef __cplusplus
double max(double value1, double value2);
double max(double value1, double value2)
{
	return ((value1 > value2) ? value1 : value2);
}
#endif
#ifdef __cplusplus
double min(double value1, double value2);
double min(double value1, double value2)
{
	return ((value1 < value2) ? value1 : value2);
}
#endif
void Convert(char* buffer, long num)
{
	sprintf(buffer, "%i", num);
}
double  err[M1 + 1][M2 + 1], err1[M1 + 1][M2 + 1], denom[M1 + 1][M2 + 1], DJak[M1 + 1][M2 + 1], Z[M1 + 1][M2 + 1], fac[M1 + 1][M2 + 1],
Vr[M1 + 1][M2 + 1], Vz[M1 + 1][M2 + 1], RS[nd + 1][M1 + 1][M2 + 1], RS_tmp[nd + 1][M1 + 1][M2 + 1],
F1[nd + 1][M1 + 1][M2 + 1], F2[nd + 1][M1 + 1][M2 + 1], F3[nd + 1][M1 + 1][M2 + 1], F4[nd + 1][M1 + 1][M2 + 1], F5[nd + 1][M1 + 1][M2 + 1],
W[M1 + 1][M2 + 1], NU[M1 + 1][M2 + 1], Q[M1 + 1][M2 + 1], W_A[M1 + 1][M2 + 1], T_G[M1 + 1][M2 + 1], Jak[M1 + 1][M2 + 1],
Pe[M1 + 1][M2 + 1], ALFA[M1 + 1][M2 + 1], ALFAz[M1 + 1][2 * (M2 + 1)], Vrz[M1 + 1][2 * (M2 + 1)], Vzz[M1 + 1][2 * (M2 + 1)], Pz[M1 + 1][2 * (M2 + 1)], RS_prev[nd + 1][M1 + 1][M2 + 1];
int  L, R;
//координаты датчиков:
int dl1 = M1 * 0.3; int dr1 = b;
int dl2 = M1 * 0.3; int dr2 = M2 - 1;
int dl3 = M1 * 0.5; int dr3 = b;
int dl4 = M1 * 0.5; int dr4 = M2 - 1;
int dl5 = M1 * 0.7; int dr5 = b;
int dl6 = M1 * 0.7; int dr6 = M2 - 1;
double TZR1, TZR2, T, TZ1, TZ2, ZN1, ZN11, ZN111;
int  Rcl, Rad_cl, Zcl_1;
double C = 1500.0, B = 129.0, RR = 8.31,
LAMBDA = 257.e-4,
hZ = 1.e-3,
hR = 1.e-3,
TAU = 1.e-8,
CP = 1006.,
PI = 3.14,
T_0 = 300.0,
R_L0 = 1000.0,
R_G = 1.29,
CAPPA = LAMBDA / (CP * R_G),
GAMMA = 1.4;
void F(double RS_next[nd + 1][M1 + 1][M2 + 1], double RS_prev[nd + 1][M1 + 1][M2 + 1],double P_0,double DEL_P0,double A0)
{
	TZR1 = 1.0E-4;
	TZR2 = 5.0E-5 + TZR1;
	TZ1 = TZR1 / 2.0;
	ZN11 = 1.0 / (sqrt(4.0 * log(10.0)));
	for (int R = b; R <= M2; R++){
		if (T < TZR1 / 2.0) RS_prev[3][0][R] = DEL_P0 * exp(-((T - TZ1) / (TZ1 * ZN11)) * ((T - TZ1) / (TZ1 * ZN11))) + P_0;
		else RS_prev[3][0][R] = DEL_P0 + P_0;
	}
	for (int L = 1; L <= M1 - 1; L++)
		for (int R = b; R <= M2 - 1; R++){
			RS_next[1][L][R] = -(1 / (R_L0 * Jak[L][R])) * (RS_prev[3][L][R] - RS_prev[3][L - 1][R]) / hR;
			RS_next[2][L][R] = -(1 / (R_L0 * Jak[L][R])) * (RS_prev[3][L][R] - RS_prev[3][L][R - 1]) / hZ;
		}
	for (int L = 1; L <= M1 - 1; L++)
		for (int R = b; R <= M2 - 1; R++)
		{
			if (
				((L - Zcl_1) * (L - Zcl_1) + (R - Rcl) * (R - Rcl)) <= (Rad_cl * Rad_cl) // сфера
				){
				T_G[L][R] = (RS_prev[4][L][R] * RS_prev[4][L][R] * RS_prev[4][L][R] * RS_prev[6][L][R] * T_0) / (A0 * A0 * A0 * P_0);
				Pe[L][R] = 12. * (GAMMA - 1.0) * (T_0 * RS_prev[4][L][R] * fabs(W[L][R])) / (CAPPA * fabs(T_G[L][R] - (T_0 + 0.0000001)));
				if (Pe[L][R] > 100.) NU[L][R] = sqrt(Pe[L][R]);
				else NU[L][R] = 10.;
				Q[L][R] = NU[L][R] * LAMBDA * (T_G[L][R] - T_0) / (2. * RS_prev[4][L][R]);
				RS_next[5][L][R] = (1.0 / ((1.0 - RS_prev[5][L][R] / C) * RS_prev[4][L][R])) * ((RS_prev[4][L][R] / (R_L0 * C)) * (-(3. * GAMMA * RS_prev[6][L][R] * RS_prev[5][L][R] +
					(3. * (GAMMA - 1.0) * Q[L][R])) / RS_prev[4][L][R] - (C * C * R_L0 * (3.0 * RS_prev[7][L][R] * RS_prev[5][L][R] / RS_prev[4][L][R] - DJak[L][R] * (RS_prev[7][L][R] / Jak[L][R])) / (1.0 - RS_prev[7][L][R]))) +
					(1.0 + RS_prev[5][L][R] / C) * (RS_prev[6][L][R] - RS_prev[3][L][R]) / R_L0 - (1.0 - RS_prev[5][L][R] / (3.0 * C)) * 1.5 * RS_prev[5][L][R] * RS_prev[5][L][R]);
				W_A[L][R] = 0.0;
				W[L][R] = RS_prev[5][L][R] + W_A[L][R];
				RS_next[4][L][R] = RS_prev[5][L][R] + W_A[L][R];
				RS_next[6][L][R] = (-3.0 * GAMMA * RS_prev[6][L][R] * W[L][R] + 3. * (GAMMA - 1.0) * Q[L][R]) / RS_prev[4][L][R];
				DJak[L][R] = RS_prev[2][L][R] / (R * hR) + (RS_prev[1][L + 1][R] - RS_prev[1][L][R]) / hZ + (RS_prev[2][L][R + 1] - RS_prev[2][L][R]) / hR;
				RS_next[7][L][R] = 3.0 * RS_prev[7][L][R] * W[L][R] / RS_prev[4][L][R] - DJak[L][R] * RS_prev[7][L][R] / Jak[L][R];
				RS_next[3][L][R] = C * C * R_L0 * (3.0 * RS_prev[7][L][R] * W[L][R] / RS_prev[4][L][R] - DJak[L][R] + RS_prev[7][L][R]) / (1.0 - RS_prev[7][L][R]);
			}
			else {
				Jak[L][R] = 1.0;
				DJak[L][R] = RS_prev[2][L][R] / (R * hR) + (RS_prev[1][L + 1][R] - RS_prev[1][L][R]) / hZ + (RS_prev[2][L][R + 1] - RS_prev[2][L][R]) / hR;
				RS_next[3][L][R] = C * C * R_L0 * (-DJak[L][R]);
			}
			RS_next[3][M1 - 1][R] = C * R_L0 * (RS_prev[1][M1 - 1][R] - RS[1][M1 - 1][R]) / TAU; //условие протекания
			//RS_next[6][M1 - 1][R] = 0; // жесткая стенка
		}
}
void func(double alfa, double r,double P_0,double DEL_P0,double A0) {
	system("cls");
	cout << "\nAlfa = " << alfa << " Radius = " << r << " DelP0 = " << DEL_P0 << " A0 = " << A0 <<  endl;
	TAU = 1.e-8;
	double uround = 5.e-7,
		eps = 1.e-12,
		TAUmax = 1.e-6,
		TAUnew,
		ERRmin = 0.0,
		FACmax = 0.0;
	int L_max = L, R_max = R, K_max = 1;
	double T_max = 0, P_max = 1.e+5;
	int K, KK = 0, N;
	CAPPA = LAMBDA / (R_G * CP);
	int numer = 1;
	int numer1 = 1;
	int numer2 = 1;
	for (int L = 0; L <= M1; L++)
	{
		for (int R = 0; R <= M2; R++)
		{
			T_G[L][R] = 300;
			Jak[L][R] = 1.0;
			DJak[L][R] = 0.0;
			RS[1][L][R] = 0.0;
			RS[2][L][R] = 0.0;
			RS[3][L][R] = 1.e+5;
			RS[4][L][R] = A0;
			RS[5][L][R] = 0.0;
			RS[6][L][R] = 1.e+5;
			Jak[L][R] = 1.0;
			DJak[L][R] = 0.0;
		}

	}
	ofstream file2("alfa.dat", ios::out);
	for (int L = 0; L <= M1; L++)
	{
		for (int R = b; R <= M2; R++)
		{
			Rcl = 0;
			Rad_cl = r * 1000.;
			Zcl_1 = M1 / 2;//M1 * 2 * 0.025;
			if (((L - Zcl_1) * (L - Zcl_1) + (R - Rcl) * (R - Rcl)) <= (Rad_cl * Rad_cl)) RS[7][L][R] = alfa;
			else RS[7][L][R] = 1.e-8;
		}
	}
	for (int L = 0; L <= M1; L++) { RS[7][L][0] = RS[7][L][b]; }

	for (int L = 0; L <= M1; L++)
		for (int R = (b - 1); R <= 2 * (M2); R++)
		{
			if (R <= (M2)) ALFAz[L][R] = RS[7][L][(M2)-R];
			else if (R >= (M2)) ALFAz[L][R] = RS[7][L][R - (M2)];
		}
	for (int L = 0; L <= M1; L++)
		for (int R = (b - 1); R <= 2 * (M2); R++)
			file2 << "\n" << L * hZ << "\t" << -((M2)-R) * hR << "\t" << ALFAz[L][R];

	ofstream file3("Datchiki.dat", ios::out);
	ofstream file6("Time.dat", ios::out);
	T = 0.0;
	eps = max(1.e-12, 7.0 * uround);
	TAU = min(max(1.e-4, fabs(TAU)), TAUmax);
	eps = max(eps, 7.0 * uround);
	bool flag;
	TZR1 = 1.0E-4;
	if (T < TZR1)
		for (int R = 0; R <= M1; R++)
		{
			RS[3][0][R] = P_0;
			RS_tmp[3][0][R] = P_0;
		}
	F(F1, RS, P_0, DEL_P0, A0);
	while (T <= 0.001)
	{
		flag = true;
		if ((ceill(KK / 10) - numer1) == 0)
		{
			cout << "\n" << KK << "\t" << TAU << "\t";
			numer1 = numer1 + 1;
		}
		if (T < TZR1)
			for (int R = b; R <= M2; R++)
			{
				RS[3][0][R] = P_0;
				RS_tmp[3][0][R] = P_0;
			}
		while (flag)
		{
			flag = false;
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 1; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * 0.2 * F1[yn][L][R];
			//второй этап
			F(F2, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * ((3. / 40.) * F1[yn][L][R] + (9. / 40.) * F2[yn][L][R]);
			//третий этап
			F(F3, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * ((44. / 45.) * F1[yn][L][R] - (56. / 15.) * F2[yn][L][R] + (32. / 9.) * F3[yn][L][R]);
			//четвертый этап
			F(F4, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * ((19372. / 6561.) * F1[yn][L][R] - (25360. / 2187.) * F2[yn][L][R] + (64448. / 6561.) * F3[yn][L][R] - (212. / 729.) * F4[yn][L][R]);
			//пятый этап
			F(F5, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * ((9017. / 3168.) * F1[yn][L][R] - (355. / 33.) * F2[yn][L][R] + (46732. / 5247.) * F3[yn][L][R] + (49. / 176.) * F4[yn][L][R] - (5103. / 18656.) * F5[yn][L][R]);
			//шестой этап
			F(F2, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						RS_tmp[yn][L][R] = RS[yn][L][R] + TAU * ((35. / 384.) * F1[yn][L][R] + (500. / 1113.) * F3[yn][L][R] + (125. / 192.) * F4[yn][L][R] - (2187. / 6784.) * F5[yn][L][R] + (11. / 84.) * F2[yn][L][R]);
			//находим F2
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						F2[yn][L][R] = (71. / 57600.) * F1[yn][L][R] - (71. / 16695.) * F3[yn][L][R] + (71. / 1920.) * F4[yn][L][R] - (17253. / 339200.) * F5[yn][L][R] + (22. / 525.) * F2[yn][L][R];
			F(F3, RS_tmp, P_0, DEL_P0, A0);
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++)
						F4[yn][L][R] = (F2[yn][L][R] - (1. / 40.) * F3[yn][L][R]) * TAU;
			//upravlenie dlinoy shaga
			for (int L = 0; L <= M1; L++)
				for (int R = b; R <= M2; R++)
					err[L][R] = 0.;
			for (int yn = 1; yn <= nd; yn++)
				for (int L = 0; L <= M1; L++)
					for (int R = b; R <= M2; R++) {
						denom[L][R] = max(max(1.e-5, fabs(RS_tmp[yn][L][R])), max(fabs(RS[yn][L][R]), (2.0 * uround) / eps));
						err[L][R] = err[L][R] + (F4[yn][L][R] / denom[L][R]) * (F4[yn][L][R] / denom[L][R]);
					}

			for (int L = 0; L <= M1; L++)
				for (int R = b; R <= M2; R++)
				{
					err[L][R] = sqrt(err[L][R] / 5.);
					fac[L][R] = max(0.1, min(5., (pow(err[L][R] / eps, 0.2)) / 0.9));
				}
			FACmax = fac[0][0];
			ERRmin = err[0][0];
			for (int L = 0; L <= M1; L++)
				for (int R = b; R <= M2; R++)
				{
					if (fac[L][R] > FACmax) FACmax = fac[L][R];
					if (err[L][R] < ERRmin) ERRmin = err[L][R];
				}
			TAUnew = TAU / FACmax;
			if (ERRmin < eps)
			{
				for (int yn = 1; yn <= nd; yn++)
				{
					for (int L = 0; L <= M1; L++)
						for (int R = b; R <= M2; R++)
						{
							F1[yn][L][R] = F3[yn][L][R];
							RS[yn][L][R] = RS_tmp[yn][L][R];
							if (RS[3][L][R] > P_max)
							{
								P_max = RS[3][L][R];
								L_max = L;
								R_max = R;
								K_max = KK;
								T_max = T;
							}

						}
				}
				T = T + TAU;
				if (TAUnew > TAUmax) TAUnew = TAUmax;
				if (flag) TAUnew = min(TAUnew, TAU);
				flag = false;
			}
			else flag = true;
			if (TAUnew >= 1e-8) TAU = TAUnew;
			if (div(KK, 5).rem == 1)
			{ // ДАТЧИКИ
				//file3 << T << "\t" << (RS[3][dl1][dr1] / P_0) << "\t" << (RS[3][dl2][dr2] / P_0) << "\t" << (RS[3][dl3][dr3] / P_0) << "\t" << (RS[3][dl4][dr4] / P_0) << "\t" << (RS[3][dl5][dr5] / P_0) << "\t" << (RS[3][dl6][dr6] / P_0) << endl;
				//file3 << T << "\t" << RS[3][dl1][dr1]  << "\t" << RS[3][dl2][dr2] << "\t" << RS[3][dl3][dr3] << "\t" << RS[3][dl4][dr4] << "\t" << RS[3][dl5][dr5] << "\t" << RS[3][dl6][dr6] << endl;
				file3 << T << "\t" << RS[3][dl1][dr1] - P_0 << "\t" << RS[3][dl2][dr2] - P_0 << "\t" << RS[3][dl3][dr3] - P_0 << "\t" << RS[3][dl4][dr4] - P_0 << "\t" << RS[3][dl5][dr5] - P_0 << "\t" << RS[3][dl6][dr6] - P_0 << endl;
			}
			for (int L = 0; L <= M1 - 1; L++) RS[3][L][0] = RS[3][L][b];
			for (int L = 0; L <= M1 - 1; L++)
				for (int R = (b - 1); R <= 2 * (M2 - 1); R++)
				{
					if (R <= (M2 - 1)) Pz[L][R] = RS[3][L][(M2 - 1) - R];
					else if (R >= (M2 - 1)) Pz[L][R] = RS[3][L][R - (M2 - 1)];
				}
			for (int L = 0; L <= M1 - 1; L++) { RS[1][L][0] = RS[1][L][b]; }
			for (int L = 0; L <= M1 - 1; L++)
				for (int R = (b - 1); R <= 2 * (M2 - 1); R++)
				{
					if (R <= (M2 - 1)) Vzz[L][R] = RS[1][L][(M2 - 1) - R];
					else if (R >= (M2 - 1)) 
							Vzz[L][R] = RS[1][L][R - (M2 - 1)];
				}
			for (int L = 0; L <= M1 - 1; L++) { RS[2][L][0] = RS[2][L][b]; }
			for (int L = 0; L <= M1 - 1; L++)
				for (int R = (b - 1); R <= 2 * (M2 - 1); R++)
				{
					if (R <= (M2 - 1)) Vrz[L][R] = RS[2][L][(M2 - 1) - R];
					else if (R >= (M2 - 1)) Vrz[L][R] = RS[2][L][R - (M2 - 1)];
				}
			if (KK >= 1) {
				if (div(KK, 10).rem == 1)
				{
					char fname[40] = "FNAME";
					char fname2[40] = "FNAME1";
					Convert(fname, numer);
					Convert(fname2, numer2);
					ofstream file4(strcat(fname, ".dat"), ios::out);
					ofstream file4_2(strcat(fname, "_2d.dat"), ios::out);
					for (L = 0; L <= M1 - 1; L++) {
						for (R = (b - 1); R <= 2 * (M2 - 1); R++)
						{
							file4 << "\n" << L * hZ << "\t" << -((M2 - 1) - R) * hR << "\t" << (Pz[L][R] - P_0) / 1000000 << "\t" << T * 1000 << "\t" << (P_max - P_0) / P_0 << "\t" << P_max << "\t" << L_max << "\t" << R_max << "\t" << K_max << "\t" << T_max;
						}
					}
					for (L = 0; L <= M1 - 1; L++) {
						file4_2 << "\n" << L * hZ << "\t" << (Pz[L][b] - P_0) / 1000000 << "\t" << (Pz[L][M2-1] - P_0) / 1000000 << "\t" << T * 1000 << "\t" << (P_max - P_0) / P_0 << "\t" << P_max << "\t" << L_max << "\t" << R_max << "\t" << K_max << "\t" << T_max;
					}
					file6 << "\n" << numer << "\t" << TAU << "\t" << P_max << "\t" << T * 1000;
					numer = numer + 1;
					numer2 = numer2 + 1;
				}
			}
			KK = KK + 1;
		}
	}
}
vector<pair<double, double>> Pairs() {
	const double firstMin = 5e5;
	const double firstMax = 15e5;
	const double firstStep = 1e5;
	const double secondMin = 0.005;
	const double secondMax = 0.02;
	const double secondStep = 0.005;

	int firstSize = (firstMax - firstMin) / firstStep + 1;
	int secondSize = (secondMax - secondMin) / secondStep + 1;

	vector<pair<double, double>> pairs;
	for (double i = firstMin; i <= firstMax; i += firstStep) {
		for (double j = secondMin; j <= secondMax; j += secondStep) {
			pairs.push_back(std::make_pair(i, j));
		}
	}
	return pairs;
}
//=============================================================================
int _tmain(int argc, _TCHAR* argv[])
{
	/*
	_wchdir(L"Z:\\Программирование\\Diplom\\Diplom");
	ifstream Count("Count.txt");
	int Iter;
	Count >> Iter;
	vector<pair<double, double>> pairs = Pairs();
	double alfa = pairs[Iter].second;
	double delp = pairs[Iter].first;
	wstring dir = L"A" + to_wstring(alfa) + L"_R0.03_DP" + to_wstring(delp) + L"Mpa_A01e-3";
	_wmkdir(dir.c_str());
	_wchdir(dir.c_str());
	func(alfa, 0.03, 1e5, delp, 1e-3);*/
	/*
	double alfa = 0.0005;
	double rad = 0.01;
	wstring dir = L"A" + to_wstring(alfa) + L"_R" + to_wstring(rad);
	_wmkdir(dir.c_str());
	_wchdir(dir.c_str());
	func(alfa, rad, 1e5, 5e5, 1e-3);*/

	func(0.005, 0.03, 1e5, 5e5, 1e-3);
	return 0;
}
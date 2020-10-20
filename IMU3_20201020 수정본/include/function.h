//#pragma once
#ifndef FUNCTION_H
#define FUNCTION_H
#include "common.h"

// function related to linear algebra
void Velwisesquare(double v3in[], int size, double v3out[]); // Matlab .^2
void diag(double *vecin, int size, double *matout);
void transpose(double matin[], double matout[], int row, int col);
void MatMult(double matin[], double matin2[], double matout[], int l, int m, int n);
void MatMult31(double matin[][3], double vecin[], double vecout[]);
void MatMult41(double matin[][4], double vecin[], double vecout[]);
void MatSum(double matin[], double matin2[], int row, int col, double matout[]);
void MatSubtract(double matin[], double matin2[], int row, int col, double matout[]);
void MatSubstitute(double matin[], double matout[], int matincolsize, int matoutcolsize, int startrowout, int endrowout, int startcolout, int endcolout, int startrowin, int endrowin, int startcolin, int endcolin);
void MatScalarMult15(double matin[], double matout[], double scalar);
void MatScalarMult3(double matin[], double matout[], double scalar);
void MatScalarMult(double matin[], double matout[], double scalar, int row, int col);
void VecSubstitute(double *vin, double *vout, int size);
void VecScalarMult(double vin[], double scalarin, double vout[], int size);
void VecSum(double vin[], double vin2[], int size, double vout[]);
void VecSubtract(double vin[], double vin2[], int size, double vout[]);
void CrossProduct(double vin[], double vin2[], double vout[]);
void GetSkew(double vin[], double skewout[][3]);
void InvMat22(double A[][2], double invA[][2]);
void InvMat33(double A[][3], double invA[][3]);
double DetMat22(double A[][2]);
double DetMat33(double A[][3]);
void inv_mat(double *a, double *result, int np);

// function related to kinematics
void eulr2dcm(double r, double p, double y, double Cbnout[][3]);
void eulr2dcmtr(double r, double p, double y, double Cbnout[][3]);
void eulr2qua(double r, double p, double y, double quatout[4]);
void qua2dcm(double qbnin[], double Cbnout[][3]);
void dcm2eulr(double Cbnin[][3], double eulrout[]);
void QuatUpdate(double *wbnbin, double dtin, double *quatout);

// function related to filter
void GenINSFmatrix15(double *llh, double *vn, double Cbn[][3], double *fn, double *wnin, double *wnie, double *wnen, double F[][15]);
void InsKf15_predict(double F[][15], double *estX, double P_DB[][15], double *Qvec, double Cbn[][3], double Cbn_T[][3], double PHI[][15], double Pout[][15]);
void InsKf15_updateV(double P_DB[][15], double *estX, double R[][3], double H[][15], double *measurement, double Pout[][15], double *updatedX);
void InsKf15_updatePV(double P_DB[][15], double *estX, double R[][6], double H[][15], double *measurement, double Pout[][15], double *updatedX);
void InsKf15_updateHgt(double P_DB[][15], double *estX, double R[][1], double H[][15], double *measurement, double Pout[][15], double *updatedX);
void InsKf15_updatelatlon(double P_DB[][15], double *estX, double R[][2], double H[][15], double *measurement, double Pout[][15], double *updatedX);
void Correction(double *estX, double *llh, double *vn, double *quat, double Cbn[][3], double *ba, double *bg);

// function related to earth parameter
void Radicurv(double latin, double *rmout, double *rpout);
void delRadicurv(double latin, double *rmmout, double *rttout);
void NavRotation(double *vin, double *llhin, double *wninout, double *wnieout, double *wnenout);
double gravity_wgs84(double lat, double h);
double km2rad(double km);

#endif

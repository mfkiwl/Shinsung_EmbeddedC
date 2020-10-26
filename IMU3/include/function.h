//#pragma once
#ifndef FUNCTION_H
#define FUNCTION_H
#include "common.h"
#include "DSP28x_Project.h"

// function related to linear algebra
void Velwisesquare(float64 v3in[], int size, float64 v3out[]); // Matlab .^2
void diag(float64 *vecin, int size, float64 *matout);
void transpose(float64 matin[], float64 matout[], int row, int col);
void MatMult(float64 matin[], float64 matin2[], float64 matout[], int l, int m, int n);
void MatMult31(float64 matin[][3], float64 vecin[], float64 vecout[]);
void MatMult41(float64 matin[][4], float64 vecin[], float64 vecout[]);
void MatSum(float64 matin[], float64 matin2[], int row, int col, float64 matout[]);
void MatSubtract(float64 matin[], float64 matin2[], int row, int col, float64 matout[]);
void MatSubstitute(float64 matin[], float64 matout[], int matincolsize, int matoutcolsize, int startrowout, int endrowout, int startcolout, int endcolout, int startrowin, int endrowin, int startcolin, int endcolin);
void MatScalarMult15(float64 matin[], float64 matout[], float64 scalar);
void MatScalarMult3(float64 matin[], float64 matout[], float64 scalar);
void MatScalarMult(float64 matin[], float64 matout[], float64 scalar, int row, int col);
void VecSubstitute(float64 *vin, float64 *vout, int size);
void VecScalarMult(float64 vin[], float64 scalarin, float64 vout[], int size);
void VecSum(float64 vin[], float64 vin2[], int size, float64 vout[]);
void VecSubtract(float64 vin[], float64 vin2[], int size, float64 vout[]);
void CrossProduct(float64 vin[], float64 vin2[], float64 vout[]);
void GetSkew(float64 vin[], float64 skewout[][3]);
void InvMat22(float64 A[][2], float64 invA[][2]);
void InvMat33(float64 A[][3], float64 invA[][3]);
float64 DetMat22(float64 A[][2]);
float64 DetMat33(float64 A[][3]);
void inv_mat(float64 *a, float64 *result, int np);

// function related to kinematics
void eulr2dcm(float64 r, float64 p, float64 y, float64 Cbnout[][3]);
void eulr2dcmtr(float64 r, float64 p, float64 y, float64 Cbnout[][3]);
void eulr2qua(float64 r, float64 p, float64 y, float64 quatout[4]);
void qua2dcm(float64 qbnin[], float64 Cbnout[][3]);
void dcm2eulr(float64 Cbnin[][3], float64 eulrout[]);
void QuatUpdate(float64 *wbnbin, float64 dtin, float64 *quatout);

// function related to filter
void GenINSFmatrix15(float64 *llh, float64 *vn, float64 Cbn[][3], float64 *fn, float64 *wnin, float64 *wnie, float64 *wnen, float64 F[][15]);
void InsKf15_predict(float64 F[][15], float64 *estX, float64 P_DB[][15], float64 *Qvec, float64 Cbn[][3], float64 Cbn_T[][3], float64 PHI[][15], float64 Pout[][15]);
void InsKf15_updateV(float64 P_DB[][15], float64 *estX, float64 R[][3], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX);
void InsKf15_updatePV(float64 P_DB[][15], float64 *estX, float64 R[][6], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX);
void InsKf15_updateHgt(float64 P_DB[][15], float64 *estX, float64 R[][1], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX);
void InsKf15_updatelatlon(float64 P_DB[][15], float64 *estX, float64 R[][2], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX);
void InsKf15_updatezupt(float64 P_DB[][15], float64 *estX, float64 R[][5], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX);
void InsKf15_updatelatlonvnve(float64 P_DB[][15], float64 *estX, float64 R[][4], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX); ////////////////////**
void Correction(float64 *estX, float64 *llh, float64 *vn, float64 *quat, float64 Cbn[][3], float64 *ba, float64 *bg);

// function related to earth parameter
void Radicurv(float64 latin, float64 *rmout, float64 *rpout);
void delRadicurv(float64 latin, float64 *rmmout, float64 *rttout);
void NavRotation(float64 *vin, float64 *llhin, float64 *wninout, float64 *wnieout, float64 *wnenout);
float64 gravity_wgs84(float64 lat, float64 h);
float64 km2rad(float64 km);

#endif

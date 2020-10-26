#include "include/function.h"
#include "include/common.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

float64 km2rad(float64 km)
{
    return km / 6371.;
}

void Velwisesquare(float64 v3in[], int size, float64 v3out[])
{
    int i = 0;

    for (i = 0; i < size; i++)
    {
        v3out[i] = pow(v3in[i], 2);
    }
}

void diag(float64 *vecin, int size, float64 *matout)
{
    int i, j;

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (i == j)
                matout[i * size + j] = vecin[i];
            else
                matout[i * size + j] = 0.0;
        }
    }
}

void transpose(float64 matin[], float64 matout[], int row, int col)
{
    int i, j;

    for (i = 0; i < col; i++)
        for (j = 0; j < row; j++)
            matout[i * row + j] = matin[i + j * col];
}

// LxM * MxN
void MatMult(float64 matin[], float64 matin2[], float64 matout[], int L, int M, int N)
{
    int i, j, k;

    for (i = 0; i < L; i++)
    {
        for (k = 0; k < N; k++)
        {
            matout[i * N + k] = 0.0;

            for (j = 0; j < M; j++)
            {
                // cout << matin[i * M + j] * matin2[j * N + k] << endl;
                matout[i * N + k] += matin[i * M + j] * matin2[j * N + k];
            }
        }
    }
}

void MatMult31(float64 matin[][3], float64 vecin[], float64 vecout[])
{
    MatMult(matin[0], vecin, vecout, 3, 3, 1);
}

void MatMult41(float64 matin[][4], float64 vecin[], float64 vecout[])
{
    MatMult(matin[0], vecin, vecout, 4, 4, 1);
}

void MatSum(float64 matin[], float64 matin2[], int row, int col, float64 matout[])
{
    int i, j;

    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] + matin2[i * col + j];
}

void MatSubtract(float64 matin[], float64 matin2[], int row, int col, float64 matout[])
{
    int i, j;

    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] - matin2[i * col + j];
}

// matlab matout(startrowout:endrowout,startcolout:endcolout) = matin(startrowin:endrowin,startcolin:endcolin);
void MatSubstitute(float64 matin[], float64 matout[], int matincolsize, int matoutcolsize, int startrowout, int endrowout, int startcolout, int endcolout, int startrowin, int endrowin, int startcolin, int endcolin)
{
    int i, j;
    int rowsizeout = endrowout - startrowout;
    int colsizeout = endcolout - startcolout;
    for (i = 0; i < rowsizeout + 1; i++)
        for (j = 0; j < colsizeout + 1; j++)
            matout[(startrowout + i) * matoutcolsize + (startcolout + j)] = matin[(startrowin + i) * matincolsize + (startcolin + j)];
}

// 15x15 matrix * scalar
void MatScalarMult15(float64 matin[], float64 matout[], float64 scalar)
{
    int i, j;

    for (i = 0; i < 15; i++)
        for (j = 0; j < 15; j++)
            matout[i * 15 + j] = matin[i * 15 + j] * scalar;
}

// 3x3 matrix * scalar
void MatScalarMult3(float64 matin[], float64 matout[], float64 scalar)
{
    int i, j;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            matout[i * 3 + j] = matin[i * 3 + j] * scalar;
}

// general matrix * scalar
void MatScalarMult(float64 matin[], float64 matout[], float64 scalar, int row, int col)
{
    int i, j;

    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] * scalar;
}

void VecSubstitute(float64 *vin, float64 *vout, int size)
{
    int i;

    for (i = 0; i < size; i++)
        vout[i] = vin[i];
}

void VecScalarMult(float64 vin[], float64 scalarin, float64 vout[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        vout[i] = vin[i] * scalarin;
}

void VecSum(float64 *vin, float64 *vin2, int size, float64 *vout)
{
    int i;
    for (i = 0; i < size; i++)
        vout[i] = vin[i] + vin2[i];
}

void VecSubtract(float64 vin[], float64 vin2[], int size, float64 vout[])
{
    int i;
    for (i = 0; i < size; i++)
        vout[i] = vin[i] - vin2[i];
}

void CrossProduct(float64 vin[], float64 vin2[], float64 vout[])
{
    vout[0] = vin[1] * vin2[2] - vin[2] * vin2[1];
    vout[1] = vin[2] * vin2[0] - vin[0] * vin2[2];
    vout[2] = vin[0] * vin2[1] - vin[1] * vin2[0];
}

void GetSkew(float64 vin[], float64 skewout[][3])
{
    skewout[0][0] = 0.0;
    skewout[0][1] = -vin[2];
    skewout[0][2] = vin[1];
    skewout[1][0] = vin[2];
    skewout[1][1] = 0.0;
    skewout[1][2] = -vin[0];
    skewout[2][0] = -vin[1];
    skewout[2][1] = vin[0];
    skewout[2][2] = 0.0;
}

//Inverse of 2x2 matrix
void InvMat22(float64 A[][2], float64 invA[][2])
{
    float64 detA = DetMat22(A);

    invA[0][0] = A[1][1] / detA;
    invA[0][1] = -A[0][1] / detA;
    invA[1][0] = -A[1][0] / detA;
    invA[1][1] = A[0][0] / detA;
}

//Inverse of 3x3 matrix
void InvMat33(float64 A[][3], float64 invA[][3])
{
    float64 detA = DetMat33(A);

    invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / detA;
    invA[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / detA;
    invA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / detA;
    invA[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / detA;
    invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / detA;
    invA[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / detA;
    invA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detA;
    invA[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / detA;
    invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / detA;
}

//Determinant of 2x2 matrix
float64 DetMat22(float64 A[][2])
{
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
}

//Determinant of 3x3 matrix
float64 DetMat33(float64 A[][3])
{
    return (A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1] -
            A[0][1] * A[1][0] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
            A[0][2] * A[1][0] * A[2][1] - A[0][2] * A[1][1] * A[2][0]);
}

#define DBL_EPSILON 2.2204460492503131e-018 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
// 2.2204460492503131e-016�� ������
#define POS_EPS DBL_EPSILON
#define mat(name, col_size, r, c) (*(name + r * col_size + c))

void inv_mat(float64 *a, float64 *result, int np)
{
    int i, j, ii, maxi;
    float64 p[15][15];
    float64 temp;
    //float64 EPS = POS_EPS;

    for (i = 0; i < np; i++)
    {
        for (j = 0; j < np; j++)
        {
            if (i == j)
                mat(result, np, i, j) = 1.0;
            else
                mat(result, np, i, j) = 0.0;
            p[i][j] = mat(a, np, i, j);
        }
    }

    for (ii = 0; ii < np; ii++)
    {
        maxi = ii;
        for (i = ii; i < np; i++)
        {
            if (fabs(p[maxi][ii]) < fabs(p[i][ii]))
                maxi = i;
        }

        if (fabs(p[maxi][ii]) < POS_EPS)
        {
            printf("\n MATRIX INVERSION ERROR !!! \n");
        }

        for (j = 0; j < np; j++)
        {
            temp = p[ii][j];
            p[ii][j] = p[maxi][j];
            p[maxi][j] = temp;

            temp = mat(result, np, ii, j);
            mat(result, np, ii, j) = mat(result, np, maxi, j);
            mat(result, np, maxi, j) = temp;
        }

        temp = p[ii][ii];
        for (i = 0; i < np; i++)
        {
            mat(result, np, ii, i) = mat(result, np, ii, i) / temp;
            p[ii][i] = p[ii][i] / temp;
        }

        for (i = 0; i < np; i++)
        {
            if (i != ii)
            {
                temp = p[i][ii];

                for (j = 0; j < np; j++)
                {
                    mat(result, np, i, j) = mat(result, np, i, j) - temp * mat(result, np, ii, j);
                    p[i][j] = p[i][j] - temp * p[ii][j];
                }
            }
        }
    }
}

//Cbn matrix
void eulr2dcm(float64 r, float64 p, float64 y, float64 Cbnout[][3])
{
    float64 Cpsi = cos(y);
    float64 Spsi = sin(y);
    float64 Cthe = cos(p);
    float64 Sthe = sin(p);
    float64 Cphi = cos(r);
    float64 Sphi = sin(r);

    Cbnout[0][0] = Cthe * Cpsi;
    Cbnout[1][0] = Cthe * Spsi;
    Cbnout[2][0] = -Sthe;
    Cbnout[0][1] = -Cphi * Spsi + Sphi * Sthe * Cpsi;
    Cbnout[1][1] = Cphi * Cpsi + Sphi * Sthe * Spsi;
    Cbnout[2][1] = Sphi * Cthe;
    Cbnout[0][2] = Sphi * Spsi + Cphi * Sthe * Cpsi;
    Cbnout[1][2] = -Sphi * Cpsi + Cphi * Sthe * Spsi;
    Cbnout[2][2] = Cphi * Cthe;
}

void eulr2dcmtr(float64 r, float64 p, float64 y, float64 Cnbout[][3])
{
    float64 Cpsi = cos(y);
    float64 Spsi = sin(y);
    float64 Cthe = cos(p);
    float64 Sthe = sin(p);
    float64 Cphi = cos(r);
    float64 Sphi = sin(r);

    Cnbout[0][0] = Cthe * Cpsi;
    Cnbout[0][1] = Cthe * Spsi;
    Cnbout[0][2] = -Sthe;
    Cnbout[1][0] = -Cphi * Spsi + Sphi * Sthe * Cpsi;
    Cnbout[1][1] = Cphi * Cpsi + Sphi * Sthe * Spsi;
    Cnbout[1][2] = Sphi * Cthe;
    Cnbout[2][0] = Sphi * Spsi + Cphi * Sthe * Cpsi;
    Cnbout[2][1] = -Sphi * Cpsi + Cphi * Sthe * Spsi;
    Cnbout[2][2] = Cphi * Cthe;
}

//qbn quaternion
void eulr2qua(float64 r, float64 p, float64 y, float64 quatout[4])
{
    float64 cpsi2 = cos(y / 2);
    float64 spsi2 = sin(y / 2);
    float64 cthe2 = cos(p / 2);
    float64 sthe2 = sin(p / 2);
    float64 cphi2 = cos(r / 2);
    float64 sphi2 = sin(r / 2);
    quatout[0] = cphi2 * cthe2 * cpsi2 + sphi2 * sthe2 * spsi2;
    quatout[1] = sphi2 * cthe2 * cpsi2 - cphi2 * sthe2 * spsi2;
    quatout[2] = cphi2 * sthe2 * cpsi2 + sphi2 * cthe2 * spsi2;
    quatout[3] = cphi2 * cthe2 * spsi2 - sphi2 * sthe2 * cpsi2;
}

//qbn to Cbn
void qua2dcm(float64 qbnin[], float64 Cbnout[][3])
{
    float64 a = qbnin[0];
    float64 b = qbnin[1];
    float64 c = qbnin[2];
    float64 d = qbnin[3];

    Cbnout[0][0] = a * a + b * b - c * c - d * d;
    Cbnout[0][1] = 2 * (b * c - a * d);
    Cbnout[0][2] = 2 * (b * d + a * c);
    Cbnout[1][0] = 2 * (b * c + a * d);
    Cbnout[1][1] = a * a - b * b + c * c - d * d;
    Cbnout[1][2] = 2.0 * (c * d - a * b);
    Cbnout[2][0] = 2.0 * (b * d - a * c);
    Cbnout[2][1] = 2.0 * (c * d + a * b);
    Cbnout[2][2] = a * a - b * b - c * c + d * d;
}

void dcm2eulr(float64 Cbnin[][3], float64 eulrout[])
{
    float64 roll = atan(Cbnin[2][1] / Cbnin[2][2]);
    float64 pitch = asin(-Cbnin[2][0]);
    float64 yaw = atan2(Cbnin[1][0], Cbnin[0][0]);
    float64 rpy[3] = {roll, pitch, yaw};
    VecSubstitute(rpy, eulrout, 3);
}

void QuatUpdate(float64 *wbnbin, float64 dtin, float64 *quatout)
{
    float64 dtheta[3];
    VecScalarMult(wbnbin, dtin, dtheta, 3);
    float64 Wx = dtheta[0];
    float64 Wy = dtheta[1];
    float64 Wz = dtheta[2];
    float64 theta = sqrt(Wx * Wx + Wy * Wy + Wz * Wz);
    float64 ac = cos(theta / 2);
    float64 as;
    if (theta < 1e-15)
        as = 0;
    else
        as = sin(theta / 2) / theta;
    float64 rk[4] = {ac, as * Wx, as * Wy, as * Wz};
    float64 Rk[4][4] = {{rk[0], -rk[1], -rk[2], -rk[3]}, {rk[1], rk[0], rk[3], -rk[2]}, {rk[2], -rk[3], rk[0], rk[1]}, {rk[3], rk[2], -rk[1], rk[0]}};
    float64 qtmp[4];
    MatMult41(Rk, quatout, qtmp);
    float64 mag = sqrt(pow(qtmp[0], 2) + pow(qtmp[1], 2) + pow(qtmp[2], 2) + pow(qtmp[3], 2));
    float64 q_normal = 1.5 - 0.5 * mag * mag;
    VecScalarMult(qtmp, q_normal, qtmp, 4);
    VecSubstitute(qtmp, quatout, 4);
}

void GenINSFmatrix15(float64 *llh, float64 *vn, float64 Cbn[][3], float64 *fn, float64 *wnin, float64 *wnie, float64 *wnen, float64 F[][15])
{
    memset(F, 0x0, sizeof(float64) * 15 * 15);

    float64 Rm, Rt, Rmm, Rtt;
    Radicurv(llh[0], &Rm, &Rt);
    delRadicurv(llh[0], &Rmm, &Rtt);
    float64 CL = cos(llh[0]);
    float64 TL = tan(llh[0]);
    float64 h = llh[2];
    float64 vn_n = vn[0];
    float64 vn_e = vn[1];
    float64 vn_d = vn[2];
    float64 wnie_n = wnie[0];
    float64 wnie_d = wnie[2];
    float64 wnen_n = wnen[0];
    float64 wnen_e = wnen[1];
    float64 wnen_d = wnen[2];

    float64 F11[3][3] = {{(wnen_e * Rmm) / (Rm + h), 0, wnen_e / (Rm + h)},
                        {wnen_n * (TL - Rtt / (Rt + h)) / CL, 0, -wnen_n / (Rt + h) / CL},
                        {0, 0, 0}};
    MatSubstitute((float64 *)F11, (float64 *)F, 3, 15, 0, 2, 0, 2, 0, 2, 0, 2);

    float64 F12[3][3] = {{1 / (Rm + h), 0, 0},
                        {0, 1 / (Rt + h) / CL, 0},
                        {0, 0, -1}};
    MatSubstitute((float64 *)F12, (float64 *)F, 3, 15, 0, 2, 3, 5, 0, 2, 0, 2);

    float64 F21[3][3] = {{(-vn_e * (2 * wnie_n + wnen_n / CL / CL) + (vn_d * wnen_e * Rmm) / (Rm + h) - wnen_n * wnen_d * Rtt), 0, ((wnen_e * vn_d) / (Rm + h) - wnen_n * wnen_d)},
                        {-vn_d * (wnen_n * Rtt / (Rt + h) - 2 * wnie_d) + vn_n * (2 * wnie_n + wnen_n / CL / CL + wnen_d * Rtt / (Rt + h)), 0, vn_n * wnen_d / (Rt + h) - wnen_n * vn_d / (Rt + h)},
                        {-2 * wnie_d * vn_e + (pow(wnen_n, 2.)) * Rtt + (pow(wnen_e, 2.)) * Rmm, 0, pow(wnen_n, 2.) + pow(wnen_e, 2.)}};

    MatSubstitute((float64 *)F21, (float64 *)F, 3, 15, 3, 5, 0, 2, 0, 2, 0, 2);

    float64 F221[3][3] = {{vn_d / (Rm + h), wnen_d, 0},
                         {0, (vn_d + vn_n * TL) / (Rt + h), 0},
                         {wnen_e, -wnen_n, 0}};

    float64 F222[3][3];
    memset(F222, 0x0, sizeof(float64) * 3 * 3);

    float64 wnsubtract[3];
    memset(wnsubtract, 0x0, sizeof(float64) * 3);

    float64 wnietmp[3];
    memset(wnietmp, 0x0, sizeof(float64) * 3);

    VecScalarMult(wnie, -2, wnietmp, 3);
    MatSubtract(wnietmp, wnen, 3, 1, wnsubtract);
    GetSkew(wnsubtract, F222);

    float64 F22[3][3];
    memset(F22, 0x0, sizeof(float64) * 3 * 3);

    MatSum((float64 *)F221, (float64 *)F222, 3, 3, (float64 *)F22);
    MatSubstitute((float64 *)F22, (float64 *)F, 3, 15, 3, 5, 3, 5, 0, 2, 0, 2);

    float64 F23[3][3];
    memset(F23, 0x0, sizeof(float64) * 3 * 3);

    GetSkew(fn, F23);
    MatSubstitute((float64 *)F23, (float64 *)F, 3, 15, 3, 5, 6, 8, 0, 2, 0, 2);

    float64 F24[3][3];
    memset(F24, 0x0, sizeof(float64) * 3 * 3);

    MatSubstitute((float64 *)Cbn, (float64 *)F24, 3, 3, 0, 2, 0, 2, 0, 2, 0, 2);
    MatSubstitute((float64 *)F24, (float64 *)F, 3, 15, 3, 5, 9, 11, 0, 2, 0, 2);

    float64 F31[3][3] = {{wnie_d - wnen_n * Rtt / (Rt + h), 0, -wnen_n / (Rt + h)},
                        {-wnen_e * Rmm / (Rm + h), 0, -wnen_e / (Rm + h)},
                        {-wnie_n - wnen_n / CL / CL - wnen_d * Rtt / (Rt + h), 0, -wnen_d / (Rt + h)}};
    MatSubstitute((float64 *)F31, (float64 *)F, 3, 15, 6, 8, 0, 2, 0, 2, 0, 2);

    float64 F32[3][3] = {{0, 1 / (Rt + h), 0},
                        {-1 / (Rm + h), 0, 0},
                        {0, -TL / (Rt + h), 0}};
    MatSubstitute((float64 *)F32, (float64 *)F, 3, 15, 6, 8, 3, 5, 0, 2, 0, 2);

    float64 F33[3][3];
    memset(F33, 0x0, sizeof(float64) * 3 * 3);

    float64 minuswnin[3];
    memset(minuswnin, 0x0, sizeof(float64) * 3);

    VecScalarMult(wnin, -1, minuswnin, 3);
    GetSkew(minuswnin, F33);
    MatSubstitute((float64 *)F33, (float64 *)F, 3, 15, 6, 8, 6, 8, 0, 2, 0, 2);

    float64 F35[3][3];
    memset(F35, 0x0, sizeof(float64) * 3 * 3);

    MatScalarMult3((float64 *)Cbn, (float64 *)F35, -1);
    MatSubstitute((float64 *)F35, (float64 *)F, 3, 15, 6, 8, 12, 14, 0, 2, 0, 2);
}

void InsKf15_predict(float64 F[][15], float64 *estX, float64 P_DB[][15], float64 *Qvec, float64 Cbn[][3], float64 Cbn_T[][3], float64 PHI[][15], float64 Pout[][15])
{
    // for (int i = 0; i < 15; i++)
    //     for (int j = 0; j < 15; j++)
    //         PHI[i][j] = 0; // ���ϸ� �����Ⱚ�� ���� ��� �߻�

    //	memset(

    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15], Qm[15][15] = {0};
    float64 A[3][3], A_tmp[3][3], G[3][3], G_tmp[3][3], Qmat1[3][3], Qmat2[3][3] = {0};
    float64 Qvec1[3] = {Qvec[0], Qvec[1], Qvec[2]};
    float64 Qvec2[3] = {Qvec[3], Qvec[4], Qvec[5]};
    float64 PHI_tmp1[15][15], PHI_tmp2[15][15], PHI_tmp3[15][15], PHI_tmp4[15][15] = {0};
    float64 P_tmp1[15][15], P_tmp2[15][15], P_tmp3[15][15], P_tmp4[15][15] = {0};

    diag(eyevec, 15, (float64 *)eyemat);
    diag(Qvec1, 3, (float64 *)Qmat1);
    diag(Qvec2, 3, (float64 *)Qmat2);
    MatMult((float64 *)Cbn, (float64 *)Qmat1, (float64 *)A_tmp, 3, 3, 3);
    MatMult((float64 *)A_tmp, (float64 *)Cbn_T, (float64 *)A, 3, 3, 3);
    MatMult((float64 *)Cbn, (float64 *)Qmat2, (float64 *)G_tmp, 3, 3, 3);
    MatMult((float64 *)G_tmp, (float64 *)Cbn_T, (float64 *)G, 3, 3, 3);
    MatSubstitute((float64 *)A, (float64 *)Qm, 3, 15, 3, 5, 3, 5, 0, 2, 0, 2);
    MatSubstitute((float64 *)G, (float64 *)Qm, 3, 15, 6, 8, 6, 8, 0, 2, 0, 2);

    // discretize F to PHI
    float64 F_tmp[15][15] = {0};
    MatSubstitute((float64 *)F, (float64 *)F_tmp, 15, 15, 0, 14, 0, 14, 0, 14, 0, 14);
    MatScalarMult15((float64 *)F, (float64 *)PHI_tmp1, dt);
    MatMult((float64 *)F, (float64 *)F_tmp, (float64 *)PHI_tmp2, 15, 15, 15);
    MatScalarMult15((float64 *)PHI_tmp2, (float64 *)PHI_tmp3, pow(dt, 2) / 2);
    MatSum((float64 *)eyemat, (float64 *)PHI_tmp1, 15, 15, (float64 *)PHI_tmp4);
    MatSum((float64 *)PHI_tmp4, (float64 *)PHI_tmp3, 15, 15, (float64 *)PHI);

    // Time update of error covariance
    MatMult((float64 *)PHI, (float64 *)P_DB, (float64 *)P_tmp1, 15, 15, 15);
    transpose((float64 *)PHI, (float64 *)P_tmp2, 15, 15);
    MatMult((float64 *)P_tmp1, (float64 *)P_tmp2, (float64 *)P_tmp3, 15, 15, 15);
    MatScalarMult15((float64 *)Qm, (float64 *)P_tmp4, dt);
    MatSum((float64 *)P_tmp3, (float64 *)P_tmp4, 15, 15, (float64 *)Pout);
}

// 10.20.10:30
void InsKf15_updateV(float64 P_DB[][15], float64 *estX, float64 R[][3], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX)
{
    static float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static float64 eyemat[15][15] = {0};
    static float64 K[15][3], Ktmp1[15][3], K_T[3][15] = {0};
    static float64 HPb[3][15], H_T[15][3] = {0};
    static float64 invtmp1[3][3], invtmp2[3][3], inverse[3][3] = {0};
    static float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    static float64 KRK[15][15], KRK_tmp1[15][3] = {0};
    static float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    static float64 Hx[3], residual[3], updatedX_tmp[15] = {0};
    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 3, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 3, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 3, 15, 3);
    MatSum((float64 *)invtmp1, (float64 *)R, 3, 3, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 3);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 3);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 3, 3);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 3, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 3);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 3, 3);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 3, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 3, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 3, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 3, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void InsKf15_updatePV(float64 P_DB[][15], float64 *estX, float64 R[][6], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX)
{
    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15] = {0};
    float64 K[15][6], Ktmp1[15][6], K_T[6][15] = {0};
    float64 HPb[6][15], H_T[15][6] = {0};
    float64 invtmp1[6][6], invtmp2[6][6], inverse[6][6] = {0};
    float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    float64 KRK[15][15], KRK_tmp1[15][6] = {0};
    float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    float64 Hx[6], residual[6], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 6, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 6, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 6, 15, 6);
    MatSum((float64 *)invtmp1, (float64 *)R, 6, 6, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 6);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 6);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 6, 6);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 6, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 6);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 6, 6);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 6, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 6, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 6, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 6, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void InsKf15_updateHgt(float64 P_DB[][15], float64 *estX, float64 R[][1], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX)
{
    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15] = {0};
    float64 K[15][1], Ktmp1[15][1], K_T[1][15] = {0};
    float64 HPb[1][15], H_T[15][1] = {0};
    float64 invtmp1[1][1], invtmp2[1][1], inverse[1][1] = {0};
    float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    float64 KRK[15][15], KRK_tmp1[15][1] = {0};
    float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    float64 Hx[1], residual[1], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 1, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 1, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 1, 15, 1);
    MatSum((float64 *)invtmp1, (float64 *)R, 1, 1, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 1);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 1);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 1, 1);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 1, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 1);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 1, 1);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 1, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 1, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 1, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 1, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void InsKf15_updatelatlon(float64 P_DB[][15], float64 *estX, float64 R[][2], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX)
{
    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15] = {0};
    float64 K[15][2], Ktmp1[15][2], K_T[2][15] = {0};
    float64 HPb[2][15], H_T[15][2] = {0};
    float64 invtmp1[2][2], invtmp2[2][2], inverse[2][2] = {0};
    float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    float64 KRK[15][15], KRK_tmp1[15][2] = {0};
    float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    float64 Hx[2], residual[2], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 2, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 2, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 2, 15, 2);
    MatSum((float64 *)invtmp1, (float64 *)R, 2, 2, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 2);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 2);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 2, 2);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 2, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 2);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 2, 2);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 2, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 2, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 2, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 2, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void InsKf15_updatelatlonvnve(float64 P_DB[][15], float64 *estX, float64 R[][4], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX) ////////////////////**
{
    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15] = {0};
    float64 K[15][4], Ktmp1[15][4], K_T[4][15] = {0};
    float64 HPb[4][15], H_T[15][4] = {0};
    float64 invtmp1[4][4], invtmp2[4][4], inverse[4][4] = {0};
    float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    float64 KRK[15][15], KRK_tmp1[15][4] = {0};
    float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    float64 Hx[4], residual[4], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 4, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 4, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 4, 15, 4);
    MatSum((float64 *)invtmp1, (float64 *)R, 4, 4, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 4);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 4);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 4, 4);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 4, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 4);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 4, 4);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 4, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 4, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 4, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 4, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void InsKf15_updatezupt(float64 P_DB[][15], float64 *estX, float64 R[][5], float64 H[][15], float64 *measurement, float64 Pout[][15], float64 *updatedX) 
{
    float64 eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float64 eyemat[15][15] = {0};
    float64 K[15][5], Ktmp1[15][5], K_T[5][15] = {0};
    float64 HPb[5][15], H_T[15][5] = {0};
    float64 invtmp1[5][5], invtmp2[5][5], inverse[5][5] = {0};
    float64 eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    float64 KRK[15][15], KRK_tmp1[15][5] = {0};
    float64 updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    float64 Hx[5], residual[5], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((float64 *)eyevec, 15, (float64 *)eyemat);
    MatMult((float64 *)H, (float64 *)P_DB, (float64 *)HPb, 5, 15, 15);
    transpose((float64 *)H, (float64 *)H_T, 5, 15);
    MatMult((float64 *)HPb, (float64 *)H_T, (float64 *)invtmp1, 5, 15, 5);
    MatSum((float64 *)invtmp1, (float64 *)R, 5, 5, (float64 *)invtmp2);
    inv_mat((float64 *)invtmp2, (float64 *)inverse, 5);
    MatMult((float64 *)P_DB, (float64 *)H_T, (float64 *)Ktmp1, 15, 15, 5);
    MatMult((float64 *)Ktmp1, (float64 *)inverse, (float64 *)K, 15, 5, 5);

    // calculating updated P
    MatMult((float64 *)K, (float64 *)H, (float64 *)eyeminusKH_tmp, 15, 5, 15);
    MatSubtract((float64 *)eyemat, (float64 *)eyeminusKH_tmp, 15, 15, (float64 *)eyeminusKH);
    transpose((float64 *)eyeminusKH, (float64 *)eyeminusKH_T, 15, 15);
    transpose((float64 *)K, (float64 *)K_T, 15, 5);
    MatMult((float64 *)K, (float64 *)R, (float64 *)KRK_tmp1, 15, 5, 5);
    MatMult((float64 *)KRK_tmp1, (float64 *)K_T, (float64 *)KRK, 15, 5, 15);
    MatMult((float64 *)eyeminusKH, (float64 *)P_DB, (float64 *)updatedP_tmp, 15, 15, 15);
    MatMult((float64 *)updatedP_tmp, (float64 *)eyeminusKH_T, (float64 *)updatedP_tmp2, 15, 15, 15);
    MatSum((float64 *)updatedP_tmp2, (float64 *)KRK, 15, 15, (float64 *)Pout);

    // calculating updated state
    MatMult((float64 *)H, (float64 *)estX, (float64 *)Hx, 5, 15, 1);
    MatSubtract((float64 *)measurement, (float64 *)Hx, 5, 1, (float64 *)residual);
    MatMult((float64 *)K, (float64 *)residual, (float64 *)updatedX_tmp, 15, 5, 1);
    VecSum((float64 *)estX, (float64 *)updatedX_tmp, 15, (float64 *)updatedX);
}

void Correction(float64 *estX, float64 *llh, float64 *vn, float64 *quat, float64 Cbn[][3], float64 *ba, float64 *bg)
{
    // bias ���Ǹ� �ݴ�� ����� ������(estimated F,Wbib = true - bias -> true = estimated + bias) correction�� ���� �ٿ���
    estX[9] = -estX[9];
    estX[10] = -estX[10];
    estX[11] = -estX[11];
    estX[12] = -estX[12];
    estX[13] = -estX[13];
    estX[14] = -estX[14];
    float64 estPos[3] = {estX[0], estX[1], estX[2]};
    float64 estVel[3] = {estX[3], estX[4], estX[5]};
    float64 estAtt[3] = {estX[6], estX[7], estX[8]};
    float64 estBa[3] = {estX[9], estX[10], estX[11]};
    float64 estBg[3] = {estX[12], estX[13], estX[14]};

    VecSubtract(llh, estPos, 3, llh);
    VecSubtract(vn, estVel, 3, vn);
    VecSubtract(ba, estBa, 3, ba);
    VecSubtract(bg, estBg, 3, bg);

    float64 Q_tmp[4][3] = {{-quat[1], -quat[2], -quat[3]},
                          {quat[0], quat[3], -quat[2]},
                          {-quat[3], quat[0], quat[1]},
                          {quat[2], -quat[1], quat[0]}};

    float64 del_Qtmp[4][3], del_Q[4] = {0};

    MatScalarMult((float64 *)Q_tmp, (float64 *)del_Qtmp, -0.5, 4, 3);
    MatMult((float64 *)del_Qtmp, (float64 *)estAtt, (float64 *)del_Q, 4, 3, 1);
    VecSubtract((float64 *)quat, (float64 *)del_Q, 4, (float64 *)quat);
    float64 qn = 1.5 - 0.5 * (pow(quat[0], 2) + pow(quat[1], 2) + pow(quat[2], 2) + pow(quat[3], 2));
    VecScalarMult(quat, qn, quat, 4);
    qua2dcm(quat, Cbn);
    VecScalarMult(estX, 0, estX, 15); // correction ���� ����ġ�� ������� ����ȭ �ϱ� ������ ������ ���� ������ ����(EKF)
}

void Radicurv(float64 latin, float64 *rmout, float64 *rpout)
{
    float64 e2 = pow(e, 2);
    float64 den = 1 - e2 * pow(sin(latin), 2);
    *rmout = R0 * (1 - e2) / pow(den, 1.5);
    *rpout = R0 / sqrt(den);
}

void delRadicurv(float64 latin, float64 *rmmout, float64 *rttout)
{
    float64 e2 = pow(e, 2);
    float64 den = 1 - e2 * pow(sin(latin), 2);
    *rmmout = (3 * R0 * (1 - e2) * e2 * sin(latin) * cos(latin)) / (pow(den, 2.5));
    *rttout = (R0 * e2 * sin(latin) * cos(latin)) / (pow(den, 1.5));
}

void NavRotation(float64 *vin, float64 *llhin, float64 *wninout, float64 *wnieout, float64 *wnenout)
{
    int i;
    float64 rm;
    float64 rp;
    Radicurv(llhin[0], &rm, &rp);
    float64 dlat = vin[0] / (rm + llhin[2]);
    float64 dlon = vin[1] / (rp + llhin[2]) / cos(llhin[0]);
    float64 clat = cos(llhin[0]);
    float64 slat = sin(llhin[0]);

    wnieout[0] = We * clat;
    wnieout[1] = 0;
    wnieout[2] = -We * slat;
    wnenout[0] = dlon * clat;
    wnenout[1] = -dlat;
    wnenout[2] = -dlon * slat;
    for (i = 0; i < 3; i++)
        wninout[i] = wnieout[i] + wnenout[i];
}

float64 gravity_wgs84(float64 lat, float64 h)
{
    float64 a = 6378137.0;          // Semi - major axis[m] : equatorial radius of the earth; WGS - 84
    float64 ee = 0.081819190842622; // eccentricity of the earth ellipsoid
    float64 ge = 9.7803253359;      // Theoretical gravity at the Equator[m / s2]
    float64 k = 0.00193185265241;   // Theoretical gravity formula constant
    float64 m = 0.00344978650684;   // m = ER ^ 2 * a ^ 2 * b / GM;

    float64 e2 = pow(ee, 2.0);
    float64 a2 = pow(a, 2.0);
    float64 slat = sin(lat);
    float64 slat2 = pow(slat, 2.0);
    float64 slat4 = pow(slat2, 2.0);
    float64 h2 = pow(h, 2.0);

    float64 g1 = k + e2 / 2;
    float64 g2 = -2.0 * (1.0 + f + m) / a;
    float64 g3 = g1 * g2 + 4.0 * f / a;
    float64 g4 = 3.0 / a2;
    float64 g5 = g1 * g4;

    return ge * (1.0 + g1 * slat2 + g2 * h + g3 * slat2 * h + g4 * h2 + g5 * slat2 * h2);
}

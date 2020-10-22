#include "function.h"
#include "common.h"
#include <iostream>
#include <cmath>

// Functions
// vector<double> &Arr2vec(double arrin[], int length, vector<double> &vecout)
// {
//     for (int i = 0; i < length; i++)
//     {
//         vecout.push_back(arrin[i]);
//     }
//     return vecout;
// }

double km2rad(double km)
{
    return km / 6371.;
}

void Velwisesquare(double v3in[], int size, double v3out[])
{
    for (int i = 0; i < size; i++)
    {
        v3out[i] = pow(v3in[i], 2);
    }
}

void diag(double *vecin, int size, double *matout)
{
    int i = 0;
    int j = 0;
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

double mean(double *&vecin, int l)
{
    double tmp = 0.0;
    int idx = 0;
    for (idx = 0; idx < l; idx++)
        tmp = tmp + vecin[idx];
    return tmp / (double)l;
}

void transpose(double matin[], double matout[], int row, int col)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < col; i++)
        for (j = 0; j < row; j++)
            matout[i * row + j] = matin[i + j * col];
}

// LxM * MxN
void MatMult(double matin[], double matin2[], double matout[], int L, int M, int N)
{
    int i = 0;
    int j = 0;
    int k = 0;

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

void MatMult31(double matin[][3], double vecin[], double vecout[])
{
    MatMult(matin[0], vecin, vecout, 3, 3, 1);
}

void MatMult41(double matin[][4], double vecin[], double vecout[])
{
    MatMult(matin[0], vecin, vecout, 4, 4, 1);
}

void MatSum(double matin[], double matin2[], int row, int col, double matout[])
{
    int i, j;

    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] + matin2[i * col + j];
}

void MatSubtract(double matin[], double matin2[], int row, int col, double matout[])
{
    int i, j;

    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] - matin2[i * col + j];
}

// matlab matout(startrowout:endrowout,startcolout:endcolout) = matin(startrowin:endrowin,startcolin:endcolin);
void MatSubstitute(double matin[], double matout[], int matincolsize, int matoutcolsize, int startrowout, int endrowout, int startcolout, int endcolout, int startrowin, int endrowin, int startcolin, int endcolin)
{
    int rowsizeout = endrowout - startrowout;
    int colsizeout = endcolout - startcolout;
    for (int i = 0; i < rowsizeout + 1; i++)
        for (int j = 0; j < colsizeout + 1; j++)
            matout[(startrowout + i) * matoutcolsize + (startcolout + j)] = matin[(startrowin + i) * matincolsize + (startcolin + j)];
}

// 15x15 matrix * scalar
void MatScalarMult15(double matin[], double matout[], double scalar)
{
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++)
            matout[i * 15 + j] = matin[i * 15 + j] * scalar;
}

// 3x3 matrix * scalar
void MatScalarMult3(double matin[], double matout[], double scalar)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            matout[i * 3 + j] = matin[i * 3 + j] * scalar;
}

// general matrix * scalar
void MatScalarMult(double matin[], double matout[], double scalar, int row, int col)
{
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matout[i * col + j] = matin[i * col + j] * scalar;
}

void VecSubstitute(double *vin, double *vout, int size)
{
    for (int i = 0; i < size; i++)
        vout[i] = vin[i];
}

void VecScalarMult(double vin[], double scalarin, double vout[], int size)
{
    for (int i = 0; i < size; i++)
        vout[i] = vin[i] * scalarin;
}

void VecSum(double *vin, double *vin2, int size, double *vout)
{
    for (int i = 0; i < size; i++)
        vout[i] = vin[i] + vin2[i];
}

void VecSubtract(double vin[], double vin2[], int size, double vout[])
{
    for (int i = 0; i < size; i++)
        vout[i] = vin[i] - vin2[i];
}

void CrossProduct(double vin[], double vin2[], double vout[])
{
    vout[0] = vin[1] * vin2[2] - vin[2] * vin2[1];
    vout[1] = vin[2] * vin2[0] - vin[0] * vin2[2];
    vout[2] = vin[0] * vin2[1] - vin[1] * vin2[0];
}

void GetSkew(double vin[], double skewout[][3])
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
void InvMat22(double A[][2], double invA[][2])
{
    double detA = DetMat22(A);

    invA[0][0] = A[1][1] / detA;
    invA[0][1] = -A[0][1] / detA;
    invA[1][0] = -A[1][0] / detA;
    invA[1][1] = A[0][0] / detA;
}

//Inverse of 3x3 matrix
void InvMat33(double A[][3], double invA[][3])
{
    double detA = DetMat33(A);

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
double DetMat22(double A[][2])
{
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
}

//Determinant of 3x3 matrix
double DetMat33(double A[][3])
{
    return (A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1] -
            A[0][1] * A[1][0] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
            A[0][2] * A[1][0] * A[2][1] - A[0][2] * A[1][1] * A[2][0]);
}

#define DBL_EPSILON 1.0204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
// #define DBL_EPSILON 2.2204460492503131e-016 //가 원래값
#define POS_EPS DBL_EPSILON
#define mat(name, col_size, r, c) (*(name + r * col_size + c))

void inv_mat(double *a, double *result, int np)
{
    int i, j, ii, maxi;
    double p[15][15];
    double temp;
    //double EPS = POS_EPS;

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
            exit(1);
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
void eulr2dcm(double r, double p, double y, double Cbnout[][3])
{
    double Cpsi = cos(y);
    double Spsi = sin(y);
    double Cthe = cos(p);
    double Sthe = sin(p);
    double Cphi = cos(r);
    double Sphi = sin(r);

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

void eulr2dcmtr(double r, double p, double y, double Cnbout[][3])
{
    double Cpsi = cos(y);
    double Spsi = sin(y);
    double Cthe = cos(p);
    double Sthe = sin(p);
    double Cphi = cos(r);
    double Sphi = sin(r);

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
void eulr2qua(double r, double p, double y, double quatout[4])
{
    double cpsi2 = cos(y / 2);
    double spsi2 = sin(y / 2);
    double cthe2 = cos(p / 2);
    double sthe2 = sin(p / 2);
    double cphi2 = cos(r / 2);
    double sphi2 = sin(r / 2);
    quatout[0] = cphi2 * cthe2 * cpsi2 + sphi2 * sthe2 * spsi2;
    quatout[1] = sphi2 * cthe2 * cpsi2 - cphi2 * sthe2 * spsi2;
    quatout[2] = cphi2 * sthe2 * cpsi2 + sphi2 * cthe2 * spsi2;
    quatout[3] = cphi2 * cthe2 * spsi2 - sphi2 * sthe2 * cpsi2;
}

//qbn to Cbn
void qua2dcm(double qbnin[], double Cbnout[][3])
{
    double a = qbnin[0];
    double b = qbnin[1];
    double c = qbnin[2];
    double d = qbnin[3];

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

void dcm2eulr(double Cbnin[][3], double eulrout[])
{
    double roll = atan(Cbnin[2][1] / Cbnin[2][2]);
    double pitch = asin(-Cbnin[2][0]);
    double yaw = atan2(Cbnin[1][0], Cbnin[0][0]);
    double rpy[3] = {roll, pitch, yaw};
    VecSubstitute(rpy, eulrout, 3);
}

void QuatUpdate(double *wbnbin, double dtin, double *quatout)
{
    double dtheta[3];
    VecScalarMult(wbnbin, dtin, dtheta, 3);
    double Wx = dtheta[0];
    double Wy = dtheta[1];
    double Wz = dtheta[2];
    double theta = sqrt(Wx * Wx + Wy * Wy + Wz * Wz);
    double ac = cos(theta / 2);
    double as;
    if (theta < 1e-15)
        as = 0;
    else
        as = sin(theta / 2) / theta;
    double rk[4] = {ac, as * Wx, as * Wy, as * Wz};
    double Rk[4][4] = {{rk[0], -rk[1], -rk[2], -rk[3]}, {rk[1], rk[0], rk[3], -rk[2]}, {rk[2], -rk[3], rk[0], rk[1]}, {rk[3], rk[2], -rk[1], rk[0]}};
    double qtmp[4];
    MatMult41(Rk, quatout, qtmp);
    double mag = sqrt(pow(qtmp[0], 2) + pow(qtmp[1], 2) + pow(qtmp[2], 2) + pow(qtmp[3], 2));
    double q_normal = 1.5 - 0.5 * mag * mag;
    VecScalarMult(qtmp, q_normal, qtmp, 4);
    VecSubstitute(qtmp, quatout, 4);
}

void GenINSFmatrix15(double *llh, double *vn, double Cbn[][3], double *fn, double *wnin, double *wnie, double *wnen, double F[][15])
{
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 15; j++)
            F[i][j] = 0; // 안하면 쓰레기값이 들어가는 경우 발생
    static double Rm, Rt, Rmm, Rtt;
    Radicurv(llh[0], &Rm, &Rt);
    delRadicurv(llh[0], &Rmm, &Rtt);
    static double CL = cos(llh[0]);
    static double TL = tan(llh[0]);
    static double h = llh[2];
    static double vn_n = vn[0];
    static double vn_e = vn[1];
    static double vn_d = vn[2];
    static double wnie_n = wnie[0];
    static double wnie_d = wnie[2];
    static double wnen_n = wnen[0];
    static double wnen_e = wnen[1];
    static double wnen_d = wnen[2];
    static double F11[3][3] = {{(wnen_e * Rmm) / (Rm + h), 0, wnen_e / (Rm + h)},
                               {wnen_n * (TL - Rtt / (Rt + h)) / CL, 0, -wnen_n / (Rt + h) / CL},
                               {0, 0, 0}};
    MatSubstitute((double *)F11, (double *)F, 3, 15, 0, 2, 0, 2, 0, 2, 0, 2);

    static double F12[3][3] = {{1 / (Rm + h), 0, 0},
                               {0, 1 / (Rt + h) / CL, 0},
                               {0, 0, -1}};
    MatSubstitute((double *)F12, (double *)F, 3, 15, 0, 2, 3, 5, 0, 2, 0, 2);

    static double F21[3][3] = {{(-vn_e * (2 * wnie_n + wnen_n / CL / CL) + (vn_d * wnen_e * Rmm) / (Rm + h) - wnen_n * wnen_d * Rtt), 0, ((wnen_e * vn_d) / (Rm + h) - wnen_n * wnen_d)},
                               {-vn_d * (wnen_n * Rtt / (Rt + h) - 2 * wnie_d) + vn_n * (2 * wnie_n + wnen_n / CL / CL + wnen_d * Rtt / (Rt + h)), 0, vn_n * wnen_d / (Rt + h) - wnen_n * vn_d / (Rt + h)},
                               {-2 * wnie_d * vn_e + (pow(wnen_n, 2.)) * Rtt + (pow(wnen_e, 2.)) * Rmm, 0, pow(wnen_n, 2.) + pow(wnen_e, 2.)}};

    MatSubstitute((double *)F21, (double *)F, 3, 15, 3, 5, 0, 2, 0, 2, 0, 2);

    static double F221[3][3] = {{vn_d / (Rm + h), wnen_d, 0},
                                {0, (vn_d + vn_n * TL) / (Rt + h), 0},
                                {wnen_e, -wnen_n, 0}};
    static double F222[3][3] = {0};
    static double wnsubtract[3] = {0};
    static double wnietmp[3] = {0};
    VecScalarMult(wnie, -2, wnietmp, 3);
    MatSubtract(wnietmp, wnen, 3, 1, wnsubtract);
    GetSkew(wnsubtract, F222);
    static double F22[3][3] = {0};
    MatSum((double *)F221, (double *)F222, 3, 3, (double *)F22);
    MatSubstitute((double *)F22, (double *)F, 3, 15, 3, 5, 3, 5, 0, 2, 0, 2);

    static double F23[3][3] = {0};
    GetSkew(fn, F23);
    MatSubstitute((double *)F23, (double *)F, 3, 15, 3, 5, 6, 8, 0, 2, 0, 2);

    static double F24[3][3] = {0};
    MatSubstitute((double *)Cbn, (double *)F24, 3, 3, 0, 2, 0, 2, 0, 2, 0, 2);
    MatSubstitute((double *)F24, (double *)F, 3, 15, 3, 5, 9, 11, 0, 2, 0, 2);

    static double F31[3][3] = {{wnie_d - wnen_n * Rtt / (Rt + h), 0, -wnen_n / (Rt + h)},
                               {-wnen_e * Rmm / (Rm + h), 0, -wnen_e / (Rm + h)},
                               {-wnie_n - wnen_n / CL / CL - wnen_d * Rtt / (Rt + h), 0, -wnen_d / (Rt + h)}};
    MatSubstitute((double *)F31, (double *)F, 3, 15, 6, 8, 0, 2, 0, 2, 0, 2);

    static double F32[3][3] = {{0, 1 / (Rt + h), 0},
                               {-1 / (Rm + h), 0, 0},
                               {0, -TL / (Rt + h), 0}};
    MatSubstitute((double *)F32, (double *)F, 3, 15, 6, 8, 3, 5, 0, 2, 0, 2);

    static double F33[3][3] = {0};
    static double minuswnin[3] = {0};
    VecScalarMult(wnin, -1, minuswnin, 3);
    GetSkew(minuswnin, F33);
    MatSubstitute((double *)F33, (double *)F, 3, 15, 6, 8, 6, 8, 0, 2, 0, 2);

    static double F35[3][3] = {0};
    MatScalarMult3((double *)Cbn, (double *)F35, -1);
    MatSubstitute((double *)F35, (double *)F, 3, 15, 6, 8, 12, 14, 0, 2, 0, 2);
}

void InsKf15_predict(double F[][15], double *estX, double P_DB[][15], double *Qvec, double Cbn[][3], double Cbn_T[][3], double PHI[][15], double Pout[][15])
{
    // for (int i = 0; i < 15; i++)
    //     for (int j = 0; j < 15; j++)
    //         PHI[i][j] = 0; // 안하면 쓰레기값이 들어가는 경우 발생
    static double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static double eyemat[15][15], Qm[15][15] = {0};
    static double A[3][3], A_tmp[3][3], G[3][3], G_tmp[3][3], Qmat1[3][3], Qmat2[3][3] = {0};
    static double Qvec1[3] = {Qvec[0], Qvec[1], Qvec[2]};
    static double Qvec2[3] = {Qvec[3], Qvec[4], Qvec[5]};
    static double PHI_tmp1[15][15], PHI_tmp2[15][15], PHI_tmp3[15][15], PHI_tmp4[15][15] = {0};
    static double P_tmp1[15][15], P_tmp2[15][15], P_tmp3[15][15], P_tmp4[15][15] = {0};
    diag(eyevec, 15, (double *)eyemat);
    diag(Qvec1, 3, (double *)Qmat1);
    diag(Qvec2, 3, (double *)Qmat2);
    MatMult((double *)Cbn, (double *)Qmat1, (double *)A_tmp, 3, 3, 3);
    MatMult((double *)A_tmp, (double *)Cbn_T, (double *)A, 3, 3, 3);
    MatMult((double *)Cbn, (double *)Qmat2, (double *)G_tmp, 3, 3, 3);
    MatMult((double *)G_tmp, (double *)Cbn_T, (double *)G, 3, 3, 3);
    MatSubstitute((double *)A, (double *)Qm, 3, 15, 3, 5, 3, 5, 0, 2, 0, 2);
    MatSubstitute((double *)G, (double *)Qm, 3, 15, 6, 8, 6, 8, 0, 2, 0, 2);
    // discretize F to PHI
    static double F_tmp[15][15] = {0};
    MatSubstitute((double *)F, (double *)F_tmp, 15, 15, 0, 14, 0, 14, 0, 14, 0, 14);
    MatScalarMult15((double *)F, (double *)PHI_tmp1, dt);
    MatMult((double *)F, (double *)F_tmp, (double *)PHI_tmp2, 15, 15, 15);
    MatScalarMult15((double *)PHI_tmp2, (double *)PHI_tmp3, pow(dt, 2) / 2);
    MatSum((double *)eyemat, (double *)PHI_tmp1, 15, 15, (double *)PHI_tmp4);
    MatSum((double *)PHI_tmp4, (double *)PHI_tmp3, 15, 15, (double *)PHI);
    // Time update of error covariance
    MatMult((double *)PHI, (double *)P_DB, (double *)P_tmp1, 15, 15, 15);
    transpose((double *)PHI, (double *)P_tmp2, 15, 15);
    MatMult((double *)P_tmp1, (double *)P_tmp2, (double *)P_tmp3, 15, 15, 15);
    MatScalarMult15((double *)Qm, (double *)P_tmp4, dt);
    MatSum((double *)P_tmp3, (double *)P_tmp4, 15, 15, (double *)Pout);
}

void InsKf15_updatePV(double P_DB[][15], double *estX, double R[][6], double H[][15], double *measurement, double Pout[][15], double *updatedX)
{
    static double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static double eyemat[15][15] = {0};
    static double K[15][6], Ktmp1[15][6], K_T[6][15] = {0};
    static double HPb[6][15], H_T[15][6] = {0};
    static double invtmp1[6][6], invtmp2[6][6], inverse[6][6] = {0};
    static double eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    static double KRK[15][15], KRK_tmp1[15][6] = {0};
    static double updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    static double Hx[6], residual[6], updatedX_tmp[15] = {0};
    // calculating Kalman gain
    diag((double *)eyevec, 15, (double *)eyemat);
    MatMult((double *)H, (double *)P_DB, (double *)HPb, 6, 15, 15);
    transpose((double *)H, (double *)H_T, 6, 15);
    MatMult((double *)HPb, (double *)H_T, (double *)invtmp1, 6, 15, 6);
    MatSum((double *)invtmp1, (double *)R, 6, 6, (double *)invtmp2);
    inv_mat((double *)invtmp2, (double *)inverse, 6);
    MatMult((double *)P_DB, (double *)H_T, (double *)Ktmp1, 15, 15, 6);
    MatMult((double *)Ktmp1, (double *)inverse, (double *)K, 15, 6, 6);

    // calculating updated P
    MatMult((double *)K, (double *)H, (double *)eyeminusKH_tmp, 15, 6, 15);
    MatSubtract((double *)eyemat, (double *)eyeminusKH_tmp, 15, 15, (double *)eyeminusKH);
    transpose((double *)eyeminusKH, (double *)eyeminusKH_T, 15, 15);
    transpose((double *)K, (double *)K_T, 15, 6);
    MatMult((double *)K, (double *)R, (double *)KRK_tmp1, 15, 6, 6);
    MatMult((double *)KRK_tmp1, (double *)K_T, (double *)KRK, 15, 6, 15);
    MatMult((double *)eyeminusKH, (double *)P_DB, (double *)updatedP_tmp, 15, 15, 15);
    MatMult((double *)updatedP_tmp, (double *)eyeminusKH_T, (double *)updatedP_tmp2, 15, 15, 15);
    MatSum((double *)updatedP_tmp2, (double *)KRK, 15, 15, (double *)Pout);

    // calculating updated state
    MatMult((double *)H, (double *)estX, (double *)Hx, 6, 15, 1);
    MatSubtract((double *)measurement, (double *)Hx, 6, 1, (double *)residual);
    MatMult((double *)K, (double *)residual, (double *)updatedX_tmp, 15, 6, 1);
    VecSum((double *)estX, (double *)updatedX_tmp, 15, (double *)updatedX);
}

void InsKf15_updateV(double P_DB[][15], double *estX, double R[][3], double H[][15], double *measurement, double Pout[][15], double *updatedX)
{
    static double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static double eyemat[15][15] = {0};
    static double K[15][3], Ktmp1[15][3], K_T[3][15] = {0};
    static double HPb[3][15], H_T[15][3] = {0};
    static double invtmp1[3][3], invtmp2[3][3], inverse[3][3] = {0};
    static double eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    static double KRK[15][15], KRK_tmp1[15][3] = {0};
    static double updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    static double Hx[3], residual[3], updatedX_tmp[15] = {0};
    // calculating Kalman gain
    diag((double *)eyevec, 15, (double *)eyemat);
    MatMult((double *)H, (double *)P_DB, (double *)HPb, 3, 15, 15);
    transpose((double *)H, (double *)H_T, 3, 15);
    MatMult((double *)HPb, (double *)H_T, (double *)invtmp1, 3, 15, 3);
    MatSum((double *)invtmp1, (double *)R, 3, 3, (double *)invtmp2);
    inv_mat((double *)invtmp2, (double *)inverse, 3);
    MatMult((double *)P_DB, (double *)H_T, (double *)Ktmp1, 15, 15, 3);
    MatMult((double *)Ktmp1, (double *)inverse, (double *)K, 15, 3, 3);

    // calculating updated P
    MatMult((double *)K, (double *)H, (double *)eyeminusKH_tmp, 15, 3, 15);
    MatSubtract((double *)eyemat, (double *)eyeminusKH_tmp, 15, 15, (double *)eyeminusKH);
    transpose((double *)eyeminusKH, (double *)eyeminusKH_T, 15, 15);
    transpose((double *)K, (double *)K_T, 15, 3);
    MatMult((double *)K, (double *)R, (double *)KRK_tmp1, 15, 3, 3);
    MatMult((double *)KRK_tmp1, (double *)K_T, (double *)KRK, 15, 3, 15);
    MatMult((double *)eyeminusKH, (double *)P_DB, (double *)updatedP_tmp, 15, 15, 15);
    MatMult((double *)updatedP_tmp, (double *)eyeminusKH_T, (double *)updatedP_tmp2, 15, 15, 15);
    MatSum((double *)updatedP_tmp2, (double *)KRK, 15, 15, (double *)Pout);

    // calculating updated state
    MatMult((double *)H, (double *)estX, (double *)Hx, 3, 15, 1);
    MatSubtract((double *)measurement, (double *)Hx, 3, 1, (double *)residual);
    MatMult((double *)K, (double *)residual, (double *)updatedX_tmp, 15, 3, 1);
    VecSum((double *)estX, (double *)updatedX_tmp, 15, (double *)updatedX);
}

void InsKf15_updateHgt(double P_DB[][15], double *estX, double R[][1], double H[][15], double *measurement, double Pout[][15], double *updatedX)
{
    static double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static double eyemat[15][15] = {0};
    static double K[15][1], Ktmp1[15][1], K_T[1][15] = {0};
    static double HPb[1][15], H_T[15][1] = {0};
    static double invtmp1[1][1], invtmp2[1][1], inverse[1][1] = {0};
    static double eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    static double KRK[15][15], KRK_tmp1[15][1] = {0};
    static double updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    static double Hx[1], residual[1], updatedX_tmp[15] = {0};
    // calculating Kalman gain
    diag((double *)eyevec, 15, (double *)eyemat);
    MatMult((double *)H, (double *)P_DB, (double *)HPb, 1, 15, 15);
    transpose((double *)H, (double *)H_T, 1, 15);
    MatMult((double *)HPb, (double *)H_T, (double *)invtmp1, 1, 15, 1);
    MatSum((double *)invtmp1, (double *)R, 1, 1, (double *)invtmp2);
    inv_mat((double *)invtmp2, (double *)inverse, 1);
    MatMult((double *)P_DB, (double *)H_T, (double *)Ktmp1, 15, 15, 1);
    MatMult((double *)Ktmp1, (double *)inverse, (double *)K, 15, 1, 1);

    // calculating updated P
    MatMult((double *)K, (double *)H, (double *)eyeminusKH_tmp, 15, 1, 15);
    MatSubtract((double *)eyemat, (double *)eyeminusKH_tmp, 15, 15, (double *)eyeminusKH);
    transpose((double *)eyeminusKH, (double *)eyeminusKH_T, 15, 15);
    transpose((double *)K, (double *)K_T, 15, 1);
    MatMult((double *)K, (double *)R, (double *)KRK_tmp1, 15, 1, 1);
    MatMult((double *)KRK_tmp1, (double *)K_T, (double *)KRK, 15, 1, 15);
    MatMult((double *)eyeminusKH, (double *)P_DB, (double *)updatedP_tmp, 15, 15, 15);
    MatMult((double *)updatedP_tmp, (double *)eyeminusKH_T, (double *)updatedP_tmp2, 15, 15, 15);
    MatSum((double *)updatedP_tmp2, (double *)KRK, 15, 15, (double *)Pout);

    // calculating updated state
    MatMult((double *)H, (double *)estX, (double *)Hx, 1, 15, 1);
    MatSubtract((double *)measurement, (double *)Hx, 1, 1, (double *)residual);
    MatMult((double *)K, (double *)residual, (double *)updatedX_tmp, 15, 1, 1);
    VecSum((double *)estX, (double *)updatedX_tmp, 15, (double *)updatedX);
}

void InsKf15_updatelatlon(double P_DB[][15], double *estX, double R[][2], double H[][15], double *measurement, double Pout[][15], double *updatedX)
{
    static double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static double eyemat[15][15] = {0};
    static double K[15][2], Ktmp1[15][2], K_T[2][15] = {0};
    static double HPb[2][15], H_T[15][2] = {0};
    static double invtmp1[2][2], invtmp2[2][2], inverse[2][2] = {0};
    static double eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    static double KRK[15][15], KRK_tmp1[15][2] = {0};
    static double updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    static double Hx[2], residual[2], updatedX_tmp[15] = {0};
    // calculating Kalman gain
    diag((double *)eyevec, 15, (double *)eyemat);
    MatMult((double *)H, (double *)P_DB, (double *)HPb, 2, 15, 15);
    transpose((double *)H, (double *)H_T, 2, 15);
    MatMult((double *)HPb, (double *)H_T, (double *)invtmp1, 2, 15, 2);
    MatSum((double *)invtmp1, (double *)R, 2, 2, (double *)invtmp2);
    inv_mat((double *)invtmp2, (double *)inverse, 2);
    MatMult((double *)P_DB, (double *)H_T, (double *)Ktmp1, 15, 15, 2);
    MatMult((double *)Ktmp1, (double *)inverse, (double *)K, 15, 2, 2);

    // calculating updated P
    MatMult((double *)K, (double *)H, (double *)eyeminusKH_tmp, 15, 2, 15);
    MatSubtract((double *)eyemat, (double *)eyeminusKH_tmp, 15, 15, (double *)eyeminusKH);
    transpose((double *)eyeminusKH, (double *)eyeminusKH_T, 15, 15);
    transpose((double *)K, (double *)K_T, 15, 2);
    MatMult((double *)K, (double *)R, (double *)KRK_tmp1, 15, 2, 2);
    MatMult((double *)KRK_tmp1, (double *)K_T, (double *)KRK, 15, 2, 15);
    MatMult((double *)eyeminusKH, (double *)P_DB, (double *)updatedP_tmp, 15, 15, 15);
    MatMult((double *)updatedP_tmp, (double *)eyeminusKH_T, (double *)updatedP_tmp2, 15, 15, 15);
    MatSum((double *)updatedP_tmp2, (double *)KRK, 15, 15, (double *)Pout);

    // calculating updated state
    MatMult((double *)H, (double *)estX, (double *)Hx, 2, 15, 1);
    MatSubtract((double *)measurement, (double *)Hx, 2, 1, (double *)residual);
    MatMult((double *)K, (double *)residual, (double *)updatedX_tmp, 15, 2, 1);
    VecSum((double *)estX, (double *)updatedX_tmp, 15, (double *)updatedX);
}

void InsKf15_updatelatlonvnve(double P_DB[][15], double *estX, double R[][4], double H[][15], double *measurement, double Pout[][15], double *updatedX) ////////////////////**
{
    double eyevec[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double eyemat[15][15] = {0};
    double K[15][4], Ktmp1[15][4], K_T[4][15] = {0};
    double HPb[4][15], H_T[15][4] = {0};
    double invtmp1[4][4], invtmp2[4][4], inverse[4][4] = {0};
    double eyeminusKH[15][15], eyeminusKH_tmp[15][15], eyeminusKH_T[15][15] = {0};
    double KRK[15][15], KRK_tmp1[15][4] = {0};
    double updatedP_tmp[15][15], updatedP_tmp2[15][15] = {0};
    double Hx[4], residual[4], updatedX_tmp[15] = {0};

    // calculating Kalman gain
    diag((double *)eyevec, 15, (double *)eyemat);
    MatMult((double *)H, (double *)P_DB, (double *)HPb, 4, 15, 15);
    transpose((double *)H, (double *)H_T, 4, 15);
    MatMult((double *)HPb, (double *)H_T, (double *)invtmp1, 4, 15, 4);
    MatSum((double *)invtmp1, (double *)R, 4, 4, (double *)invtmp2);
    inv_mat((double *)invtmp2, (double *)inverse, 4);
    MatMult((double *)P_DB, (double *)H_T, (double *)Ktmp1, 15, 15, 4);
    MatMult((double *)Ktmp1, (double *)inverse, (double *)K, 15, 4, 4);

    // calculating updated P
    MatMult((double *)K, (double *)H, (double *)eyeminusKH_tmp, 15, 4, 15);
    MatSubtract((double *)eyemat, (double *)eyeminusKH_tmp, 15, 15, (double *)eyeminusKH);
    transpose((double *)eyeminusKH, (double *)eyeminusKH_T, 15, 15);
    transpose((double *)K, (double *)K_T, 15, 4);
    MatMult((double *)K, (double *)R, (double *)KRK_tmp1, 15, 4, 4);
    MatMult((double *)KRK_tmp1, (double *)K_T, (double *)KRK, 15, 4, 15);
    MatMult((double *)eyeminusKH, (double *)P_DB, (double *)updatedP_tmp, 15, 15, 15);
    MatMult((double *)updatedP_tmp, (double *)eyeminusKH_T, (double *)updatedP_tmp2, 15, 15, 15);
    MatSum((double *)updatedP_tmp2, (double *)KRK, 15, 15, (double *)Pout);

    // calculating updated state
    MatMult((double *)H, (double *)estX, (double *)Hx, 4, 15, 1);
    MatSubtract((double *)measurement, (double *)Hx, 4, 1, (double *)residual);
    MatMult((double *)K, (double *)residual, (double *)updatedX_tmp, 15, 4, 1);
    VecSum((double *)estX, (double *)updatedX_tmp, 15, (double *)updatedX);
}

void Correction(double *estX, double *llh, double *vn, double *quat, double Cbn[][3], double *ba, double *bg)
{
    // bias 정의를 반대로 해줬기 때문에(estimated F,Wbib = true - bias -> true = estimated + bias) correction도 음수 붙여서
    estX[9] = -estX[9];
    estX[10] = -estX[10];
    estX[11] = -estX[11];
    estX[12] = -estX[12];
    estX[13] = -estX[13];
    estX[14] = -estX[14];
    double estPos[3] = {estX[0], estX[1], estX[2]};
    double estVel[3] = {estX[3], estX[4], estX[5]};
    double estAtt[3] = {estX[6], estX[7], estX[8]};
    double estBa[3] = {estX[9], estX[10], estX[11]};
    double estBg[3] = {estX[12], estX[13], estX[14]};
    VecSubtract(llh, estPos, 3, llh);
    VecSubtract(vn, estVel, 3, vn);
    VecSubtract(ba, estBa, 3, ba);
    VecSubtract(bg, estBg, 3, bg);
    double Q_tmp[4][3] = {{-quat[1], -quat[2], -quat[3]},
                          {quat[0], quat[3], -quat[2]},
                          {-quat[3], quat[0], quat[1]},
                          {quat[2], -quat[1], quat[0]}};
    double del_Qtmp[4][3], del_Q[4] = {0};
    MatScalarMult((double *)Q_tmp, (double *)del_Qtmp, -0.5, 4, 3);
    MatMult((double *)del_Qtmp, (double *)estAtt, (double *)del_Q, 4, 3, 1);
    VecSubtract((double *)quat, (double *)del_Q, 4, (double *)quat);
    double qn = 1.5 - 0.5 * (pow(quat[0], 2) + pow(quat[1], 2) + pow(quat[2], 2) + pow(quat[3], 2));
    VecScalarMult(quat, qn, quat, 4);
    qua2dcm(quat, Cbn);
    VecScalarMult(estX, 0, estX, 15); // correction 이후 추정치를 기반으로 선형화 하기 때문에 오차가 없는 것으로 간주(EKF)
}

void Radicurv(double latin, double *rmout, double *rpout)
{
    double e2 = pow(e, 2);
    double den = 1 - e2 * pow(sin(latin), 2);
    *rmout = R0 * (1 - e2) / pow(den, 1.5);
    *rpout = R0 / sqrt(den);
}

void delRadicurv(double latin, double *rmmout, double *rttout)
{
    double e2 = pow(e, 2);
    double den = 1 - e2 * pow(sin(latin), 2);
    *rmmout = (3 * R0 * (1 - e2) * e2 * sin(latin) * cos(latin)) / (pow(den, 2.5));
    *rttout = (R0 * e2 * sin(latin) * cos(latin)) / (pow(den, 1.5));
}

void NavRotation(double *vin, double *llhin, double *wninout, double *wnieout, double *wnenout)
{
    double rm;
    double rp;
    Radicurv(llhin[0], &rm, &rp);
    double dlat = vin[0] / (rm + llhin[2]);
    double dlon = vin[1] / (rp + llhin[2]) / cos(llhin[0]);
    double dh = -vin[2];
    double clat = cos(llhin[0]);
    double slat = sin(llhin[0]);
    wnieout[0] = We * clat;
    wnieout[1] = 0;
    wnieout[2] = -We * slat;
    wnenout[0] = dlon * clat;
    wnenout[1] = -dlat;
    wnenout[2] = -dlon * slat;
    for (int i = 0; i < 3; i++)
        wninout[i] = wnieout[i] + wnenout[i];
}

double gravity_wgs84(double lat, double h)
{
    double a = 6378137.0;          // Semi - major axis[m] : equatorial radius of the earth; WGS - 84
    double b = 6356752.3142;       // Semi - minor axis[m]
    double ee = 0.081819190842622; // eccentricity of the earth ellipsoid

    double ER = 7.292115e-5;     // earth's rotation rate [rad/sec]
    double GM = 3.986004418e14;  // earth's gravitational constant [m3/s2]
    double ge = 9.7803253359;    // Theoretical gravity at the Equator[m / s2]
    double gp = 9.8321849378;    // Theoretical gravity at the Pole[m / s2]
    double k = 0.00193185265241; // Theoretical gravity formula constant
    double m = 0.00344978650684; // m = ER ^ 2 * a ^ 2 * b / GM;

    double e2 = pow(ee, 2.0);
    double a2 = pow(a, 2.0);
    double slat = sin(lat);
    double slat2 = pow(slat, 2.0);
    double slat4 = pow(slat2, 2.0);
    double h2 = pow(h, 2.0);

    double g1 = k + e2 / 2;
    double g2 = -2.0 * (1.0 + f + m) / a;
    double g3 = g1 * g2 + 4.0 * f / a;
    double g4 = 3.0 / a2;
    double g5 = g1 * g4;
    double g6 = e2 * (3.0 * e2 / 4.0 + k) / 2.0;
    double g7 = 4.0 * g1 * f / a + g2 * g6;
    double g8 = g4 * g6;

    return ge * (1.0 + g1 * slat2 + g2 * h + g3 * slat2 * h + g4 * h2 + g5 * slat2 * h2);
}

dynamicArray::dynamicArray() : m_size(1), m_used(0)
{
    m_array = new double[m_size];
}
dynamicArray::dynamicArray(int size) : m_size(size), m_used(0)
{
    m_array = new double[m_size];
}
dynamicArray::~dynamicArray()
{
    delete[] m_array;
}
void dynamicArray::pushback(double data)
{
    if (this->m_used < this->m_size)
    {
        this->m_array[this->m_used] = data;
    }
    else
    {
        double *temp = new double[m_size];
        for (int i = 0; i < this->m_used; i++)
            temp[i] = this->m_array[i];
        delete[] this->m_array; // 원래남아있는 거 지움
        this->m_size *= 2;
        this->m_array = new double[m_size];
        for (int i = 0; i < this->m_used; i++)
        {
            this->m_array[i] = temp[i];
        }
        delete[] temp;
        this->m_array[this->m_used] = data;
    }
    this->m_used++;
}

void dynamicArray::popback()
{
    if (m_used == 0)
        return;
    m_used--;
    double *temp = new double[m_used];
    for (int i = 0; i < m_used; i++)
        temp[i] = m_array[i];
    delete[] m_array;
    m_array = new double[m_size];
    for (int i = 0; i < m_used; i++)
        m_array[i] = temp[i];
    delete[] temp;
    return;
}

double &dynamicArray::operator[](int index)
{
    return this->m_array[index];
}
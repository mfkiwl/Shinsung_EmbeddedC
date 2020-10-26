#include "Initialize.h"
#include "common.h"
#include "function.h"
#include "ReadIMU.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "math.h"

using namespace std;
int main()
{
    // Initializing(filter) P, Q, R, H ,F declaration
    double pos_err_deviation[3] = {km2rad(5 * 1e-3), km2rad(5 * 1e-3), 10};
    double vel_err_deviation[3] = {noise_v, noise_v, noise_v};
    double att_err_deviation[3] = {a_repeatability / 9.81, a_repeatability / 9.81, g_repeatability * mg / We / cos(36 * d2r)}; //07.27.19:26 debug complete.
    double hgt_err_deviation = noise_baro_h;
    double Ppos[3];
    double Pvel[3];
    double Patt[3];
    double PBa[3] = {pow(a_repeatability, 2), pow(a_repeatability, 2), pow(a_repeatability, 2)};
    double PBg[3] = {pow(g_repeatability, 2), pow(g_repeatability, 2), pow(g_repeatability, 2)};
    double Qwa[3] = {pow(a_ARW, 2), pow(a_ARW, 2), pow(a_ARW, 2)};
    double Qwg[3] = {pow(g_ARW, 2), pow(g_ARW, 2), pow(g_ARW, 2)};
    // VecScalarMult(Qwg, 0.1, Qwg, 3); // tuning parameter
    double Qvec[6] = {Qwa[0], Qwa[1], Qwa[2], Qwg[0], Qwg[1], Qwg[2]};
    double P_DB[15][15] = {0};
    double H_latlon[2][15], H_NHC[2][15], H_PV[6][15], H_Zupt[3][15], H_Zupt3[5][15], H_yaw[1][15], H_hgt[1][15], H_latlonVnVe[4][15] = {0}; // 위도경도가 측정치일 때, 위치속도가 측정치일 때, Zupt, 고도가 측정치일 때 자코비안
    double gpsvel_N = 0;
    double gpsvel_E = 0;
    double subcount = 0;
    H_latlon[0][0] = 1;
    H_latlon[1][1] = 1;
    H_NHC[0][0] = 1;
    H_NHC[1][1] = 1;
    H_hgt[0][2] = 1;
    H_yaw[0][8] = -1;
    H_PV[0][0] = 1;
    H_PV[1][1] = 1;
    H_PV[2][2] = 1;
    H_PV[3][3] = 1;
    H_PV[4][4] = 1;
    H_PV[5][5] = 1;
    H_latlonVnVe[0][0] = 1;
    H_latlonVnVe[1][1] = 1;
    H_latlonVnVe[2][3] = 1;
    H_latlonVnVe[3][4] = 1;
    H_Zupt[0][3] = 1;
    H_Zupt[1][4] = 1;
    H_Zupt[2][5] = 1;
    H_Zupt3[0][0] = 1;
    H_Zupt3[1][1] = 1;
    H_Zupt3[2][3] = 1;
    H_Zupt3[3][4] = 1;
    H_Zupt3[4][5] = 1;
    double R_latlon[2][2];
    double R_NHC[2][2];
    double R_latlonvnve[4][4];
    double R_PV[6][6];
    double R_vec[2] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000)};
    double R_hgt[1][1] = {pow(noise_baro_h, 2)};
    double R_yaw[1][1] = {pow((0.5 * d2r), 2)};
    double R_zupt[3][3] = {0};
    double R_zupt2[3][3] = {0};
    double R_zupt3[5][5] = {0};
    double F[15][15], PHI[15][15] = {0};
    double R_vec2[6] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_vp, noise_v, noise_v, noise_v};
    double R_vec3[4] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_v, noise_v};
    double R_vec4[3] = {noise_v, noise_v, noise_v};
    double R_vec5[2] = {noise_v, noise_v};
    double R_vec6[3] = {noise_v * 2, noise_v * 2, noise_v * 2};
    double R_vec7[5] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_v, noise_v, noise_v};
    Velwisesquare(R_vec2, 6, R_vec2);
    Velwisesquare(R_vec, 2, R_vec);
    Velwisesquare(R_vec3, 4, R_vec3);
    Velwisesquare(R_vec4, 3, R_vec4);
    Velwisesquare(R_vec5, 2, R_vec5);
    Velwisesquare(R_vec6, 3, R_vec6);
    Velwisesquare(R_vec7, 5, R_vec7);
    diag(R_vec, 2, (double *)R_latlon);
    diag(R_vec2, 6, (double *)R_PV);
    diag(R_vec3, 4, (double *)R_latlonvnve);
    diag(R_vec4, 3, (double *)R_zupt);
    diag(R_vec5, 2, (double *)R_NHC);
    diag(R_vec6, 3, (double *)R_zupt2);
    diag(R_vec7, 5, (double *)R_zupt3);
    Velwisesquare(pos_err_deviation, 3, Ppos);
    Velwisesquare(vel_err_deviation, 3, Pvel);
    Velwisesquare(att_err_deviation, 3, Patt);
    double P_INS[15] = {Ppos[0], Ppos[1], Ppos[2], Pvel[0], Pvel[1], Pvel[2], Patt[0], Patt[1], Patt[2], PBa[0], PBa[1], PBa[2], PBg[0], PBg[1], PBg[2]};
    diag(P_INS, 15, (double *)P_DB);
    //
    // Initializing (position, velocity, attitude)
    double INSRadPos[3], INSDegPos[3], INSRadEulr[3], INSDegEulr[3], INSQT[4], INSDCMCnb[3][3], INSDCMCbn[3][3], INSMpsVned[3], INSMpsVb[3] = {0};
    double Ba_init[3], Bg_init[3], Ba[3], Bg[3], wnin[3], wnie[3], wnen[3], wbin[3], wbnb[3] = {0};
    double estX[15] = {0};                                                   //추정 오차
    double Fb[3], Wbib[3], fn[3] = {0};                                      // fb, wbib에서 바이어스를 뺀 값, Fb에 Cbn곱한 값
    double MeasurementPV[6], Measurementlatlon[2] = {0}, MeasurementHgt = 0; // 위치,속도만 사용하는 경우, 위도,경도만 사용하는 경우, 고도만 사용하는 경우 EKF의 측정치 (delta 측정치)
    double Measurementlatlonvnve[4] = {0};
    double MeasurementNHC[2] = {0};
    double MeasurementZupt[5] = {0};
    double MeasurementV[3] = {0};
    double MeasurementYaw = 0;
    double g[3] = {0}; // 중력 벡터
    double rm, rp = 0; // meridian radius of curvature, transverse radius of curvature
    double firsttime = 0;
    //
    // read data
    CheckInputFile();
    ifstream IMUGPSfile; // IMU와 GPS정보가 같이 들어오는 파일
    // ifstream RTKGPSfile;      // RTK(reference)정보가 들어오는 파일
    string str, str2;
    double validrow = 0;                                                                         // GPS 첫 데이터가 들어온 시점부터 정렬을 시작하기위해 invalidrow를 정의, validrow는 정렬을 시작하는 첫번째 열을 뜻함.
    bool dataisincorrect = true;                                                                 // GPS 데이터 기준으로 사용 가능한 데이터가 들어오는지 판단
    bool aligned = false;                                                                        // 정렬 여부 판단
    double GPS_count = 0;                                                                        // 정렬 시간 까지 유효한 GPS데이터 개수 세는데 사용
    int linenum = 0;                                                                             // 한 epoch마다 count
    double time, fb[3], wbib[3], gps_data[7], Baro_data;                                         // Data save용 변수, previous gps 는 zupt 실행을 위해 직전 gps data와 비교
    double mean_acc[3], mean_gyr[3] = {0}, mean_BaroOffset, roll_init, pitch_init, yaw_init = 0; // 정렬에 사용되는 변수, 기압계 데이터가 전체적으로 bias되어 있어 이를 고려함.
    double gps_tmp, baro_tmp = 0;                                                                // mean_BaroOffset을 계산하기 위해 사용되는 변수
    mean_BaroOffset = 0;
    vector<double> parsedIMUGPS;          // 동적 array, GPS 데이터 섞일 때마다 데이터의 길이가 달라지므로 사용함
    vector<double> parsedPrevious(17, 0); // Zupt실행을 위해 이전 GPS값끼리 비교하기 위한 변수, 0으로 초기화해서 정의.
    // IMUGPSfile.open("../imu20201016_v2.txt");
    // IMUGPSfile.open("../imu_V2_201021.txt");
    IMUGPSfile.open("../imu20201023_BMI.txt");
    // IMUGPSfile.open("../parsed.txt");
    // IMUGPSfile.open("../imu20200923.txt");
    // RTKGPSfile.open("20191023_RTKGPS.txt");
    //
    // write data
    ofstream EstResult("../result_posvel.txt");
    // ofstream myGPS("../mygps.txt");
    // ofstream myhgt("../altimeter.txt");
    // ofstream myvel("../myvel.txt");
    if (!EstResult.is_open())
        cout << "output file is not opened!" << endl;
    if (IMUGPSfile.is_open())
    {
        while (getline(IMUGPSfile, str))
        {
            linenum++; //1부터 셈
            if ((linenum % 3) == 1)
                ReadData(str, parsedIMUGPS);
            else
                continue;
            if (dataisincorrect) // IMU 정보가 GPS정보와 같이 들어온 최초의 경우 GPS정보가 부정확할 수 있음, 이를 걸러냄
            {
                if (parsedIMUGPS.size() == 17)
                {
                    if (floor(parsedIMUGPS[11]) == 37)
                    {
                        dataisincorrect = false;
                        cout << "The first non-trash value is at " << linenum << "th row" << endl;
                        firsttime = parsedIMUGPS[0];
                        validrow = linenum;
                    }
                }
                else
                {
                    cout << "Passing trash value" << endl;
                }
            }
            if (!dataisincorrect) //첫줄 부터 정확한 신호가 들어오는 경우가 존재할 수 있으므로 if문 필요
            {
                time = parsedIMUGPS[0];    // 시간
                wbib[0] = parsedIMUGPS[1]; // 자이로스코프 입력
                wbib[1] = parsedIMUGPS[2];
                wbib[2] = parsedIMUGPS[3];
                fb[0] = parsedIMUGPS[4]; // 가속도계 입력
                fb[1] = parsedIMUGPS[5];
                fb[2] = parsedIMUGPS[6];
                Baro_data = parsedIMUGPS[10]; // 기압계 입력
                // myhgt << Baro_data << endl;
                if (parsedIMUGPS.size() == 17)
                {
                    gps_data[0] = parsedIMUGPS[0];                  // GPS 시간
                    gps_data[1] = parsedIMUGPS[11] * d2r;           // 위도
                    gps_data[2] = parsedIMUGPS[12] * d2r;           // 경도
                    gps_data[3] = parsedIMUGPS[13];                 // 고도
                    gps_data[4] = parsedIMUGPS[14] * (1.852 / 3.6); // speed
                    gps_data[5] = parsedIMUGPS[15];                 // 가시위성 개수
                    gps_data[6] = parsedIMUGPS[16];                 // 방향각
                }
                gpsvel_N = gps_data[4] * cos(gps_data[6] * d2r);
                gpsvel_E = gps_data[4] * sin(gps_data[6] * d2r);
                // myvel << gpsvel_N << " " << gpsvel_E << endl;
                if (!aligned)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        mean_acc[i] += fb[i];
                        mean_gyr[i] += wbib[i];
                    }
                    if (parsedIMUGPS.size() == 17)
                    {
                        gps_tmp += gps_data[3];
                        baro_tmp += Baro_data;
                        GPS_count++;
                    }
                    else
                        baro_tmp += Baro_data;
                    // if (((linenum - 1) / 3 + 1) - validrow == sampleHZ * t_align) // 유효한 데이터 받고나서 t_align만큼 정렬
                    if (gps_data[0] >= firsttime + t_align) // 유효한 데이터 받고나서 t_align만큼 정렬  10.21.11:13
                    {
                        double align_count = (linenum + 2) / 3 - (validrow - 1); // ex)300초 데이터(25*300 = 7500개) 에 쓰레기 데이터가 1~50열 까지 있으면 정렬에 사용된 데이터는 51~7501 열(7501열에 GPS들어오는 첫 순간까지 정렬로 사용).
                        // double align_count = 1;
                        for (int i = 0; i < 3; i++)
                        {
                            mean_acc[i] /= align_count;
                            mean_gyr[i] /= align_count;
                            Bg_init[i] = mean_gyr[i]; // 초기 바이어스를 정렬값으로 설정
                            mean_BaroOffset = gps_tmp / GPS_count - baro_tmp / align_count;
                        }
                        cout << "\n\n***Alignment Data***" << endl;
                        cout << align_count << "number of IMU data is used to align." << endl;
                        cout << GPS_count << "number of GPS data is used to calculate altitude offset between GPS and Barometer." << endl;
                        cout << "gyr_x: " << mean_gyr[0] << "rad/s\tgyr_y: " << mean_gyr[1] << "rad/s\tgyr_z: " << mean_gyr[2] << "rad/s" << endl;
                        cout << "acc_x: " << mean_acc[0] << "m/s^2\tacc_y: " << mean_acc[1] << "m/s^2\tacc_z: " << mean_acc[2] << "m/s^2" << endl;
                        roll_init = atan2(-mean_acc[1], -mean_acc[2]);
                        pitch_init = atan2(mean_acc[0], sqrt(pow(mean_acc[1], 2) + pow(mean_acc[2], 2)));
                        double C1T[3][3] = {0};
                        double C1T_T[3][3] = {0};
                        eulr2dcm(roll_init, pitch_init, 0, C1T);
                        transpose((double *)C1T, (double *)C1T_T, 3, 3);
                        double omega_tmp[3] = {0};
                        MatMult31(C1T_T, mean_gyr, omega_tmp);
                        yaw_init = atan2(-omega_tmp[1], omega_tmp[0]);
                        // yaw_init = 60 * d2r;
                        // 추정치 초기화
                        INSRadEulr[0] = roll_init;
                        INSRadEulr[1] = pitch_init;
                        INSRadEulr[2] = yaw_init; //기존코드: yaw init을 사용하지 않고 나침판어플(?)로 측정한 heading 사용함
                        // INSRadEulr[2] = 8 * d2r;
                        for (int i = 0; i < 3; i++)
                        {
                            INSRadPos[i] = gps_data[i + 1];
                            if (i < 2)
                                INSDegPos[i] = INSRadPos[i] * r2d; //높이값은 r2d곱하면 안됨.
                            else
                                INSDegPos[i] = INSRadPos[i];
                            INSDegEulr[i] = INSRadEulr[i] * r2d;
                            INSMpsVned[i] = 0;
                        }
                        VecSubstitute(Ba_init, Ba, 3);
                        VecSubstitute(Bg_init, Bg, 3);
                        // // /* For specific experiment
                        // double bgtmp[3] = {0.088926 * d2r, 0.15873 * d2r, 0.084921 * d2r};
                        // VecSubstitute(bgtmp, Bg, 3);
                        // eulr2qua(roll_init, pitch_init, 8 * d2r, INSQT);
                        // eulr2dcm(roll_init, pitch_init, 8 * d2r, INSDCMCbn);
                        // // */
                        eulr2qua(roll_init, pitch_init, yaw_init, INSQT);
                        eulr2dcm(roll_init, pitch_init, yaw_init, INSDCMCbn);
                        transpose((double *)INSDCMCbn, (double *)INSDCMCnb, 3, 3); // 대칭성 확보를 위한 trick
                        MatMult31(INSDCMCnb, INSMpsVned, INSMpsVb);
                        double zero6[6] = {0};
                        EstResult << setprecision(15) << fixed << scientific << INSRadPos[0] << " " << INSRadPos[1] << " " << INSRadPos[2] << " "
                                  << INSMpsVned[0] << " " << INSMpsVned[1] << " " << INSMpsVned[2] << " "
                                  << INSRadEulr[0] << " " << INSRadEulr[1] << " " << INSRadEulr[2] << " " << zero6[0] << " "
                                  << zero6[1] << " " << zero6[2] << " " << zero6[3] << " " << zero6[4] << " " << zero6[5] << endl;
                        aligned = true;
                    }
                }
                if (aligned)
                {
                    // Inertial Navigation
                    // Attitude
                    MatSubtract(fb, Ba, 3, 1, Fb);
                    MatSubtract(wbib, Bg, 3, 1, Wbib);
                    NavRotation(INSMpsVned, INSRadPos, wnin, wnie, wnen);
                    MatMult31(INSDCMCnb, wnin, wbin);
                    MatSubtract(Wbib, wbin, 3, 1, wbnb);
                    // if (((linenum - 1) / 3 + 1) - validrow <= t_finealign * sampleHZ)
                    // if (gps_data[4] <= 0.2 || ((linenum - 1) / 3 + 1) - validrow <= t_finealign * sampleHZ)
                    if (gps_data[0] <= firsttime + t_finealign)
                    {
                    }
                    // QuatUpdate(wbnb, 0, INSQT);
                    else
                        QuatUpdate(wbnb, dt, INSQT);
                    qua2dcm(INSQT, INSDCMCbn);
                    dcm2eulr(INSDCMCbn, INSRadEulr);
                    transpose((double *)INSDCMCbn, (double *)INSDCMCnb, 3, 3); // 대칭성 확보를 위한 trick
                    // Velocity
                    MatMult31(INSDCMCbn, Fb, fn);
                    double gtmp[3] = {0, 0, gravity_wgs84(INSRadPos[0], INSRadPos[2])};
                    VecSubstitute(gtmp, g, 3);
                    double delMpsVned[3] = {0};
                    double wnie2tmp[3] = {0};
                    double coriolitmp[3] = {0};
                    VecScalarMult(wnie, 2, wnie2tmp, 3);
                    MatSum(wnie2tmp, wnen, 3, 1, coriolitmp);
                    double coriolitmp2[3] = {0};
                    CrossProduct(coriolitmp, INSMpsVned, coriolitmp2);
                    MatSubtract(fn, coriolitmp2, 3, 1, delMpsVned);
                    MatSum(delMpsVned, g, 3, 1, delMpsVned);
                    double delvdt[3] = {0};
                    VecScalarMult(delMpsVned, dt, delvdt, 3);
                    MatSum(INSMpsVned, delvdt, 3, 1, INSMpsVned);
                    // Position
                    Radicurv(INSRadPos[0], &rm, &rp);
                    double delRadLat = INSMpsVned[0] / (rm + INSRadPos[2]) * dt;
                    // !! 함수 리턴 값이 degree로 처리되는 오류 발생할 수 있어 if else 구문 추가
                    double cosval = 0;
                    if ((1 - cos(INSRadPos[0])) < 1e-5)
                    {
                        cosval = cos(INSRadPos[0] * d2r);
                        // cout<<"cosine value is : "<< cosval << endl;
                    }
                    else
                    {
                        cosval = cos(INSRadPos[0]);
                        // cout<< setprecision(15) << "cosine value is : "<< cosval << endl;
                    }
                    double delRadLon = INSMpsVned[1] / ((rp + INSRadPos[2]) * cosval) * dt;
                    double delMetHei = -INSMpsVned[2] * dt;
                    double delpos[3] = {delRadLat, delRadLon, delMetHei};
                    MatSum(INSRadPos, delpos, 3, 1, INSRadPos);

                    // NHC _ 10.19.17:17
                    // MatMult31(INSDCMCnb, INSMpsVned, INSMpsVb);
                    // double measurement3[2] = {0, 0};
                    // double measest3[2] = {INSMpsVb[0], INSMpsVb[1]};
                    // MatSubtract((double *)measest3, (double *)measurement3, 2, 1, (double *)MeasurementNHC);
                    // InsKf15_updatelatlon(P_DB, estX, R_NHC, H_NHC, MeasurementNHC, P_DB, estX);
                    // Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
                    ///

                    // Filtering
                    GenINSFmatrix15(INSRadPos, INSMpsVned, INSDCMCbn, fn, wnin, wnie, wnen, F);
                    InsKf15_predict(F, estX, P_DB, Qvec, INSDCMCbn, INSDCMCnb, PHI, P_DB);
                    if (parsedIMUGPS.size() == 17 && gps_data[5] >= 6) // GPS정보가 같이 들어오고, 그 개수가 6개 이상일 때
                    {
                        // vector<double> PV_now(&parsedIMUGPS[12], &parsedIMUGPS[14]);
                        // vector<double> PV_previous(&parsedPrevious[12], &parsedPrevious[14]);
                        // if (((linenum - 1) / 3 + 1) - validrow <= t_finealign * sampleHZ) // finealign 실행하는 정지구간, 데이터의 유효한 데이터 들어오고 나서 t_finealign만큼
                        if (gps_data[0] <= firsttime + t_finealign) // finealign 실행하는 정지구간, 데이터의 유효한 데이터 들어오고 나서 t_finealign만큼
                        {
                            if ((linenum + 2) % 60 == 0)
                            {
                                double measurement2[5] = {gps_data[1], gps_data[2], 0, 0, 0};
                                double measest2[5] = {INSRadPos[0], INSRadPos[1], INSMpsVned[0], INSMpsVned[1], INSMpsVned[2]};
                                MatSubtract((double *)measest2, (double *)measurement2, 5, 1, (double *)MeasurementZupt);
                                InsKf15_updatezupt(P_DB, estX, R_zupt3, H_Zupt3, MeasurementZupt, P_DB, estX);
                                Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
                            }
                        }
                        else
                        {
                            if ((linenum + 2) % 60 != 0)
                            {
                            }
                            else // GPS의 위치, 속도만 사용 후 고도 정보는 barometer 사용
                            {
                                double measurement4[4] = {gps_data[1], gps_data[2], gpsvel_N, gpsvel_E};
                                double measest4[4] = {INSRadPos[0], INSRadPos[1], INSMpsVned[0], INSMpsVned[1]};
                                MatSubtract((double *)measest4, (double *)measurement4, 4, 1, (double *)Measurementlatlonvnve);
                                InsKf15_updatelatlonvnve(P_DB, estX, R_latlonvnve, H_latlonVnVe, Measurementlatlonvnve, P_DB, estX);
                                Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
                                // }
                                if ((sqrt(pow(Wbib[0], 2) + pow(Wbib[1], 2) + pow(Wbib[2], 2)) < 2 * d2r) && (gps_data[4] >= 5))
                                {
                                    double MeasurementYaw[1] = {INSRadEulr[2] - gps_data[6] * d2r};
                                    for (int yidx = 0; yidx < 3; yidx++)
                                    {
                                        MeasurementYaw[0] = MeasurementYaw[0] - (MeasurementYaw[0] > PI) * (2 * PI) - (MeasurementYaw[0] < -PI) * (-2 * PI);
                                    }
                                    InsKf15_updateHgt(P_DB, estX, R_yaw, H_yaw, MeasurementYaw, P_DB, estX);
                                    Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
                                }
                            }
                        }
                        parsedPrevious.assign(parsedIMUGPS.begin(), parsedIMUGPS.end());
                    }
                    // if (((linenum - 1) / 3 + 1) % int(BaroIMURatio) == 0 && (linenum - 1) / 3 - validrow > t_finealign * sampleHZ)
                    if (((linenum) + 1) % int(BaroIMURatio) == 0)
                    {
                        Baro_data += mean_BaroOffset;
                        double MeasurementHgt[1] = {INSRadPos[2] - Baro_data};
                        InsKf15_updateHgt(P_DB, estX, R_hgt, H_hgt, MeasurementHgt, P_DB, estX);
                        Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
                    } //barometer update
                    EstResult << setprecision(15) << fixed << scientific << INSRadPos[0] << " " << INSRadPos[1] << " " << INSRadPos[2] << " "
                              << INSMpsVned[0] << " " << INSMpsVned[1] << " " << INSMpsVned[2] << " " << INSRadEulr[0] << " "
                              << " " << INSRadEulr[1] << " " << INSRadEulr[2] << " " << Ba[0] << " "
                              << Ba[1] << " " << Ba[2] << " " << Bg[0] << " " << Bg[1] << " " << Bg[2] << endl;
                }
            }
            parsedIMUGPS.clear();
        }
    }
    EstResult.close();
    IMUGPSfile.close();
    // myGPS.close();
    // myhgt.close();
    // myvel.close();
    return 0;
}
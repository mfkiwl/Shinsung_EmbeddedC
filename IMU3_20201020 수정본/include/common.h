//#pragma once
#ifndef COMMON_H
#define COMMON_H

// Define Unit=========================================================
#define PI 3.141592653589793238463
#define d2r (double)(PI / 180.)
#define r2d (double)(180. / PI)
#define hr2s (double)(3600.0)
#define km2m (double)(1000.0)
#define m2km (double)(0.001)
#define dph2rps (double)(d2r / hr2s)
#define rps2dph (double)(1 / dph2rps)
#define ft2m (double)(0.3048)
#define earth_g (double)(9.781)
#define mg (double)(1e-3 * 9.781)
#define ug (double)(1e-6 * 9.781)
#define R0 (double)(6378137.0)
#define r0 (double)(6356752.3142)
#define f (double)(1 / 298.257223563)
#define e (double)(0.0818191908426)
#define We (double)(7.292115 * 1e-5)
#define light_speed (double)(299792458.0)

// IMU, aiding sensor(mag,baro,gps) specification ==================
// IMU
#define sampleHZ (double)(20.0)
#define g_repeatability (double)(1 * d2r)
#define g_ARW (double)(0.014 * d2r)
#define a_repeatability (double)(30 * mg)
#define a_ARW (double)(0.23 * mg)
#define dt (double)(1. / sampleHZ)
// Magnetometer
#define MagHZ (double)(20.0)
#define secDelt (double)(1. / smapleHZ)
#define MagIMURatio (double)(sampleHZ / MagHZ)
// Barometer
#define BaroHZ (double)(20.0)
#define noise_baro_h (double)(0.02)
#define BaroIMURatio (double)(sampleHZ / BaroHZ)
// GPS
// #define noise_hp (double)(0.01)
// #define noise_vp (double)(0.02)
// #define noise_v (double)(0.05)
#define noise_hp (double)(0.05) //10.20.10:38
#define noise_vp (double)(1)
#define noise_v (double)(0.1)

//======================================================================
#define t_align 120     // coarse align 시간 2분   //
#define t_finealign 300 // 정지구간 GPS 결합 (전원 인가 후 2분에서 5분),

#endif

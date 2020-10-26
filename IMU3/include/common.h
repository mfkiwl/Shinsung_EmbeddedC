//#pragma once
#ifndef COMMON_H
#define COMMON_H

// Define Unit=========================================================
#define PI 3.141592653589793238463
#define d2r (float64)(PI / 180.)
#define r2d (float64)(180. / PI)
#define hr2s (float64)(3600.0)
#define km2m (float64)(1000.0)
#define m2km (float64)(0.001)
#define dph2rps (float64)(d2r / hr2s)
#define rps2dph (float64)(1 / dph2rps)
#define ft2m (float64)(0.3048)
#define earth_g (float64)(9.781)
#define mg (float64)(1e-3 * 9.781)
#define ug (float64)(1e-6 * 9.781)
#define R0 (float64)(6378137.0)
#define r0 (float64)(6356752.3142)
#define f (float64)(1 / 298.257223563)
#define e (float64)(0.0818191908426)
#define We (float64)(7.292115 * 1e-5)
#define light_speed (float64)(299792458.0)

// IMU, aiding sensor(mag,baro,gps) specification ==================
// IMU
#define sampleHZ (float64)(20.0)
#define g_repeatability (float64)(1 * d2r)
#define g_ARW (float64)(0.014 * d2r)
#define a_repeatability (float64)(30 * mg)
#define a_ARW (float64)(0.23 * mg)
#define dt (float64)(1. / sampleHZ)
// Magnetometer
#define MagHZ (float64)(20.0)
#define secDelt (float64)(1. / smapleHZ)
#define MagIMURatio (float64)(sampleHZ / MagHZ)
// Barometer
#define BaroHZ (float64)(20.0)
#define noise_baro_h (float64)(0.02)
#define BaroIMURatio (float64)(sampleHZ / BaroHZ)
// GPS
// #define noise_hp (float64)(0.01)
// #define noise_vp (float64)(0.02)
// #define noise_v (float64)(0.05)
#define noise_hp (float64)(0.05) //10.20.10:38
#define noise_vp (float64)(1)
#define noise_v (float64)(0.1)

//======================================================================
#define t_align 120     // coarse align 시간 2분   //
#define t_finealign 300 // 정지구간 GPS 결합 (전원 인가 후 2분에서 5분),

#endif

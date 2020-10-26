#include "include/function.h"
#include "DSP28x_Project.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define SPI_RX_BUFF_LEN 34

// GLOBAL VAL
float64 R_latlon[2][2]; //32bit
float64 R_PV[6][6];
float64 R_zupt[3][3];	   // 10.19.13:22
float64 R_zupt2[3][3];	   //10.20.10:24
float64 R_zupt3[5][5];
float64 R_yaw[1][1];
float64 R_latlonvnve[4][4]; //
float64 P_DB[15][15];
float64 Qvec[6];
float64 H_latlon[2][15];
float64 H_PV[6][15];
float64 H_latlonVnVe[4][15]; // 10.20.10:43
float64 H_Zupt[3][15];		// 10.19.13:22
float64 H_hgt[1][15];
float64 H_Zupt3[5][15];
float64 H_yaw[1][15];
float64 F[15][15];
float64 PHI[15][15];
float64 R_hgt[1][1];

// Initializing (position, velocity, attitude)
float64 INSRadPos[3], INSDegPos[3], INSRadEulr[3], INSDegEulr[3], INSQT[4], INSDCMCnb[3][3], INSDCMCbn[3][3], INSMpsVned[3], INSMpsVb[3] = {0};
float64 Ba_init[3], Bg_init[3], Ba[3], Bg[3], wnin[3], wnie[3], wnen[3], wbin[3], wbnb[3] = {0};

//추정 오차
float64 estX[15] = {0};
// fb, wbib에서 바이어스를 뺀 값, Fb에 Cbn곱한 값
float64 Fb[3], Wbib[3], fn[3] = {0};
// 위치,속도만 사용하는 경우, 위도,경도만 사용하는 경우, 고도만 사용하는 경우 EKF의 측정치 (delta 측정치)
float64 MeasurementPV[6], Measurementlatlon[2] = {0}; // 중력 벡터   // 10.20.10:40
float64 Measurementlatlonvnve[4] = {0};
float64 MeasurementV[3] = {0}; //10.19.13:28
float64 MeasurementZupt[5] = {0};
float64 g[3] = {0};
// meridian radius of curvature, transverse radius of curvature
float64 rm, rp = 0;
// GPS 첫 데이터가 들어온 시점부터 정렬을 시작하기위해 invalidrow를 정의, validrow는 정렬을 시작하는 첫번째 열을 뜻함.
float64 validrow = 0;
// GPS 데이터 기준으로 사용 가능한 데이터가 들어오는지 판단
bool dataisincorrect = true;
// 정렬 여부 판단
bool aligned = false;
// 정렬 시간 까지 유효한 GPS데이터 개수 세는데 사용
int GPS_count = 0;
// 한 epoch마다 count
int linenum = 0;
int flagGPS = 0;
// Data save용 변수, previous gps 는 zupt 실행을 위해 직전 gps data와 비교
float64 time_p, fb[3], wbib[3], gps_data[7], Baro_data;
// 정렬에 사용되는 변수, 기압계 데이터가 전체적으로 bias되어 있어 이를 고려함.
float64 mean_acc[3], mean_gyr[3] = {0}, mean_BaroOffset, roll_init, pitch_init, yaw_init = 0;
// mean_BaroOffset을 계산하기 위해 사용되는 변수
float64 gps_tmp = 0;
float64 baro_tmp = 0;

float32 parsedIMUGPS[SPI_RX_BUFF_LEN / 2];
float32 parsedPrevious[SPI_RX_BUFF_LEN / 2];
Uint16 prevGPS[2];
float64 gpsvel_N = 0; //
float64 gpsvel_E = 0; //

char uart_buf[384];

// Test 1,SCIA  DLB, 8-bit word, baud rate 0x000F, default, 1 STOP bit, no parity
void scia_echoback_init()
{
	// Note: Clocks were turned on to the SCIA peripheral
	// in the InitSysCtrl() function
	SciaRegs.SCICCR.all = 0x0007;  // 1 stop bit,  No loopback
								   // No parity,8 char bits,
								   // async mode, idle-line protocol
	SciaRegs.SCICTL1.all = 0x0003; // enable TX, RX, internal SCICLK,
								   // Disable RX ERR, SLEEP, TXWAKE
	SciaRegs.SCICTL2.all = 0x0003;
	SciaRegs.SCICTL2.bit.TXINTENA = 1;
	SciaRegs.SCICTL2.bit.RXBKINTENA = 1;

#if 0
	SciaRegs.SCIHBAUD = 0x0003;		// 9600 baud @LSPCLK = 75.0MHz.
	SciaRegs.SCILBAUD = 0x00D0;
#else
	SciaRegs.SCIHBAUD = 0x0000; // 115200 baud @LSPCLK = 75.0MHz.
	//SciaRegs.SCILBAUD = 0x0050;		// 115200 baud @LSPCLK = 75.0MHz.
	SciaRegs.SCILBAUD = 0x0028; // 230400 baud
#endif
	SciaRegs.SCICTL1.all = 0x0023; // Relinquish SCI from Reset
}

// Initalize the UART FIFO
void scia_fifo_init()
{
	SciaRegs.SCIFFTX.all = 0xE040;
	SciaRegs.SCIFFRX.all = 0x204F;
	SciaRegs.SCIFFCT.all = 0x0;
}

// Transmit a character from the SCI
void scia_xmit(int a)
{
	while (SciaRegs.SCIFFTX.bit.TXFFST != 0)
	{
	}
	SciaRegs.SCITXBUF = a;
}

void scia_msg(char *msg)
{
	int i = 0;

	while (msg[i] != '\0')
	{
		scia_xmit(msg[i]);
		i++;
	}
}

// Initialize SPI FIFO registers
void spi_fifo_init()
{
	SpidRegs.SPIFFTX.all = 0xE040;
	SpidRegs.SPIFFRX.all = 0x204F;
	SpidRegs.SPIFFCT.all = 0x0;
}

void spi_init()
{
	//SpiaRegs.SPICCR.all =0x000F;
	SpidRegs.SPICCR.bit.SPISWRESET = 0;
	SpidRegs.SPICCR.bit.CLKPOLARITY = 0; // rising edge,  CLK_ Polarity  = 0
	SpidRegs.SPICCR.bit.SPICHAR = 0xFF;	 // 16-bit char bits
	SpidRegs.SPICCR.bit.SPILBK = 0;

	//SpidRegs.SPICTL.all =0x0006;
	SpidRegs.SPICTL.bit.MASTER_SLAVE = 1; // Master mode
	SpidRegs.SPICTL.bit.TALK = 1;		  //TX. Enable
	SpidRegs.SPICTL.bit.CLK_PHASE = 0;	  //CLK_Phase  = 0
	SpidRegs.SPICTL.bit.SPIINTENA = 1;

	SpidRegs.SPIBRR = 4;				// 18 MHz Clock setting
	SpidRegs.SPIPRI.bit.FREE = 1;		// Set so breakpoints don't disturb xmission
	SpidRegs.SPICCR.bit.SPISWRESET = 1; // reset
}

// "999"값을 찾아 spi data align를 맞춤
void process_trash()
{
	Uint16 t_buff[2] = {0, 0};
	float32 value;

	while (1)
	{
		t_buff[0] = t_buff[1];
		SpidRegs.SPITXBUF = 0x0;
		while (SpidRegs.SPIFFTX.bit.TXFFST != 0)
		{
		}

		t_buff[1] = SpidRegs.SPIRXBUF;
		DELAY_US(500);

		memcpy(&value, t_buff, 2);

		if (floor(value) == 999)
		{
			printf("find align data\n");
			break;
		}
	}
}

// read spi sensor data
void read_data()
{
	int i;
	Uint16 r_buf[SPI_RX_BUFF_LEN];

	memset(r_buf, 0x0, SPI_RX_BUFF_LEN);

	for (i = 0; i < SPI_RX_BUFF_LEN; i++)
	{
		SpidRegs.SPITXBUF = 0x00;
		while (SpidRegs.SPIFFTX.bit.TXFFST != 0)
		{
		}

		//while(SpidRegs.SPIFFRX.bit.RXFFST !=1) {} 	// abnormal
		r_buf[i] = SpidRegs.SPIRXBUF;
		DELAY_US(500);
	}

	memcpy(&parsedIMUGPS[1], r_buf, SPI_RX_BUFF_LEN - 2);

	// Raw data Gyro, Accelo Y and Z -1
	parsedIMUGPS[2] = parsedIMUGPS[2]*-1;
	parsedIMUGPS[3] = parsedIMUGPS[3]*-1;

	parsedIMUGPS[5] = parsedIMUGPS[5]*-1;
	parsedIMUGPS[6] = parsedIMUGPS[6]*-1;

	parsedIMUGPS[0]  = 300;
	flagGPS = 1;

	sprintf(uart_buf, "R,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%014.10f,%014.10f,%e,%e,%e,%e\r\n",
			parsedIMUGPS[1], parsedIMUGPS[2], parsedIMUGPS[3],
			parsedIMUGPS[4], parsedIMUGPS[5], parsedIMUGPS[6],
			parsedIMUGPS[7], parsedIMUGPS[8], parsedIMUGPS[9],
			parsedIMUGPS[10], parsedIMUGPS[11], parsedIMUGPS[12],
			parsedIMUGPS[13], parsedIMUGPS[14], parsedIMUGPS[15], parsedIMUGPS[16]);

	scia_msg(uart_buf);

	if (!aligned)
		linenum++;
}

void init_data_00()
{
	memset(R_latlon, 0x0, sizeof(float64) * 2 * 2);
	memset(R_PV, 0x0, sizeof(float64) * 6 * 6);
	memset(R_latlonvnve, 0x0, sizeof(float64) * 4 * 4); //
	memset(P_DB, 0x0, sizeof(float64) * 15 * 15);
	memset(R_zupt, 0x0, sizeof(float64) * 3 * 3);  //10.19.13:24
	memset(R_zupt2, 0x0, sizeof(float64) * 3 * 3); //10.20.10:25
	memset(R_zupt3,0x0, sizeof(float64) * 5 * 5);
	memset(R_yaw,0x0, sizeof(float64) * 2); // 246 line과 비교
	memset(Qvec, 0x0, sizeof(float64) * 6);
	memset(H_latlon, 0x0, sizeof(float64) * 2 * 15);
	memset(H_latlonVnVe, 0x0, sizeof(float64) * 4 * 15);
	memset(H_PV, 0x0, sizeof(float64) * 6 * 15);
	memset(H_Zupt, 0x0, sizeof(float64) * 3 * 15); // 10.19.13:28
	memset(H_Zupt3,0x0, sizeof(float64) * 5 * 15);
	memset(H_yaw, 0x0, sizeof(float64) * 1 * 15);
	memset(H_hgt, 0x0, sizeof(float64) * 1 * 15);
	memset(F, 0x0, sizeof(float64) * 15 * 15);
	memset(PHI, 0x0, sizeof(float64) * 15 * 15);
	memset(R_hgt, 0x0, sizeof(float64) * 2); // 10.20.10:52 2가 아니라 1아닌지?(porting 과정에서 오류?)

	memset(INSRadPos, 0x0, sizeof(float64) * 3);
	memset(INSDegPos, 0x0, sizeof(float64) * 3);
	memset(INSRadEulr, 0x0, sizeof(float64) * 3);
	memset(INSDegEulr, 0x0, sizeof(float64) * 3);

	memset(INSQT, 0x0, sizeof(float64) * 4);
	memset(INSDCMCnb, 0x0, sizeof(float64) * 3 * 3);
	memset(INSDCMCbn, 0x0, sizeof(float64) * 3 * 3);
	memset(INSMpsVned, 0x0, sizeof(float64) * 3);
	memset(INSMpsVb, 0x0, sizeof(float64) * 3);

	memset(Ba_init, 0x0, sizeof(float64) * 3);
	memset(Bg_init, 0x0, sizeof(float64) * 3);

	memset(Ba, 0x0, sizeof(float64) * 3);
	memset(Bg, 0x0, sizeof(float64) * 3);

	memset(wnin, 0x0, sizeof(float64) * 3);
	memset(wnie, 0x0, sizeof(float64) * 3);

	memset(wnen, 0x0, sizeof(float64) * 3);
	memset(wbin, 0x0, sizeof(float64) * 3);
	memset(wbnb, 0x0, sizeof(float64) * 3);

	memset(estX, 0x0, sizeof(float64) * 15);

	memset(Fb, 0x0, sizeof(float64) * 3);
	memset(Wbib, 0x0, sizeof(float64) * 3);
	memset(Measurementlatlonvnve, 0x0, sizeof(float64) * 4); //
	memset(fn, 0x0, sizeof(float64) * 3);
	memset(MeasurementV, 0x0, sizeof(float64) * 3); // 10.19.13:30
	memset(MeasurementPV, 0x0, sizeof(float64) * 6);
	memset(Measurementlatlon, 0x0, sizeof(float64) * 2); // 10.20.10:42
	memset(MeasurementZupt, 0x0, sizeof(float64) * 5);

	memset(g, 0x0, sizeof(float64) * 3);

	rm = 0;
	rp = 0;
	validrow = 0;
	dataisincorrect = true;
	aligned = false;
	GPS_count = 0;
	linenum = 0;
	time_p = 0;

	Baro_data = 0;

	memset(fb, 0x0, sizeof(float64) * 3);
	memset(wbib, 0x0, sizeof(float64) * 3);
	memset(gps_data, 0x0, sizeof(float64) * 7);

	memset(mean_acc, 0x0, sizeof(float64) * 3);
	memset(mean_gyr, 0x0, sizeof(float64) * 3);

	mean_BaroOffset = 0;
	roll_init = 0;
	pitch_init = 0;
	yaw_init = 0;

	gps_tmp = 0;
	baro_tmp = 0;

	memset(parsedIMUGPS, 0x0, sizeof(float64) * (SPI_RX_BUFF_LEN / 2));
	memset(parsedPrevious, 0x0, sizeof(float64) * (SPI_RX_BUFF_LEN / 2));
	memset(prevGPS, 0x0, sizeof(Uint16) * 2);
}

void init_data_01()
{
	float64 R_vec[2] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000)};
	float64 R_vec2[6] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_vp, noise_v, noise_v, noise_v};
	float64 R_vec3[4] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_v, noise_v}; //
	float64 R_vec4[3] = {noise_v, noise_v, noise_v};											 //10.20.10:26
	float64 R_vec5[3] = {noise_v * 2, noise_v * 2, noise_v * 2};								 //10.19.13:30
	float64 R_vec6[5] = {km2rad(noise_hp / 1000), km2rad(noise_hp / 1000), noise_v, noise_v, noise_v};
	Velwisesquare(R_vec2, 6, R_vec2);
	Velwisesquare(R_vec, 2, R_vec);
	Velwisesquare(R_vec3, 4, R_vec3); //
	Velwisesquare(R_vec4, 3, R_vec4); //10.19.13:30
	Velwisesquare(R_vec5, 3, R_vec5); //10.20.10:26
	Velwisesquare(R_vec6, 5, R_vec6);
	diag(R_vec, 2, (float64 *)R_latlon);
	diag(R_vec2, 6, (float64 *)R_PV);
	diag(R_vec3, 4, (float64 *)R_latlonvnve);
	diag(R_vec4, 3, (float64 *)R_zupt);	
	diag(R_vec5, 3, (float64 *)R_zupt2); 
	diag(R_vec6, 5, (float64 *)R_zupt3);
}

void init_data_02()
{
	float64 Ppos[3];
	float64 Pvel[3];
	float64 Patt[3];

	float64 pos_err_deviation[3] = {km2rad(5 * 1e-3), km2rad(5 * 1e-3), 10};
	float64 vel_err_deviation[3] = {noise_v, noise_v, noise_v};
	float64 att_err_deviation[3] = {a_repeatability / 9.81, a_repeatability / 9.81, 5 * d2r}; //07.27.19:26 debug complete.

	float64 PBa[3] = {pow(a_repeatability, 2), pow(a_repeatability, 2), pow(a_repeatability, 2)};
	float64 PBg[3] = {pow(g_repeatability, 2), pow(g_repeatability, 2), pow(g_repeatability, 2)};

	Velwisesquare(pos_err_deviation, 3, Ppos);
	Velwisesquare(vel_err_deviation, 3, Pvel);
	Velwisesquare(att_err_deviation, 3, Patt);

	float64 P_INS[15] = {Ppos[0], Ppos[1], Ppos[2], Pvel[0], Pvel[1], Pvel[2], Patt[0], Patt[1], Patt[2], PBa[0], PBa[1], PBa[2], PBg[0], PBg[1], PBg[2]};

	diag(P_INS, 15, (float64 *)P_DB);
}

void init_data_03()
{
	float64 Qwa[3] = {pow(a_ARW, 2), pow(a_ARW, 2), pow(a_ARW, 2)};
	float64 Qwg[3] = {pow(g_ARW, 2), pow(g_ARW, 2), pow(g_ARW, 2)};

	Qvec[0] = Qwa[0];
	Qvec[1] = Qwa[1];
	Qvec[2] = Qwa[2];
	Qvec[3] = Qwg[0];
	Qvec[4] = Qwg[1];
	Qvec[5] = Qwg[2];
}

void init_data_04()
{
	R_hgt[0][0] = pow(noise_baro_h, 2);
	R_yaw[0][0] = pow((0.5*d2r),2);

	H_latlon[0][0] = 1;
	H_latlon[1][1] = 1;
	H_hgt[0][2] = 1;
	H_PV[0][0] = 1;
	H_PV[1][1] = 1;
	H_PV[2][2] = 1;
	H_PV[3][3] = 1;
	H_PV[4][4] = 1;
	H_PV[5][5] = 1;

	H_Zupt[0][3] = 1; ///
	H_Zupt[1][4] = 1; ///10.19.13:32
	H_Zupt[2][5] = 1; ///

	H_latlonVnVe[0][0] = 1; //////////////////////////////////////////****
	H_latlonVnVe[1][1] = 1; //
	H_latlonVnVe[2][3] = 1; //
	H_latlonVnVe[3][4] = 1; //

	H_Zupt3[0][0] = 1;
    H_Zupt3[1][1] = 1;
    H_Zupt3[2][3] = 1;
    H_Zupt3[3][4] = 1;
    H_Zupt3[4][5] = 1;

	H_yaw[0][8] = -1;
}

// align 전 데이터 처리
int pre_process()
{
	int i;
	int align_count;

	float64 C1T[3][3];
	float64 C1T_T[3][3];
	float64 omega_tmp[3];

	memset(C1T, 0x0, sizeof(float64) * 3 * 3);
	memset(C1T_T, 0x0, sizeof(float64) * 3 * 3);
	memset(omega_tmp, 0x0, sizeof(float64) * 3);

	for (i = 0; i < 3; i++)
	{
		mean_acc[i] += fb[i];
		mean_gyr[i] += wbib[i];
	}

	if (flagGPS == 1)
	{
		gps_tmp += gps_data[3];
		baro_tmp += Baro_data;
		GPS_count++;
	}
	else
		baro_tmp += Baro_data;

	// 유효한 데이터 받고나서 t_align만큼 정렬
	if (linenum - validrow == sampleHZ * t_align)
	//if (gps_data[0] >= t_align) // pre_process를 gps 시간을 기준으로 판단
	{
		// ex)300초 데이터(25*300 = 7500개) 에 쓰레기 데이터가 1~50열 까지 있으면 정렬에 사용된 데이터는 51~7501 열(7501열에 GPS들어오는 첫 순간까지 정렬로 사용).
		align_count = linenum - (validrow - 1);
		for (i = 0; i < 3; i++)
		{
			mean_acc[i] /= align_count;
			mean_gyr[i] /= align_count;
			Bg_init[i] = mean_gyr[i]; // 초기 바이어스를 정렬값으로 설정
			mean_BaroOffset = gps_tmp / GPS_count - baro_tmp / align_count;
		}

		printf("%d number of IMU data is used to align \n", align_count);
		printf("%d number of GPS data is used to calculate altitude offset between GPS and Barometer. \n", GPS_count);

		roll_init = atan2(-mean_acc[1], -mean_acc[2]);
		pitch_init = atan2(mean_acc[0], sqrt(pow(mean_acc[1], 2) + pow(mean_acc[2], 2)));
		eulr2dcm(roll_init, pitch_init, 0, C1T);
		transpose((float64 *)C1T, (float64 *)C1T_T, 3, 3);
		MatMult31(C1T_T, mean_gyr, omega_tmp);
		yaw_init = atan2(-omega_tmp[1], omega_tmp[0]);

		// 추정치 초기화
		INSRadEulr[0] = roll_init;
		INSRadEulr[1] = pitch_init;
		INSRadEulr[2] = yaw_init; //기존코드: yaw init을 사용하지 않고 나침판어플(?)로 측정한 heading 사용함

		// INSRadEulr[2] = 8 * d2r;
		for (i = 0; i < 3; i++)
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
		eulr2qua(roll_init, pitch_init, yaw_init, INSQT);
		eulr2dcm(roll_init, pitch_init, yaw_init, INSDCMCbn);
		transpose((float64 *)INSDCMCbn, (float64 *)INSDCMCnb, 3, 3); // 대칭성 확보를 위한 trick
		MatMult31(INSDCMCnb, INSMpsVned, INSMpsVb);

#if 0
    	printf("%d %e %e %e %e %e %e %e %e %e \n\r", linenum, INSRadPos[0], INSRadPos[1], INSRadPos[2],
    			         	 	 	 	 	   INSMpsVned[0], INSMpsVned[1], INSMpsVned[2],
    										   INSRadEulr[0], INSRadEulr[1], INSRadEulr[2]);
#endif
		aligned = true;
	}

	return 1;
}

// align 데이터 처리
int align_process()
{
	// Inertial Navigation
	// Attitude
	float64 gtmp[3];
	float64 delMpsVned[3];
	float64 wnie2tmp[3];
	float64 coriolitmp[3];
	float64 coriolitmp2[3];
	float64 delvdt[3];
	float64 delRadLat;
	float64 cosval = 0;
	float64 delRadLon;
	float64 delMetHei;
	float64 delpos[3];

	// never referenced
	//float64 measurement1[6];
	//float64 measest1[6];
	float64 measurement2[3]; //10.20.10:47
	float64 measest2[3];		//10.20.10:47
	float64 measest4[4];		//////////////////////////////////////////////**
	float64 measurement3[2]; //
	float64 measurement4[4];
	float64 measest3[2];
	float64 measurement5[5];
	float64 measest5[5];

	float64 MeasurementHgt[1];
	float64 MeasurementYaw[1];

	MatSubtract(fb, Ba, 3, 1, Fb);
	MatSubtract(wbib, Bg, 3, 1, Wbib);
	NavRotation(INSMpsVned, INSRadPos, wnin, wnie, wnen);
	MatMult31(INSDCMCnb, wnin, wbin);
	MatSubtract(Wbib, wbin, 3, 1, wbnb);

	// 10.20.10:31
	//if (gps_data[4] <= 0.2 || (linenum - validrow) <= t_finealign * sampleHZ)
	if (linenum - validrow <= sampleHZ * t_finealign){} // finealign 동안 자세 업데이트 하지 않음
		QuatUpdate(wbnb, 0, INSQT);
	else
		QuatUpdate(wbnb, dt, INSQT);
	//
	// QuatUpdate(wbnb, dt, INSQT);


	qua2dcm(INSQT, INSDCMCbn);
	dcm2eulr(INSDCMCbn, INSRadEulr);
	transpose((float64 *)INSDCMCbn, (float64 *)INSDCMCnb, 3, 3); // 대칭성 확보를 위한 trick

	// Velocity
	MatMult31(INSDCMCbn, Fb, fn);

	gtmp[0] = 0;
	gtmp[1] = 0;
	gtmp[2] = gravity_wgs84(INSRadPos[0], INSRadPos[2]);

	VecSubstitute(gtmp, g, 3);

	memset(delMpsVned, 0x0, sizeof(float64) * 3);
	memset(wnie2tmp, 0x0, sizeof(float64) * 3);
	memset(coriolitmp, 0x0, sizeof(float64) * 3);

	VecScalarMult(wnie, 2, wnie2tmp, 3);
	MatSum(wnie2tmp, wnen, 3, 1, coriolitmp);

	memset(coriolitmp2, 0x0, sizeof(float64) * 3);
	CrossProduct(coriolitmp, INSMpsVned, coriolitmp2);
	MatSubtract(fn, coriolitmp2, 3, 1, delMpsVned);
	MatSum(delMpsVned, g, 3, 1, delMpsVned);

	memset(delvdt, 0x0, sizeof(float64) * 3);
	VecScalarMult(delMpsVned, dt, delvdt, 3);
	MatSum(INSMpsVned, delvdt, 3, 1, INSMpsVned);

	// Position
	Radicurv(INSRadPos[0], &rm, &rp);
	delRadLat = INSMpsVned[0] / (rm + INSRadPos[2]) * dt;

	// !! 함수 리턴 값이 degree로 처리되는 오류 발생할 수 있어 if else 구문 추가
	if ((1 - cos(INSRadPos[0])) < 1e-5)
	{
		cosval = cos(INSRadPos[0] * d2r);
	}
	else
	{
		cosval = cos(INSRadPos[0]);
	}

	delRadLon = INSMpsVned[1] / ((rp + INSRadPos[2]) * cosval) * dt;
	delMetHei = -INSMpsVned[2] * dt;
	delpos[0] = delRadLat;
	delpos[1] = delRadLon;
	delpos[2] = delMetHei;

	MatSum(INSRadPos, delpos, 3, 1, INSRadPos);

	// Filtering
	GenINSFmatrix15(INSRadPos, INSMpsVned, INSDCMCbn, fn, wnin, wnie, wnen, F);
	InsKf15_predict(F, estX, P_DB, Qvec, INSDCMCbn, INSDCMCnb, PHI, P_DB);
	if (flagGPS == 1 && gps_data[5] >= 6) // GPS정보가 같이 들어오고, 그 개수가 6개 이상일 때
	{
		// 10.20.10:28
		if (linenum - validrow <= t_finealign * sampleHZ) // finealign 실행하는 정지구간, 데이터의 유효한 데이터 들어오고 나서 t_finealign만큼
		//if(gps_data[0]<=t_finealign) // gps시간 기준으로 t_finealign 판단
		{
			if (linenum % 20 == 0) // 20샘플당 한 번만
			{
				measurement5[0] = gps_data[1];
				measurement5[1] = gps_data[2];
				measurement5[2] = 0;
				measurement5[3] = 0;
				measurement5[4] = 0;
				measest5[0] = INSRadPos[0];
				measest5[1] = INSRadPos[1];
				measest5[2] = INSMpsVned[0];
				measest5[1] = INSMpsVned[1];
				measest5[2] = INSMpsVned[2];
				MatSubtract((float64 *)measest5, (float64 *)measurement5, 5, 1, (float64 *)MeasurementZupt);
				InsKf15_updateV(P_DB, estX, R_zupt3, H_Zupt3, MeasurementZupt, P_DB, estX);
				Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
			}
		}
		else
		{
			gpsvel_N = gps_data[4] * cos(gps_data[6] * d2r); // deg -> rad
			gpsvel_E = gps_data[4] * sin(gps_data[6] * d2r);
			if (linenum % 20 == 0){   // 20샘플당 한번만
				//10.19.13:36   line#583 ~ 604(592line "else if")
				measurement4[0] = gps_data[1];
				measurement4[1] = gps_data[2];
				measurement4[2] = gpsvel_N;
				measurement4[3] = gpsvel_E;
				measest4[0] = INSRadPos[0];
				measest4[1] = INSRadPos[1];
				measest4[2] = INSMpsVned[0];
				measest4[3] = INSMpsVned[1];
				MatSubtract((float64 *)measest4, (float64 *)measurement4, 4, 1, (float64 *)Measurementlatlonvnve);
				InsKf15_updatelatlonvnve(P_DB, estX, R_latlonvnve, H_latlonVnVe, Measurementlatlonvnve, P_DB, estX);
				Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
				if ((sqrt(pow(Wbib[0], 2) + pow(Wbib[1], 2) + pow(Wbib[2], 2)) < 2 * d2r) && (gps_data[4] >= 5))
				{
					MeasurementYaw = INSRadEulr[2] - gps_data[6]*d2r;
					for (int yidx = 0; yidx < 3; yidx++)
					{
						MeasurementYaw[0] = MeasurementYaw[0] - (MeasurementYaw[0] > PI) * (2 * PI) - (MeasurementYaw[0] < -PI) * (-2 * PI);
					}
					InsKf15_updateHgt(P_DB, estX, R_yaw, H_yaw, MeasurementYaw, P_DB, estX);
					Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
				}
			}
			else // GPS
			{
			}
		}

		memcpy(&parsedPrevious[0], &parsedIMUGPS[0], sizeof(float64) * 16);
	}

	if ((linenum % (int)BaroIMURatio) == 0 )
	{
		Baro_data += mean_BaroOffset;
		MeasurementHgt[0] = INSRadPos[2] - Baro_data;
		InsKf15_updateHgt(P_DB, estX, R_hgt, H_hgt, MeasurementHgt, P_DB, estX);
		Correction(estX, INSRadPos, INSMpsVned, INSQT, INSDCMCbn, Ba, Bg);
	} //barometer update

	// processed data
	sprintf(uart_buf, "P,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le,%.10Le\r\n",
			INSRadPos[0], INSRadPos[1], INSRadPos[2],
			INSMpsVned[0], INSMpsVned[1], INSMpsVned[2],
			INSRadEulr[0], INSRadEulr[1], INSRadEulr[2],
			Ba[0], Ba[1], Ba[2], Bg[0], Bg[1], Bg[2]);
	scia_msg(uart_buf);

	return 1;
}

void main()
{
	float32 aa = 1.12345678912345L;
	float64 bb = 1.12345678912345L;

	init_data_00();
	init_data_01();
	init_data_02();
	init_data_03();
	init_data_04();

	InitSysCtrl();
	InitSpidGpio();

	// init spi
	spi_init();
	spi_fifo_init();

	// init UART
	InitSciaGpio();
	scia_fifo_init();
	scia_echoback_init();

	// uart buff clear
	memset(uart_buf, 0x0, 256);

	sprintf(uart_buf, "\r\n%f %.15Le\r\n", aa, bb);
	scia_msg(uart_buf);

	process_trash();

	while (1)
	{
		// receive spi data
		read_data();

		if (dataisincorrect)
		{
			// check gps condition
			if (floor(parsedIMUGPS[11]) == 37)
			{
				dataisincorrect = false;
				printf("The first non-trash value is at %d \n", linenum);
				validrow = linenum;
			}
		}

		if (!dataisincorrect)
		{
			// 시간
			time_p = parsedIMUGPS[0];
			// 자이로스코프
			wbib[0] = parsedIMUGPS[1];
			wbib[1] = parsedIMUGPS[2]; // 10.22 y,z축 센서값 부호 역전 -> 데이터 취득부분으로 옮김(200line)
			wbib[2] = parsedIMUGPS[3];
			// 가속도계 입력
			fb[0] = parsedIMUGPS[4];
			fb[1] = parsedIMUGPS[5];
			fb[2] = parsedIMUGPS[6];
			// 기압계 입력
			Baro_data = parsedIMUGPS[10];

			if (flagGPS == 1)
			{
				gps_data[0] = parsedIMUGPS[0];					// GPS 시간
				gps_data[1] = parsedIMUGPS[11] * d2r;			// 위도
				gps_data[2] = parsedIMUGPS[12] * d2r;			// 경도
				gps_data[3] = parsedIMUGPS[13];					// 고도
				gps_data[4] = parsedIMUGPS[14] * (1.852 / 3.6); // speed   //10.20:10:34
				gps_data[5] = parsedIMUGPS[15];					// 가시위성 개수
				gps_data[6] = parsedIMUGPS[16];					//
			}

			if (aligned)
				align_process();
			else
				pre_process();
		}

		memset(parsedIMUGPS, 0x0, sizeof(float64) * (SPI_RX_BUFF_LEN / 2));
	}
}

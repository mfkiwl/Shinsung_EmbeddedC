#pragma once
#ifndef READIMU_H
#define READIMU_H
#include "common.h"
#include <vector>
void CheckInputFile(); // IMU,GPS가 같이 load되는 파일에서 GPS정보 분리 기준 설정
void ReadData(string &strin, vector<double> &dataout);

#endif
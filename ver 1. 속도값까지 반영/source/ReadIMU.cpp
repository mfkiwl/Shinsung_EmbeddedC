#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "ReadIMU.h"
#include <vector>
string delimiter1 = " ";
string delimiter2 = ",";
string delimiter3 = "\r";
char *pEnd;
using namespace std;

void ReadData(string &strin, vector<double> &dataout)
{
    double data;
    double time;
    size_t pos = 0;
    string token;
    int skip = 0; // 알파벳 인덱스 버릴려고
    while (((pos = strin.find(delimiter1)) != string::npos) || (pos = strin.find(delimiter2)) != string::npos || (pos = strin.find(delimiter3)))
    {
        if (pos == strin.find(delimiter1))
        { // 시간 정보
            token = strin.substr(0, pos);
            char *str = new char[token.size() + 1];
            copy(token.begin(), token.end(), str);
            time = stod(token.substr(0, 2)) * 3600 + stod(token.substr(3, 2)) * 60 + stod(token.substr(6, 2)) + stod(token.substr(9, 3)) * 0.001;
            strin.erase(0, pos + delimiter1.length());
            dataout.push_back(time);
            delete[] str;
        }
        else if (pos == strin.find(delimiter2))
        {
            if (skip == 0)
            {
                token = strin.substr(0, pos);
                strin.erase(0, pos + delimiter2.length());
                skip++;
                continue;
            }
            else
            {
                token = strin.substr(0, pos);
                istringstream os(token);
                os >> data;
                strin.erase(0, pos + delimiter1.length());
                dataout.push_back(data);
            }
        }
        else if (pos == strin.find(delimiter3))
        {
            token = strin.substr(0, pos);
            istringstream os(token);
            os >> data;
            strin.erase(0, pos + delimiter1.length());
            dataout.push_back(data);
            break;
        }
    }
}

void CheckInputFile()
{
    ifstream IMUGPS;
    int linenum = 0;
    bool standardlength = true;
    IMUGPS.open("20191023_IMU_test.txt");
    if (IMUGPS.is_open())
    {
        cout << "File is opened successfully\n";
        while (!IMUGPS.eof())
        {
            linenum++;
            string str;
            getline(IMUGPS, str);
            if (str.length() > 90 && str.length() < 120)
            {
                cout << "In line num: " << linenum << ", its length is " << str.length() << endl;
                standardlength = false;
            }
            if (!standardlength)
            {
                cout << "choose another stadard length!" << endl;
                break;
            }
        }
        IMUGPS.close();
        if (standardlength)
        {
            cout << "will use 100 as a standard length for separating long info from short info\n";
        }
    }
}

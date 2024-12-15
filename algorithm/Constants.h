//
// Created by 王勇 on 2018/7/6.
//

#ifndef TDGT_CONSTANTS_H
#define TDGT_CONSTANTS_H

#include <string>
#include <map>
#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <cassert>


#define DE_INTV -2147483648
#define INTV_CNTED -2147483647
#define DE_W 2147483647


double EPSILON = 1e-5; //1/86400 =1.157e-5

int TMAX = 86400;
unsigned long LEAF_SIZE = 256;
unsigned long FANOUT = 4;

std::string tgraph_path = "../COL_3.txt";
std::string index_path = "../COL.idx";
std::string tdsp_query_path = "../Data/FLATDSP.demands";
std::string tisp_query_path = "../Data/FLATISP.demands";
int MAXShortcutNUM = 20;
std::string result_path = "../result.txt";

#endif //TDGT_CONSTANTS_H

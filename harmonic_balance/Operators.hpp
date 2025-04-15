#pragma once

#include <vector>
// see Liu et al 2006, square matrices, size  2*order+1

std::vector<double> getA(int order);
std::vector<double> getE(int order);
std::vector<double> getAsquared(int order);
std::vector<double> getInvE(int order);
std::vector<double> getD(int order);
std::vector<double> getDuffingR(int order);
std::vector<double> getDsquared(int order);
std::vector<double> getIdentity(int order);

/** Copyright 2016 Vitaly Vorobyev
 **
 **/

#include "../include/dlzbinstruct.h"

#include <fstream>
#include <iostream>
#include <cmath>

typedef std::vector<double> vectd;
typedef std::string str;

using std::ifstream;
using std::cout;
using std::endl;

namespace DDlz {

const str DlzBinStruct::m_format =
        "%d: C = %lf, S = %lf, K+ = %lf, K- = %lf";
const str DlzBinStruct::m_format_old =
        "%d: C = %lf, S = %lf, K = %lf, Kb = %lf, Q = %lf";

DlzBinStruct::DlzBinStruct(const unsigned nbins) :
    m_nbins(nbins) {
    m_C.resize(m_nbins, 0);
    m_S.resize(m_nbins, 0);
    m_Kp.resize(m_nbins, 0);
    m_Kn.resize(m_nbins, 0);
}

DlzBinStruct::DlzBinStruct(const str& fname) {
    m_nbins = GetParams(fname);
}

void DlzBinStruct::Clear(void) {
    m_C.clear();
    m_S.clear();
    m_Kp.clear();
    m_Kn.clear();
    m_D.clear();
    m_nbins = 0;
}

int DlzBinStruct::GetParams(const str& fname) {
    Clear();
    cout << "Getting parameters from " << fname << endl;
    ifstream ifile(fname, ifstream::in);
    if (!ifile.is_open()) return -1;
    str line;
    int bin;
    double C, S, Kp, Kn;
    int flag = 5;
    while (flag == 5) {
        int i = m_S.size();
        getline(ifile, line);
        flag = sscanf(line.c_str(), m_format.c_str(), &bin, &C, &S, &Kp, &Kn);
        if (i != (bin - 1) || flag != 5) break;
        m_nbins++;
        m_C.push_back(C);
        m_S.push_back(S);
        m_Kp.push_back(Kp);
        m_Kn.push_back(Kn);
        m_D.push_back(dilution(i + 1));
        cout << i+1 << " C: " << m_C[i] << " S: " << m_S[i] << " Kp: "
             << m_Kp[i] << " Kn: " << m_Kn[i] << " D: " << m_D[i] << endl;
    }
    cout << "Done! " << m_nbins << " bins" << endl;
    return m_nbins;
}

double DlzBinStruct::dilution(const int bin) const {
    if (!CheckBinNumber(bin)) return 0.;
    return 2. * m_C[bin-1] * sqrt(m_Kp[bin-1] * m_Kn[bin-1]) /
                                 (m_Kp[bin-1] + m_Kn[bin-1]);
}

bool DlzBinStruct::CheckBinNumber(const int bin) const {
    if (bin > 0 && bin <= m_nbins) {
        return true;
    } else {
        cout << "Wrong bin number detected: " << bin << " " << m_nbins << endl;
    }
    return false;
}

void DlzBinStruct::SetKp(const int bin, const double v) {
    if (!CheckBinNumber(bin)) return;
    m_Kp[bin-1] = v;
    m_D[bin-1] = dilution(bin);
}

void DlzBinStruct::SetKn(const int bin, const double v) {
    if (!CheckBinNumber(bin)) return;
    m_Kn[bin-1] = v;
    m_D[bin-1] = dilution(bin);
}

void DlzBinStruct::SetC(const int bin, const double v) {
    if (!CheckBinNumber(bin)) return;
    m_C[bin-1] = v;
    m_D[bin-1] = dilution(bin);
}

void DlzBinStruct::SetS(const int bin, const double v) {
    if (!CheckBinNumber(bin)) return;
    m_S[bin-1] = v;
}

void DlzBinStruct::SetD(const int bin, const double v) {
    if (!CheckBinNumber(bin)) return;
    m_D[bin-1] = v;
}

double DlzBinStruct::Kp(const int bin) const {
    if (!CheckBinNumber(abs(bin))) return 0;
    return bin > 0 ? m_Kp[bin-1] : m_Kn[-bin-1];
}

double DlzBinStruct::Kn(const int bin) const {
    if (!CheckBinNumber(abs(bin))) return 0;
    return bin > 0 ? m_Kn[bin - 1] : m_Kp[-bin - 1];
}

double DlzBinStruct::C(const int bin)  const {
    if (!CheckBinNumber(abs(bin))) return 0;
    return m_C[abs(bin) - 1];
}

double DlzBinStruct::S(const int bin)  const {
    if (!CheckBinNumber(abs(bin))) return 0;
    return bin > 0 ? m_S[bin - 1] : -m_S[-bin - 1];
}

double DlzBinStruct::D(const int bin) const {
    if (!CheckBinNumber(abs(bin))) return 0;
    return m_D[abs(bin) - 1];
}

}  // namespace DDlz

/** Copyright 2016 Vitaly Vorobyev
 **
 **/

#include "../include/ddlzbinstruct.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;
using std::fabs;

typedef std::string str;

namespace DDlz {

typedef DDlzBinStruct DD;

DD::DDlzBinStruct(const double& tau, const double& dm, const double& phi1,
                  const double& wtag, const str& fd, const str& fb) :
    m_sin2beta(sin(2. * phi1)), m_cos2beta(cos(2. * phi1)) {
    Init(tau, dm, wtag, fd, fb);
}

void DD::Init(const double& tau, const double& dm, const double& wtag,
                         const str& fd, const str& fb) {
    m_ddlz  = new DlzBinStruct(fd);
    m_bdlz  = new DlzBinStruct(fb);
    m_tau   = tau;
    m_dm    = dm;
    m_wrtag = wtag;
    CalcXiDelut();
}

void DD::CalcKKKK(const int binb, const int bind) {
    const double kpkp = m_ddlz->Kp(bind) * m_bdlz->Kp(binb);
    const double knkn = m_ddlz->Kn(bind) * m_bdlz->Kn(binb);
    m_kkpkk = kpkp + knkn;
    m_kkmkk = kpkp - knkn;
}

void DD::CalcXiDelut(void) {
    m_xi    = 1. / (1. + pow(m_tau * m_dm, 2));
    m_delut = 1. - 2. * m_wrtag;
}

void DD::Set_tau(const double& x)      { m_tau = x; CalcXiDelut();}
void DD::Set_dm(const double& x)       { m_dm = x; CalcXiDelut();}
void DD::Set_wrtag(const double& x)    { m_wrtag = x; CalcXiDelut();}
void DD::Set_sin2beta(const double& x) { m_sin2beta = x;}
void DD::Set_cos2beta(const double& x) { m_cos2beta = x;}
void DD::Set_phi1(const double& x)     { m_sin2beta = sin(2. * x);
                                         m_cos2beta = cos(2. * x); }
void DD::Set_C_bdlz(const int bin, const double& x) { m_bdlz->SetC(bin, x);}
void DD::Set_S_bdlz(const int bin, const double& x) { m_bdlz->SetS(bin, x);}
void DD::Set_D_bdlz(const int bin, const double& x) { m_bdlz->SetD(bin, x);}
void DD::Set_C_ddlz(const int bin, const double& x) { m_ddlz->SetC(bin, x);}
void DD::Set_S_ddlz(const int bin, const double& x) { m_ddlz->SetS(bin, x);}
void DD::Set_D_ddlz(const int bin, const double& x) { m_ddlz->SetD(bin, x);}

double DD::FractionD0pipi(const int bin, const int flv) const {
    return 0.5 * ( (m_bdlz->Kp(bin) + m_bdlz->Kn(bin)) +
                 m_delut * m_xi * flv * (m_bdlz->Kp(bin) - m_bdlz->Kn(bin)) );
}

double DD::FractionKspipi(const int bin, const int flv) const {
    return 0.5 * ( (m_ddlz->Kp(bin) + m_ddlz->Kn(bin)) +
                 m_delut * m_xi * flv * (m_ddlz->Kp(bin) - m_ddlz->Kn(bin)) );
}

double DD::FractionDblDlz(const int binb, const int bind, const int flv) {
    CalcKKKK(binb, bind);
    return 0.5 * (m_kkpkk + m_delut * flv * m_xi * m_kkmkk);
}

double DD::SinCoefFlv(void) const {  return 0;}
double DD::CosCoefFlv(const int flv) const { return m_delut*flv;}

// * Functions with full B Dalitz plot * //
double DD::SinCoefCP(const int flv, const int cp, const int bbin) const {
    const double flvfact = -2. * m_delut * flv * cp;
    const double kfact = sqrt(m_bdlz->Kp(bbin) * m_bdlz->Kn(bbin)) /
                             (m_bdlz->Kp(bbin) + m_bdlz->Kn(bbin));
    const double cpfact = m_bdlz->S(bbin) * m_cos2beta -
                          m_bdlz->C(bbin) * m_sin2beta;
    return flvfact * kfact * cpfact;
}

double DD::CosCoefCP(const int flv, const int bbin) const {
    const double coef = m_delut * flv * (m_bdlz->Kp(bbin) - m_bdlz->Kn(bbin)) /
                                        (m_bdlz->Kp(bbin) + m_bdlz->Kn(bbin));
    if (fabs(coef) > 1) {
        cerr << "|C| > 1: " << coef << " = "
             << m_delut << " * " << flv << " * ("
             << m_bdlz->Kp(bbin) << " - " << m_bdlz->Kn(bbin) << ") / ("
             << m_bdlz->Kp(bbin) << " + " << m_bdlz->Kn(bbin) << ")"
             << endl;
    }
    return coef;
}

double DD::SinCoefDD(const int flv, const int dbin, const int bbin) {
    CalcKKKK(bbin, dbin);
    const double flvfact = -2. * m_delut * flv;
    const double kfact = sqrt(m_bdlz->Kp(bbin) * m_bdlz->Kn(bbin) *
                              m_ddlz->Kp(dbin) * m_ddlz->Kn(dbin)) /
                              m_kkpkk;
    const double cspsc = m_ddlz->C(dbin) * m_bdlz->S(bbin) +
                         m_ddlz->S(dbin) * m_bdlz->C(bbin);
    const double ccmss = m_ddlz->C(dbin) * m_bdlz->C(bbin) -
                         m_ddlz->S(dbin) * m_bdlz->S(bbin);
    const double cpfact = cspsc * m_cos2beta - ccmss * m_sin2beta;
    return flvfact * kfact * cpfact;
}

double DD::CosCoefDD(const int flv, const int dbin, const int bbin) {
    CalcKKKK(bbin, dbin);
    return m_delut * flv * m_kkmkk / m_kkpkk;
}

// * Functions with dilution factors * //
double DD::DSinCoefCP(const int flv, const int cp, const int bbin) const {
    return flv * cp * m_delut * m_bdlz->D(bbin) * m_sin2beta;
}

double DD::DCosCoefCP(void) const { return 0;}

double DD::DSinCoefDD(const int flv, const int dbin, const int bbin) {
    const double flvfact = -2. * m_bdlz->D(bbin) * m_delut * flv;
    const double kfact   = sqrt(m_ddlz->Kp(dbin) * m_ddlz->Kn(dbin)) /
                               (m_ddlz->Kp(dbin) + m_ddlz->Kn(dbin));
    const double cpfact = m_ddlz->S(dbin) * m_cos2beta -
                          m_ddlz->C(dbin) * m_sin2beta;
    return flvfact * kfact * cpfact;
}

double DD::DCosCoefDD(const int flv, const int dbin) {
    const double coef = m_delut * flv * (m_ddlz->Kp(dbin) - m_ddlz->Kn(dbin)) /
                                        (m_ddlz->Kp(dbin) + m_ddlz->Kn(dbin));
    if (fabs(coef) > 1) {
        cerr << "|C| > 1: " << coef << " = "
             << m_delut << " * " << flv << " * ("
             << m_ddlz->Kp(dbin) << " - " << m_ddlz->Kn(dbin) << ") / ("
             << m_ddlz->Kp(dbin) << " + " << m_ddlz->Kn(dbin) << ")"
             << endl;
    }
    return coef;
}

void DD::print_params(void) const {
    cout << "DDPdf init info:" << endl
         << "  tau: " << m_tau << endl
         << "  dm: " << m_dm << endl
         << "  cos(beta): " << m_cos2beta
         << ", sin(beta): " << m_sin2beta << endl
         << "  wrong tag: " << m_wrtag << endl
         << "  D pars:" << endl;
    for (unsigned i = 0; i < m_ddlz->nbins(); i++) {
        cout << "    " << i+1 << ": "
             << m_ddlz->Kp(i+1) << " "
             << m_ddlz->Kn(i+1) << " "
             << m_ddlz->C(i+1) << " "
             << m_ddlz->S(i+1) << endl;
    }
    cout << "  B pars:" << endl;
    for (unsigned i = 0; i < m_bdlz->nbins(); i++) {
        cout << "    " << i+1 << ": "
             << m_bdlz->Kp(i+1) << " "
             << m_bdlz->Kn(i+1) << " "
             << m_bdlz->C(i+1) << " "
             << m_bdlz->S(i+1) << endl;
    }
}

}  // namespace DDlz

/** Copyright 2016 Vitaly Vorobyev
 **
 **/

#pragma once

#include <string>
#include <memory>

#include "./dlzbinstruct.h"

namespace DDlz {

///
/// \brief The DDlzBinStruct class
///
class DDlzBinStruct {
    mutable std::unique_ptr<DlzBinStruct> m_ddlz;
    mutable std::unique_ptr<DlzBinStruct> m_bdlz;

    double m_tau;
    double m_dm;
    double m_wrtag;

    double m_sin2beta;
    double m_cos2beta;

    double m_xi, m_delut;

    void CalcKKKK(int binb, int bind) const;
    void CalcXiDelut(void);

    mutable double m_kkpkk, m_kkmkk;

 public:
    ///
    /// \brief DDlzBinStruct
    /// \param tau
    /// \param dm
    /// \param phi1
    /// \param wtag
    /// \param fd
    /// \param fb
    ///
    DDlzBinStruct(double tau, double dm, double phi1, double wtag,
                  const std::string &fd, const std::string &fb);
    inline auto tau() const {return m_tau;}
    inline auto dm() const {return m_dm;}
  // Modificators //
    ///
    /// \brief Set_tau
    /// \param x
    ///
    void Set_tau(double x);
    ///
    /// \brief Set_dm
    /// \param x
    ///
    void Set_dm(double x);
    ///
    /// \brief Set_wrtag
    /// \param x
    ///
    void Set_wrtag(double x);
    ///
    /// \brief Set_sin2beta
    /// \param x
    ///
    void Set_sin2beta(double x);
    ///
    /// \brief Set_cos2beta
    /// \param x
    ///
    void Set_cos2beta(double x);
    ///
    /// \brief Set_phi1
    /// \param x
    ///
    void Set_phi1(double x);
    ///
    /// \brief Set_C_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_C_bdlz(int bin, double x);
    ///
    /// \brief Set_S_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_S_bdlz(int bin, double x);
    ///
    /// \brief Set_D_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_D_bdlz(int bin, double x);
    ///
    /// \brief Set_C_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_C_ddlz(int bin, double x);
    ///
    /// \brief Set_S_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_S_ddlz(int bin, double x);
    ///
    /// \brief Set_D_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_D_ddlz(int bin, double x);

  // Constant methods //
    ///
    /// \brief cos_beta
    /// \return
    ///
    double cos_beta(void) const {return m_cos2beta;}
    ///
    /// \brief sin_beta
    /// \return
    ///
    double sin_beta(void) const {return m_sin2beta;}
    ///
    /// \brief wt
    /// \return
    ///
    double wt(void) const {return m_wrtag;}
    ///
    /// \brief C_bdlz. Averaged cos of strong phase difference
    /// for B0 -> D0 pi+ pi-
    /// \param bin
    /// \return
    ///
    double C_bdlz(int bin) const {return m_bdlz->C(bin);}
    ///
    /// \brief S_bdlz. Averaged sin of strong phase difference
    /// for B0 -> D0 pi+ pi-
    /// \param bin
    /// \return
    ///
    double S_bdlz(int bin) const {return m_bdlz->S(bin);}
    ///
    /// \brief D_bdlz. Dilution factor for B0 -> D0 pi+ pi-
    /// \param bin
    /// \return
    ///
    double D_bdlz(int bin) const {return m_bdlz->D(bin);}
    ///
    /// \brief FractionKspipi
    /// \param bin
    /// \param flv
    /// \return
    ///
    double FractionKspipi(int bin, int flv) const;
    ///
    /// \brief FractionD0pipi
    /// \param bin
    /// \param flv
    /// \return
    ///
    double FractionD0pipi(int bin, int flv) const;
    ///
    /// \brief FractionDblDlz
    /// \param binb
    /// \param bind
    /// \param flv
    /// \return
    ///
    double FractionDblDlz(int binb, int bind, int flv);
    ///
    /// \brief DPars
    /// \return
    ///
    const DlzBinStruct& DPars(void) const { return *m_ddlz;}
    ///
    /// \brief BPars
    /// \return
    ///
    const DlzBinStruct& BPars(void) const { return *m_bdlz;}
    ///
    /// \brief SinCoefFlv
    /// \return
    ///
    double SinCoefFlv(void) const;
    ///
    /// \brief CosCoefFlv
    /// \param flv
    /// \return
    ///
    double CosCoefFlv(int flv) const;
    ///
    /// \brief SinCoefCP. sin(dt) factor for Dcp
    /// \param flv
    /// \param cp
    /// \param bbin
    /// \return
    ///
    double SinCoefCP(int flv, int cp, int bbin) const;
    ///
    /// \brief CosCoefCP. cos(dt) factor for Dcp
    /// \param flv
    /// \param bbin
    /// \return
    ///
    double CosCoefCP(int flv, int bbin) const;
    ///
    /// \brief SinCoefDD. sin(dt) factor for double Dalitz
    /// \param flv
    /// \param dbin
    /// \param bbin
    /// \return
    ///
    double SinCoefDD(int flv, int dbin, int bbin) const;
    ///
    /// \brief CosCoefDD. cos(dt) factor for double Dalitz
    /// \param flv
    /// \param dbin
    /// \param bbin
    /// \return
    ///
    double CosCoefDD(int flv, int dbin, int bbin) const;
    ///
    /// \brief DSinCoefCP. sin(dt) factor for Dcp and averaged B Dalitz plot
    /// \param flv
    /// \param cp
    /// \param bbin
    /// \return
    ///
    double DSinCoefCP(int flv, int cp, int bbin) const;
    ///
    /// \brief DCosCoefCP. cos(dt) factor for Dcp and averaged B Dalitz plot
    /// \return
    ///
    double DCosCoefCP(void) const;
    ///
    /// \brief DSinCoefDD. sin(dt) factor for double Dalitz and
    /// averaged B Dalitz plot
    /// \param flv
    /// \param dbin
    /// \param bbin
    /// \return
    ///
    double DSinCoefDD(int flv, int dbin, int bbin) const;
    ///
    /// \brief DCosCoefDD. cos(dt) factor for double Dalitz and averaged
    /// B Dalitz plot
    /// \param flv
    /// \param dbin
    /// \return
    ///
    double DCosCoefDD(int flv, int dbin) const;
    ///
    /// \brief print_params
    ///
    void print_params(void) const;
};

}  // namespace DDlz

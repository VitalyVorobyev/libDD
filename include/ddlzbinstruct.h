/** Copyright 2016 Vitaly Vorobyev
 **
 **/

#ifndef _HOME_VITALY_MYLIBS_LIBDD_INCLUDE_DDLZBINSTRUCT_H_
#define _HOME_VITALY_MYLIBS_LIBDD_INCLUDE_DDLZBINSTRUCT_H_

#include <string>

#include "./dlzbinstruct.h"

namespace DDlz {

///
/// \brief The DDlzBinStruct class
///
class DDlzBinStruct {
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
    DDlzBinStruct(const double &tau, const double &dm, const double& phi1,
                  const double &wtag, const std::string &fd,
                  const std::string &fb);
    ///
    /// \brief Init
    /// \param tau
    /// \param dm
    /// \param wtag
    /// \param fd
    /// \param fb
    ///
    void Init(const double &tau, const double &dm, const double &wtag,
              const std::string &fd, const std::string &fb);

  // Modificators //
    ///
    /// \brief Set_tau
    /// \param x
    ///
    void Set_tau(const double& x);
    ///
    /// \brief Set_dm
    /// \param x
    ///
    void Set_dm(const double& x);
    ///
    /// \brief Set_wrtag
    /// \param x
    ///
    void Set_wrtag(const double& x);
    ///
    /// \brief Set_sin2beta
    /// \param x
    ///
    void Set_sin2beta(const double& x);
    ///
    /// \brief Set_cos2beta
    /// \param x
    ///
    void Set_cos2beta(const double& x);
    ///
    /// \brief Set_phi1
    /// \param x
    ///
    void Set_phi1(const double& x);
    ///
    /// \brief Set_C_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_C_bdlz(const int bin, const double& x);
    ///
    /// \brief Set_S_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_S_bdlz(const int bin, const double& x);
    ///
    /// \brief Set_D_bdlz
    /// \param bin
    /// \param x
    ///
    void Set_D_bdlz(const int bin, const double& x);
    ///
    /// \brief Set_C_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_C_ddlz(const int bin, const double& x);
    ///
    /// \brief Set_S_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_S_ddlz(const int bin, const double& x);
    ///
    /// \brief Set_D_ddlz
    /// \param bin
    /// \param x
    ///
    void Set_D_ddlz(const int bin, const double& x);

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
    double C_bdlz(const int bin) const {return m_bdlz->C(bin);}
    ///
    /// \brief S_bdlz. Averaged sin of strong phase difference
    /// for B0 -> D0 pi+ pi-
    /// \param bin
    /// \return
    ///
    double S_bdlz(const int bin) const {return m_bdlz->S(bin);}
    ///
    /// \brief D_bdlz. Dilution factor for B0 -> D0 pi+ pi-
    /// \param bin
    /// \return
    ///
    double D_bdlz(const int bin) const {return m_bdlz->D(bin);}
    ///
    /// \brief FractionKspipi
    /// \param bin
    /// \param flv
    /// \return
    ///
    double FractionKspipi(const int bin, const int flv) const;
    ///
    /// \brief FractionD0pipi
    /// \param bin
    /// \param flv
    /// \return
    ///
    double FractionD0pipi(const int bin, const int flv) const;
    ///
    /// \brief FractionDblDlz
    /// \param binb
    /// \param bind
    /// \param flv
    /// \return
    ///
    double FractionDblDlz(const int binb, const int bind, const int flv);
    ///
    /// \brief DPars
    /// \return
    ///
    const DlzBinStruct* DPars(void) const { return m_ddlz;}
    ///
    /// \brief BPars
    /// \return
    ///
    const DlzBinStruct* BPars(void) const { return m_bdlz;}
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
    double CosCoefFlv(const int flv) const;
    ///
    /// \brief SinCoefCP. sin(dt) factor for Dcp
    /// \param flv
    /// \param cp
    /// \param bbin
    /// \return
    ///
    double SinCoefCP(const int flv, const int cp, const int bbin) const;
    ///
    /// \brief CosCoefCP. cos(dt) factor for Dcp
    /// \param flv
    /// \param bbin
    /// \return
    ///
    double CosCoefCP(const int flv, const int bbin) const;
    ///
    /// \brief SinCoefDD. sin(dt) factor for double Dalitz
    /// \param flv
    /// \param dbin
    /// \param bbin
    /// \return
    ///
    double SinCoefDD(const int flv, const int dbin, const int bbin);
    ///
    /// \brief CosCoefDD. cos(dt) factor for double Dalitz
    /// \param flv
    /// \param dbin
    /// \param bbin
    /// \return
    ///
    double CosCoefDD(const int flv, const int dbin, const int bbin);
    ///
    /// \brief DSinCoefCP. sin(dt) factor for Dcp and averaged B Dalitz plot
    /// \param flv
    /// \param cp
    /// \param bbin
    /// \return
    ///
    double DSinCoefCP(const int flv, const int cp, const int bbin) const;
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
    double DSinCoefDD(const int flv, const int dbin, const int bbin);
    ///
    /// \brief DCosCoefDD. cos(dt) factor for double Dalitz and averaged
    /// B Dalitz plot
    /// \param flv
    /// \param dbin
    /// \return
    ///
    double DCosCoefDD(const int flv, const int dbin);
    ///
    /// \brief print_params
    ///
    void print_params(void) const;

 protected:
    double m_tau;
    double m_dm;

 private:
    DlzBinStruct* m_ddlz;
    DlzBinStruct* m_bdlz;

    double m_wrtag;
    double m_sin2beta;
    double m_cos2beta;

    double m_xi, m_delut;

    void CalcKKKK(const int binb, const int bind);
    void CalcXiDelut(void);

    double m_kkpkk, m_kkmkk;
};

}  // namespace DDlz

#endif  // _HOME_VITALY_MYLIBS_LIBDD_INCLUDE_DDLZBINSTRUCT_H_

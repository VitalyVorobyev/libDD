/** Copyright 2016 Vitaly Vorobyev
 **
 **/

#pragma once

#include <vector>
#include <string>

namespace DDlz {

///
/// \brief The DlzBinStruct class. Class for keeping binned
/// Dalitz plot parameters
///
class DlzBinStruct {
    ///
    /// \brief dilution
    /// \param bin
    /// \return
    ///
    double dilution(int bin) const;
    ///
    /// \brief CheckBinNumber
    /// \param bin
    /// \return
    ///
    bool CheckBinNumber(int bin) const;
    ///
    /// \brief m_nbins
    ///
    int m_nbins;
    ///
    /// \brief m_C
    ///
    std::vector<double> m_C;
    ///
    /// \brief m_S
    ///
    std::vector<double> m_S;
    ///
    /// \brief m_Kp
    ///
    std::vector<double> m_Kp;
    ///
    /// \brief m_Kn
    ///
    std::vector<double> m_Kn;
    ///
    /// \brief m_D
    ///
    std::vector<double> m_D;
    /**
     * @brief m_format. Config file string format
     */
    static const std::string m_format;
    /**
     * @brief m_format_old. Deprecated
     */
    static const std::string m_format_old;

 public:
    ///
    /// \brief DlzBinStruct
    /// \param nbins
    ///
    explicit DlzBinStruct(unsigned nbins);
    ///
    /// \brief DlzBinStruct
    /// \param fname
    ///
    explicit DlzBinStruct(const std::string& fname);
    ///
    /// \brief SetKp. Set decay probability for positive bin
    /// \param bin
    /// \param v
    ///
    void SetKp(int bin, double v);
    ///
    /// \brief SetKn. Set decay probability for negative bin
    /// \param bin
    /// \param v
    ///
    void SetKn(int bin, double v);
    ///
    /// \brief SetC. Set average cos of strong phase difference
    /// \param bin
    /// \param v
    ///
    void SetC(int bin, double v);
    ///
    /// \brief SetS. Set average sin of strong phase difference
    /// \param bin
    /// \param v
    ///
    void SetS(int bin, double v);
    ///
    /// \brief SetD. Set dilution factor
    /// \param bin
    /// \param v
    ///
    void SetD(int bin, double v);
    ///
    /// \brief Kp. Decay probability for positive bins
    /// \param bin
    /// \return
    ///
    double Kp(int bin) const;
    ///
    /// \brief nbins
    /// \return
    ///
    unsigned nbins(void) const {return m_Kp.size();}
    ///
    /// \brief Kn. Decay probability for negative bins
    /// \param bin
    /// \return
    ///
    double Kn(int bin) const;
    ///
    /// \brief C. Average cos of strong phase difference
    /// \param bin
    /// \return
    ///
    double C(int bin)  const;
    ///
    /// \brief S. Average sin of strong phase difference
    /// \param bin
    /// \return
    ///
    double S(int bin)  const;
    ///
    /// \brief D. Dilution factor
    /// \param bin
    /// \return
    ///
    double D(int bin)  const;
    ///
    /// \brief GetParams
    /// \param fname
    /// \return
    ///
    int GetParams(const std::string& fname);
    ///
    /// \brief Clear
    ///
    void Clear(void);
};

}  // namespace DDlz

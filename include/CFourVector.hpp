#ifndef COMPLEX_FOUR_VECTOR_H
#define COMPLEX_FOUR_VECTOR_H

#include <complex>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <TLorentzVector.h>

using cdouble = std::complex<double>;

class CFourVector {
 public:
  CFourVector();
  CFourVector(const CFourVector&);
  explicit CFourVector(const cdouble*);
  explicit CFourVector(const double*);
  CFourVector(cdouble, cdouble, cdouble, cdouble);
  explicit CFourVector(const TLorentzVector&);
  CFourVector(const TLorentzVector&, double);
  cdouble getM2() const;
  double getCosTheta() const;
  cdouble getX() const;
  cdouble getY() const;
  cdouble getZ() const;
  cdouble getE() const;
  cdouble& operator[](std::size_t);
  const cdouble& operator[](std::size_t) const;
  friend CFourVector operator + (const CFourVector&, const CFourVector&);
  CFourVector operator - (const CFourVector&) const;
  friend cdouble operator * (const CFourVector&, const CFourVector&);
  friend CFourVector operator * (const CFourVector&, const cdouble&);
  friend CFourVector operator * (const cdouble&, const CFourVector&);
  friend CFourVector operator / (const CFourVector&, const cdouble&);
  CFourVector& operator += (const CFourVector&);
  CFourVector& operator -= (const CFourVector&);
  CFourVector& operator *= (const cdouble&);
  CFourVector operator - () const;
  friend std::ostream & operator << (std::ostream&, const CFourVector&);
  friend CFourVector epsilon_conv(
      const CFourVector&,
      const CFourVector&,
      const CFourVector&);
  friend double electron_current_conv(
      const CFourVector&,
      const CFourVector&,
      const CFourVector&);
  friend double phase_space_weight(
      const std::vector<CFourVector>&);

 private:
  static cdouble permutation_parity(
      const std::string::const_iterator&,
      const std::string::const_iterator&);
  static std::unordered_map<std::string, cdouble> getLeviCivita();
  static const std::unordered_map<std::string, cdouble> epsilon_;
  cdouble x_[4];
};

#endif

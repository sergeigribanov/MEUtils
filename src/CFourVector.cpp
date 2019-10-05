#include <algorithm>
#include <numeric>
#include <cmath>
#include "TMath.h"
#include "CFourVector.h"

cdouble CFourVector::permutation_parity(
      const std::string::const_iterator& itBegin,
      const std::string::const_iterator& itEnd) {
  int nInv = 0;
  for (auto it = itBegin; it != itEnd - 1; ++it) {
    for (auto jt = it + 1; jt != itEnd; ++jt) {
      if (*it > *jt) {
        nInv++;
      }
    }
  }
  if (nInv % 2 == 0) {
    return cdouble(1., 0);
  }
  return cdouble(-1., 0);
}

std::unordered_map<std::string, cdouble>
CFourVector::getLeviCivita() {
  std::string s = "0123";
  std::unordered_map<std::string, cdouble> res;
  do {
    res[s] = permutation_parity(s.begin(), s.end());
  } while (std::next_permutation(s.begin(), s.end()));
  return res;
}

const std::unordered_map<std::string, cdouble>
CFourVector::epsilon_ = getLeviCivita();

CFourVector epsilon_conv(
      const CFourVector& p1,
      const CFourVector& p2,
      const CFourVector& p3) {
  CFourVector res;
  std::string s = "0123";
  int i, j, k, l;
  do {
    i = s[0] - '0';
    j = s[1] - '0';
    k = s[2] - '0';
    l = s[3] - '0';
    res[i] += CFourVector::epsilon_.at(s) *
        p1[j] * p2[k] * p3[l];
  } while (std::next_permutation(s.begin(), s.end()));
  return res;
}

CFourVector::CFourVector() {
}

CFourVector::CFourVector(const CFourVector& v):
    x_(v.x_) {
}

CFourVector::CFourVector(const cdouble* array) {
  for (std::size_t i = 0; i < 4; ++i) {
    x_[i] = array[i];
  }
}

CFourVector::CFourVector(const double* array) {
  for (std::size_t i = 0; i < 4; ++i) {
    x_[i] = array[i];
  }
}

CFourVector::CFourVector(cdouble x, cdouble y, cdouble z, cdouble t):
    x_({x, y, z, t}) {
}

CFourVector::CFourVector(const TLorentzVector& lvec) {
  for (std::size_t i = 0; i < 4; ++i) {
    x_[i] = lvec[i];
  }
}

CFourVector::CFourVector(const TLorentzVector& lvec, double scale) {
  for (std::size_t i = 0; i < 4; ++i) {
    x_[i] = lvec[i] * scale;
  }
}

cdouble CFourVector::getM2() const {
  return *this * *this;
}

cdouble& CFourVector::operator[](std::size_t i) {
  return x_[i];
}

const cdouble& CFourVector::operator[](std::size_t i) const {
  return x_[i];
}

CFourVector operator + (const CFourVector& v1, const CFourVector& v2) {
  cdouble array[4];
  for (std::size_t i = 0; i < 4; ++i) {
    array[i] = v1[i] + v2[i];
  }
  return CFourVector(array);
}

CFourVector CFourVector::operator - (const CFourVector& v) const {
  cdouble array[4];
  return *this + (-v);
}

cdouble operator * (const CFourVector& v1, const CFourVector& v2) {
  cdouble res = v1.x_[3] * v2.x_[3];
  for (std::size_t i = 0; i < 3; ++i) {
    res -= v1.x_[i] * v2.x_[i];
  }
  return res;
}

CFourVector operator * (const CFourVector& v1, const cdouble& z) {
  cdouble array[4];
  for (std::size_t i = 0; i < 4; ++i) {
    array[i] = v1.x_[i] * z;
  }
  return CFourVector(array);
}

CFourVector operator * (const cdouble& z, const CFourVector& v1) {
  return v1 * z;
}

CFourVector operator / (const CFourVector& v1, const cdouble& z) {
  cdouble array[4];
  for (std::size_t i = 0; i < 4; ++i) {
    array[i] = v1.x_[i] / z;
  }
  return CFourVector(array);
}

CFourVector& CFourVector::operator += (const CFourVector& v) {
  *this = *this + v;
  return *this;
}

CFourVector& CFourVector::operator -= (const CFourVector& v) {
  *this = *this - v;
  return *this;
}

CFourVector& CFourVector::operator *= (const cdouble& z) {
  *this = *this * z;
  return *this;
}

CFourVector CFourVector::operator - () const {
  cdouble array[4];
  for (std::size_t i = 0; i < 4; ++i) {
    array[i] = -this->x_[i];
  }
  return CFourVector(array);
}

std::ostream & operator << (std::ostream& out, const CFourVector& v) {
  return
    out << "(" <<
    v.x_[0] << ", " <<
    v.x_[1] << ", " <<
    v.x_[2] << ", " <<
    v.x_[3] << ")";
}

double CFourVector::getCosTheta() const {
  return x_[2].real() /
      std::sqrt(abs(x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2]));
}

cdouble CFourVector::getX() const {
  return x_[0];
}

cdouble CFourVector::getY() const {
  return x_[1];
}

cdouble CFourVector::getZ() const {
  return x_[2];
}

cdouble CFourVector::getE() const {
  return x_[3];
}

double electron_current_conv(
    const CFourVector& p_minus,
    const CFourVector& p_plus,
    const CFourVector& v) {
  cdouble p1p2 = p_minus * p_plus;
  cdouble g;
  cdouble res = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      g = 0;
      if (i == j && i == 3) {
        g = 1;
      }
      if (i == j && i != 3) {
        g = -1;
      }
      res += (p_minus[i] * p_plus[j] +
              p_minus[j] * p_plus[i] -
              p1p2 * g) * v[i] * conj(v[j]);
    }
  }
  return res.real();
}

double phase_space_weight(
    const std::vector<CFourVector>& momenta) {
  double t_weight = 1.;
  CFourVector p_n = std::accumulate(
      momenta.begin(), momenta.end(), CFourVector(0., 0., 0., 0.));
  CFourVector p_m;
  std::vector<double> masses;
  masses.reserve(momenta.size());
  for (const auto& p : momenta) {
    masses.push_back(sqrt(p.getM2().real()));
  }
  double mu_n = std::accumulate(masses.begin(), masses.end(), 0.);
  double m_n = sqrt(p_n.getM2().real());
  double t_n = m_n - mu_n;
  double mu_m;
  double m_m;
  double t_m;
  double omega_k;
  double p_k;
  double koeff = 0.5 * pow(M_PI, 1.5 * (momenta.size() - 1));
  koeff /= TMath::Gamma(1.5 * (momenta.size() - 1));
  koeff *= pow(t_n, 1.5 * momenta.size() - 2.5) / m_n;
  int n = momenta.size() - 1;
  for (int i = n; i > 0; --i) {
    p_m = p_n - momenta[i];
    mu_m = mu_n - masses[i];
    m_m = sqrt(p_m.getM2().real());
    t_m = m_m - mu_m;
    omega_k = 0.5 * (m_n * m_n + masses[i] * masses[i] - m_m * m_m) / m_n;
    p_k = sqrt(omega_k * omega_k - masses[i] * masses[i]);
    t_weight *= p_k / sqrt(t_n - t_m);
    p_n = p_m;
    mu_n = mu_m;
    m_n = m_m;
    t_n = t_m;
  }
  return t_weight * koeff;
}

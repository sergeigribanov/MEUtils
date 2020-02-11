#ifndef MElParticle_H
#define MElParticle_H

#include <complex>

class MElParticle {
 public:
  explicit MElParticle(int);
  MElParticle(int, double (*)(double, double, double, double*));
  MElParticle(double, double);
  MElParticle(double, double, double (*)(double, double, double, double*));
  MElParticle(const MElParticle&);
  virtual ~MElParticle();
  double getMass() const;
  double getWidthOnMass() const;
  double getWidth(double, double* = nullptr) const;
  std::complex<double> getPropagator(const double&, double* = nullptr) const;

 private:
  double mass_;
  double widthOnMass_;
  double (*widthFCN_)(double, double, double, double*);
};

#endif

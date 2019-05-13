#include <TDatabasePDG.h>
#include "MElParticle.h"

MElParticle::MElParticle(int particleId):
  mass_(TDatabasePDG::Instance()->GetParticle(particleId)->Mass()),
  widthOnMass_(TDatabasePDG::Instance()->GetParticle(particleId)->Width()), 
  widthFCN_(nullptr) {
}

MElParticle::MElParticle(int particleId, double (*widthFCN)(double, double, double, double*)):
mass_(TDatabasePDG::Instance()->GetParticle(particleId)->Mass()),
widthOnMass_(TDatabasePDG::Instance()->GetParticle(particleId)->Width()),
widthFCN_(widthFCN) {
}

MElParticle::MElParticle(double mass, double widthOnMass):
mass_(mass),
widthOnMass_(widthOnMass),
widthFCN_(nullptr) {
}

MElParticle::MElParticle(double mass, double width, double (*widthFCN)(double, double, double, double*)):
mass_(mass),
widthOnMass_(width),
widthFCN_(widthFCN) {
}

MElParticle::MElParticle(const MElParticle& cmd3Particle):
  mass_(cmd3Particle.mass_),
  widthOnMass_(cmd3Particle.widthOnMass_) {
}

MElParticle::~MElParticle() {
}

double MElParticle::getMass() const {
  return mass_;
}

double MElParticle::getWidthOnMass() const {
  return widthOnMass_;
}

double MElParticle::getWidth(double s, double* params) const {
  if (widthFCN_) {
    return (*widthFCN_)(s, mass_, widthOnMass_, params);
  }
  return widthOnMass_;
}

std::complex<double> MElParticle::getPropagator(
    const double& q2, double* params) const {
  return 1. / (q2 - mass_ * mass_ + std::complex<double>(0, 1) *
               getWidth(q2, params));
}

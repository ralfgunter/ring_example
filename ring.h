#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include "pup_stl.h"
#include "liveViz.h"

using namespace std;

const double K = 1e-9;

struct Particle {
    double x, vel, charge;

    double computeForce(Particle &p);
    void pup(PUP::er &p);
};

#include "ring.decl.h"

class Main : public CBase_Main {
  public:
    Main(CkArgMsg *m);
    void cycleDone();

  private:
    CProxy_RingArc ringProxy_;
};

class RingArc : public CBase_RingArc {
  public:
    RingArc(int numSlices, int maxNumParticlesPerSlice);
    RingArc(CkMigrateMessage *m);
    ~RingArc();

    void passMsg(int srcIndex, int otherNumParticles, Particle *ps, double *forces);
    void startCycle();

    void pup(PUP::er &p);

  private:
    int numSlices_;
    vector<Particle> particles_;

    void computeExternalForces(int otherNumParticles, Particle *ps, double *forces);
    void computeMyOwnForces(double *forces);
};

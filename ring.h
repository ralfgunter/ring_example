#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include "pup_stl.h"
#include "liveViz.h"

using namespace std;

const double K = 1e-9;
const int W = 1920, H = 50;
const double convergenceThreshold = 1e-17;

struct Particle {
    double updatePosition(double force);
    double computeForce(Particle &p);
    void pup(PUP::er &p);

    double x, vel, charge;
};

#include "ring.decl.h"

class ForceMsg : public CMessage_ForceMsg {
  public:
    ForceMsg(vector<Particle> &ps, int idx);

    Particle *particles_;
    double *forces_;

    int idx_;
    int numParticles_;
};

class Main : public CBase_Main {
  public:
    Main(CkArgMsg *m);
    void cycleDone(double maxChange);

  private:
    CProxy_RingArc ringProxy_;
};

class RingArc : public CBase_RingArc {
  public:
    RingArc();
    RingArc(CkMigrateMessage *m);
    ~RingArc();

    void passMsg(ForceMsg *m);
    void startCycles();

    void recvParticles(vector<Particle> &ps);

    void requestNextFrame(liveVizRequestMsg *m);
    void pup(PUP::er &p);

  private:
    int w_;
    unsigned char *intensity_;

    double sliceOffset_;

    vector<Particle> particles_;
    vector<Particle> nextRoundParticles_;

    void buildNextFrame();
    double updateParticles(ForceMsg *m);
    void migrateParticles(ForceMsg *m);
    void computePairwiseForces(ForceMsg *m);
    void computeMyOwnForces(ForceMsg *m);
};

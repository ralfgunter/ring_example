#include "ring.h"

/* readonly */ CProxy_Main mainProxy;

/******************************* Main *******************************/
Main::Main(CkArgMsg *m)
{
    int numSlices, maxNumParticlesPerSlice;
    if (m->argc == 3) {
        numSlices = atoi(m->argv[1]);
        maxNumParticlesPerSlice = atoi(m->argv[2]);
    } else {
        numSlices = 3;
        maxNumParticlesPerSlice = 10;
    }

    mainProxy = thisProxy;

    ringProxy_ = CProxy_RingArc::ckNew(numSlices, maxNumParticlesPerSlice, numSlices);
    ringProxy_.startCycle();
}

void Main::cycleDone()
{
    ckout << "Done" << endl;
    CkExit();
}


/******************************* RingArc *******************************/
RingArc::RingArc(int numSlices, int maxNumParticlesPerSlice)
    : numSlices_(numSlices)
{
    double sliceLength = 1.0 / numSlices;
    double sliceOffset = thisIndex * sliceLength;

    // Make sure each chare has a different random seed.
    srand(time(NULL) + thisIndex);

    int numParticles = rand() % maxNumParticlesPerSlice;
    particles_.resize(numParticles);

    for (int i = 0; i < numParticles; ++i) {
        double coord  = ((double) rand()) / RAND_MAX;
        // double charge = ((double) rand()) / RAND_MAX;
        double charge = 1.0;
        double sign = 1.0 - (2*(rand() % 2));
        // double sign = 1.0;

        particles_[i].x = sliceLength * coord + sliceOffset;
        particles_[i].charge = sign * charge;
    }
}

RingArc::RingArc(CkMigrateMessage *m) {}
RingArc::~RingArc() {}

void RingArc::computeExternalForces(int otherNumParticles, Particle *ps, double *forces)
{
    for (int i = 0; i < otherNumParticles; ++i)
        for (int j = 0; j < particles_.size(); ++j)
            forces[i] += particles_[j].computeForce(ps[i]);
}

void RingArc::computeMyOwnForces(double *forces)
{
    for (int i = 0; i < particles_.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            double force = particles_[i].computeForce(particles_[j]);
            forces[j] +=  force;
            forces[i] += -force;
        }
    }
}

void RingArc::passMsg(int srcIndex, int otherNumParticles, Particle *ps, double *forces)
{
    if (srcIndex == thisIndex) {
        // My cycle is complete.
        CkCallback cb(CkReductionTarget(Main, cycleDone), mainProxy);
        contribute(0, NULL, CkReduction::nop, cb);
    } else {
        // I'm in the middle of someone else's cycle; compute my force contribution.
        computeExternalForces(otherNumParticles, ps, forces);

        int nextIndex = (thisIndex + 1) % numSlices_;
        thisProxy(nextIndex).passMsg(srcIndex, otherNumParticles, ps, forces);
    }
}

void RingArc::startCycle()
{
    double *forces = new double[particles_.size()]();
    computeMyOwnForces(forces);

    int nextIndex = (thisIndex + 1) % numSlices_;
    thisProxy(nextIndex).passMsg(thisIndex, particles_.size(), &particles_[0], forces);
}

void RingArc::pup(PUP::er &p)
{
    CBase_RingArc::pup(p);
    p | particles_;
}

/******************************* Particles *******************************/
void Particle::pup(PUP::er &p)
{
    p | x;
    p | vel;
    p | charge;
}

// Compute the force I exert on p.
double Particle::computeForce(Particle &p)
{
    double dist = abs(x - p.x);
    double force = K * charge * p.charge * (pow(dist, -2.0) - pow(1-dist, -2.0));

    int sign = (p.x > x) ? 1 : -1;

    return sign * force;
}

#include "ring.def.h"

#include "ring.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ int numSlices;
/* readonly */ int maxNumParticlesPerSlice;
/* readonly */ double sliceLength;
/* readonly */ double speedLimit;
/* readonly */ double timestep;

Main::Main(CkArgMsg *m)
{
    if (m->argc == 3) {
        numSlices = atoi(m->argv[1]);
        maxNumParticlesPerSlice = atoi(m->argv[2]);
    } else {
        numSlices = 3;
        maxNumParticlesPerSlice = 10;
    }

    timestep = 0.0001;
    sliceLength = 1.0 / numSlices;
    speedLimit = sliceLength / timestep;

    mainProxy = thisProxy;

    CkArrayOptions opts(numSlices);
    ringProxy_ = CProxy_RingArc::ckNew(opts);

    CkCallback c(CkIndex_RingArc::requestNextFrame(0), ringProxy_);
    liveVizConfig cfg(liveVizConfig::pix_color, true);
    liveVizInit(cfg, ringProxy_, c, opts);

    ringProxy_.startCycles();
}

void Main::cycleDone(double maxChange)
{
    if (maxChange > convergenceThreshold) {
        ringProxy_.startCycles();
    } else {
        ckout << "Done" << endl;
        CkExit();
    }
}


RingArc::RingArc()
{
    sliceOffset_ = thisIndex * sliceLength;

    // Make sure each chare has a different random seed.
    srand(time(NULL) + thisIndex);

    int numParticles = rand() % maxNumParticlesPerSlice;

    nextRoundParticles_.resize(numParticles);

    for (int i = 0; i < numParticles; ++i) {
        double coord  = ((double) rand()) / RAND_MAX;
        // double charge = ((double) rand()) / RAND_MAX;
        double charge = 1.0;
        double sign = 1.0 - (2*(rand() % 2));
        // double sign = 1.0;

        nextRoundParticles_[i].x = sliceLength * coord + sliceOffset_;
        nextRoundParticles_[i].charge = sign * charge;
    }

    w_ = W / numSlices;
    intensity_ = new unsigned char[3*w_*H]();
}

RingArc::RingArc(CkMigrateMessage *m) {}
RingArc::~RingArc() { delete[] intensity_; }

void RingArc::computePairwiseForces(ForceMsg *m)
{
    for (int i = 0; i < m->numParticles_; ++i)
        for (int j = 0; j < particles_.size(); ++j)
            m->forces_[i] += m->particles_[i].computeForce(particles_[j]);
}

void RingArc::computeMyOwnForces(ForceMsg *m)
{
    for (int i = 0; i < particles_.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            double force = particles_[i].computeForce(particles_[j]);
            m->forces_[j] +=  force;
            m->forces_[i] += -force;
        }
    }
}

double RingArc::updateParticles(ForceMsg *m)
{
    double maxChange = 0.0;
    for (int j = 0; j < particles_.size(); ++j) {
        maxChange = max(particles_[j].updatePosition(m->forces_[j]), maxChange);
    }
    return maxChange;
}

void RingArc::migrateParticles(ForceMsg *m)
{
    vector<Particle> prev, next;
    for (int j = 0; j < particles_.size(); ++j) {
        if (particles_[j].x < sliceOffset_) {
            // Account for wrap-around.
            if (particles_[j].x < 0.0)
                particles_[j].x += 1.0;

            prev.push_back(particles_[j]);
        } else if (particles_[j].x >= sliceOffset_ + sliceLength) {
            // Account for wrap-around.
            if (particles_[j].x >= 1.0)
                particles_[j].x -= 1.0;

            next.push_back(particles_[j]);
        } else {
            nextRoundParticles_.push_back(particles_[j]);
        }
    }

    // Send particles to my neighbors.
    int nextIndex = (thisIndex + 1) % numSlices;
    int prevIndex = (thisIndex - 1 + numSlices) % numSlices;
    thisProxy(nextIndex).recvParticles(next);
    thisProxy(prevIndex).recvParticles(prev);
}

void RingArc::passMsg(ForceMsg *m)
{
    if (m->idx_ == thisIndex) {
        // My cycle is complete.

        // Act upon particles using the accumulated forces.
        double maxChange = updateParticles(m);

        // Decide which particles stay and which migrate.
        migrateParticles(m);

        CkCallback cb(CkReductionTarget(Main, cycleDone), mainProxy);
        contribute(sizeof(double), &maxChange, CkReduction::max_double, cb);

        delete m;
    } else {
        // I'm in the middle of someone else's cycle.

        computePairwiseForces(m);

        int nextIndex = (thisIndex + 1) % numSlices;
        thisProxy(nextIndex).passMsg(m);
    }
}

void RingArc::startCycles()
{
    particles_.swap(nextRoundParticles_);
    nextRoundParticles_.clear();

    buildNextFrame();

    ForceMsg *m = new (particles_.size(), particles_.size()) ForceMsg(particles_, thisIndex);

    computeMyOwnForces(m);

    int nextIndex = (thisIndex + 1) % numSlices;
    thisProxy(nextIndex).passMsg(m);
}

void RingArc::recvParticles(vector<Particle> &ps)
{
    nextRoundParticles_.insert(nextRoundParticles_.end(), ps.begin(), ps.end());
}

void RingArc::pup(PUP::er &p)
{
    CBase_RingArc::pup(p);

    p | sliceOffset_;
    p | particles_;
    p | nextRoundParticles_;

    p | w_;
    if (p.isUnpacking())
        intensity_ = new unsigned char[3 * w_ * H];
    PUParray(p, intensity_, 3 * w_ * H);
}

void RingArc::requestNextFrame(liveVizRequestMsg *m)
{
    int sx = thisIndex * w_;
    liveVizDeposit(m, sx, 0, w_, H, intensity_, this);
}

void RingArc::buildNextFrame()
{
    fill_n(intensity_, 3 * w_ * H, 0);

    // // Chare delimiters.
    // for (int h = 0; h < H; ++h) {
    //     intensity_[3*(h*w_+0)+0] = 150;
    //     intensity_[3*(h*w_+0)+1] = 150;
    //     intensity_[3*(h*w_+0)+2] = 150;
    // }

    int h = H/2;
    for (int n = 0; n < particles_.size(); ++n) {
        int j = static_cast<int>(floor(((double) w_) * (particles_[n].x - sliceOffset_) / sliceLength));
        if (j < 0) {
            CkPrintf("(x,j) = (%f,%d)\n", particles_[n].x, j);
        }
        int sign = particles_[n].charge > 0;
        intensity_[3*(h*w_+j)+sign] = 255;

        assert(3*(h*w_+j)+sign >= 0);
        assert(3*(h*w_+j)+sign < 3*w_*H);
    }
}

double Particle::updatePosition(double force)
{
    double dx = vel * timestep + (force * timestep*timestep / 2.0);
    dx = max(-sliceLength, min(sliceLength, dx));
    vel = max(-speedLimit, min(speedLimit, vel + force * timestep));
    assert(vel <=  speedLimit);
    assert(vel >= -speedLimit);
    x += dx;
    return dx;
}

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

ForceMsg::ForceMsg(vector<Particle> &ps, int idx) : numParticles_(ps.size()), idx_(idx)
{
    copy(ps.begin(), ps.end(), particles_);
    fill_n(forces_, numParticles_, 0.0);
}

#include "ring.def.h"

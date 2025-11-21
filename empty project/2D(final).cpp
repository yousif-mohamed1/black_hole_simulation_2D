// blackhole_2d_interactive_placement_optionA.cpp
// Interactive 2D black hole lens (Option A) - placement toolbox + auto-fall particles
//
// Controls (quick):
//  P        : toggle placement mode
//  LMB      : when placement OFF -> activate next ray + trigger particle fall
//             when placement ON  -> place particle at cursor
//  RMB      : when placement OFF -> burst 10 rays
//             when placement ON  -> spawn a ray at cursor Y
//  1 / 2 / 3 / 4 : quick spawn particle at Top / Bottom / Left / Right edges
//  Arrow keys : while placement ON, nudge last placed particle
//  UP/DOWN : increase / decrease bendingScale
//  [ / ]   : decrease / increase rayCount (rebuild rays)
//  SPACE   : toggle ray animation
//  C       : clear (deactivate) all rays
//  R       : rebuild rays
//  B       : toggle beam visual mode
//  N / M   : increase / decrease ray animation speed
//  Esc     : exit

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <GLFW/glfw3.h>
#ifdef _WIN32
#pragma comment(lib, "opengl32.lib")
#endif

// ---------------- Configuration / Constants ----------------
static const float G = 1.0f;
static const float M = 0.5f;
static const float C = 1.0f;
static const float Rs = 2.0f * G * M / (C * C);
static const float PhotonSphereR = 1.5f * Rs;
static const float ShadowR = 2.6f * Rs;
static const float PI = 3.14159265358979323846f;

static float DiskInner = 3.0f * Rs;
static float DiskOuter = 7.0f * Rs;

static int gWidth = 1280;
static int gHeight = 720;
static float viewScale = 10.0f;

static float bendingScale = 1.0f;

static int rayCount = 160;
static float rayStartX = -viewScale * 0.95f;
static float rayEndX = viewScale * 0.95f;
static float rayStep = 0.008f;
static int   maxSteps = 12000;

// particle auto-fall: after this many seconds a particle will begin to fall
static float lifeToFall = 8.0f;

struct Vec2 { float x, y; Vec2() :x(0), y(0) {} Vec2(float X, float Y) :x(X), y(Y) {} };
static float lengthV(const Vec2& v) { return std::sqrt(v.x * v.x + v.y * v.y); }
static Vec2 normalizeV(const Vec2& v) { float L = lengthV(v); return (L > 1e-8f) ? Vec2(v.x / L, v.y / L) : Vec2(0, 0); }

// ---------------- Particles ----------------
struct Particle {
    float radius;
    float angularVel;
    float phase;
    float r, g, b;
    bool inside;
    bool falling;
    Vec2 pos;
    float life; // seconds since spawn
};
static std::vector<Particle> gParticles;
static size_t gNextFallIndex = 0;
static int gLastPlacedIndex = -1;

// initialize starting particles (photon ring + inner)
static void initParticles() {
    if (!gParticles.empty()) return;
    int photonRingCount = 60;
    for (int i = 0;i < photonRingCount;i++) {
        float t = (float)i / photonRingCount;
        Particle p;
        p.radius = PhotonSphereR; p.angularVel = 0.6f; p.phase = t * 2.0f * PI;
        p.r = 1.0f; p.g = 0.9f; p.b = 0.3f; p.inside = false; p.falling = false;
        p.pos = Vec2(p.radius * std::cos(p.phase), p.radius * std::sin(p.phase));
        p.life = 0.0f;
        gParticles.push_back(p);
    }
    int innerCount = 30;
    for (int i = 0;i < innerCount;i++) {
        float t = (float)i / innerCount;
        Particle p;
        p.radius = Rs * (0.2f + 0.8f * t);
        p.angularVel = 0.8f + 0.5f * t;
        p.phase = t * 2.0f * PI;
        p.r = 1.0f; p.g = 0.3f + 0.4f * (1.0f - t); p.b = 0.2f; p.inside = true; p.falling = false;
        p.pos = Vec2(p.radius * std::cos(p.phase), p.radius * std::sin(p.phase));
        p.life = 0.0f;
        gParticles.push_back(p);
    }
}

// trigger next non-inside particle to fall (manual)
static void triggerNextFall() {
    while (gNextFallIndex < gParticles.size()) {
        Particle& p = gParticles[gNextFallIndex];
        if (!p.inside && !p.falling) { p.falling = true; gNextFallIndex++; break; }
        gNextFallIndex++;
    }
}

// update particles: orbiting motion and falling behavior
static void updateParticles(float dt) {
    for (size_t i = 0;i < gParticles.size();++i) {
        Particle& p = gParticles[i];
        p.life += dt;
        // Option A: auto-fall if life exceeds threshold and not already inside/falling
        if (!p.inside && !p.falling && p.life >= lifeToFall) {
            p.falling = true;
        }

        if (p.falling) {
            // spiral inward over time
            p.radius -= dt * Rs * 0.4f; // fall speed (tweak)
            p.angularVel += dt * 0.2f;
            p.phase += p.angularVel * dt;
            if (p.radius <= Rs * 0.95f) {
                // move to inner swirl (artistic)
                p.radius = Rs * (0.2f + 0.3f * (std::fmod(p.phase, 1.0f)));
                p.inside = true;
                p.falling = false;
                p.r = 0.8f; p.g = 0.2f; p.b = 0.1f;
            }
            p.pos = Vec2(p.radius * std::cos(p.phase), p.radius * std::sin(p.phase));
        }
        else {
            // normal orbiting
            p.phase += p.angularVel * dt;
            p.pos = Vec2(p.radius * std::cos(p.phase), p.radius * std::sin(p.phase));
        }
    }
}

// ---------------- Rays ----------------
struct Ray {
    std::vector<Vec2> pts;
    std::vector<unsigned char> inDisk;
    float progress;
    bool active;
    bool graze;
    int loops;
    float grazeWeight;
};
static std::vector<Ray> gRays;
static size_t gNextRayIndex = 0;
static float rayAnimSpeed = 140.0f;
static bool animateRays = true;
static bool beamMode = true;
static bool adaptiveSampling = true;

// RNG helper
static float frand() { return (float)rand() / (float)RAND_MAX; }

// adaptive sampling near photon sphere
static float sampleY(int i, int N) {
    float halfSpan = ShadowR * 1.75f;
    int Ncore = (int)(N * 0.55f);
    int Ncluster = N - Ncore;
    if (i < Ncore) {
        float t = (float)i / std::max(1, Ncore - 1);
        return (t - 0.5f) * 2.0f * halfSpan;
    }
    else {
        int k = i - Ncore;
        int half = std::max(1, Ncluster / 2);
        float sigma = 0.12f * Rs;
        if (k < half) { float jitter = (frand() * 2.0f - 1.0f) * sigma; return +PhotonSphereR + jitter; }
        else { float jitter = (frand() * 2.0f - 1.0f) * sigma; return -PhotonSphereR + jitter; }
    }
}

// integrate ray (same pseudo-gravity integrator)
static void integrateRay(float impactY, std::vector<Vec2>& outPoints, float& outMinR, float& outPhiTravel, std::vector<unsigned char>* outInDisk) {
    outPoints.clear(); if (outInDisk) outInDisk->clear();
    Vec2 pos(rayStartX, impactY);
    Vec2 vel(1.0f, 0.0f);
    vel = normalizeV(vel);
    float b = std::abs(impactY);
    float deltaThetaApprox = (2.0f * Rs) / std::max(b, 0.001f);
    deltaThetaApprox *= bendingScale;
    outMinR = 1e9f; outPhiTravel = 0.0f;
    float prevPhi = std::atan2(pos.y, pos.x);
    for (int i = 0;i < maxSteps;i++) {
        outPoints.push_back(pos);
        if (outInDisk) { float rr = lengthV(pos); outInDisk->push_back((rr >= DiskInner && rr <= DiskOuter) ? 1u : 0u); }
        if (std::fabs(pos.x) > viewScale * 1.2f || std::fabs(pos.y) > viewScale * 1.2f) break;
        Vec2 rvec = pos; float r = lengthV(rvec);
        outMinR = std::min(outMinR, r);
        if (r < Rs) {
            for (int k = 0;k < 30;k++) { vel = normalizeV(vel); pos.x += vel.x * rayStep * viewScale * 0.12f; pos.y += vel.y * rayStep * viewScale * 0.12f; outPoints.push_back(pos); if (outInDisk) outInDisk->push_back(0); }
            break;
        }
        float phi = std::atan2(pos.y, pos.x);
        float dphi = phi - prevPhi; if (dphi > PI) dphi -= 2.0f * PI; else if (dphi < -PI) dphi += 2.0f * PI;
        outPhiTravel += std::fabs(dphi);
        prevPhi = phi;
        float proximityWeight = std::exp(-(r - b) * (r - b) / (2.0f * (b * b + 1e-4f)));
        float ringFactor = 1.0f + 4.0f * std::exp(-(r - PhotonSphereR) * (r - PhotonSphereR) / (2.0f * (0.18f * Rs) * (0.18f * Rs)));
        Vec2 toward = (r > 1e-6f) ? Vec2(-rvec.x / r, -rvec.y / r) : Vec2(0, 0);
        float localBend = deltaThetaApprox * proximityWeight * rayStep * 0.5f * ringFactor;
        Vec2 vnorm = normalizeV(vel);
        float dotVT = vnorm.x * toward.x + vnorm.y * toward.y;
        Vec2 steer = Vec2(toward.x - vnorm.x * dotVT, toward.y - vnorm.y * dotVT);
        vel.x += steer.x * localBend; vel.y += steer.y * localBend;
        vel = normalizeV(vel);
        float stepScale = viewScale * 0.5f / (1.0f + 2.0f * std::exp(-(r - PhotonSphereR) * (r - PhotonSphereR) / (2.0f * (0.25f * Rs) * (0.25f * Rs))));
        pos.x += vel.x * rayStep * stepScale; pos.y += vel.y * rayStep * stepScale;
    }
}

static Vec2 worldToNDC(const Vec2& p) { return { p.x / viewScale, p.y / viewScale }; }

// ---------- Drawing helpers (GL) ----------
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// Minimal UI: toolbox as rectangles at top-right. We'll compute hit tests using pixel coords.
struct Button {
    float x0, y0, x1, y1; // NDC coords [-1..1]
    std::string label;
    bool active;
};
static std::vector<Button> gButtons;

static void buildToolboxButtons() {
    gButtons.clear();
    // Create a column of small buttons near top-right in NDC coordinates
    float right = 0.92f;
    float top = 0.90f;
    float w = 0.12f; // width in NDC
    float h = 0.06f; // height
    float gap = 0.01f;
    // Place buttons downward
    auto add = [&](const char* label) {
        Button b; b.x0 = right - w; b.x1 = right; b.y1 = top; b.y0 = top - h; b.label = label; b.active = false;
        gButtons.push_back(b);
        top -= (h + gap);
        };
    add("PLACE");   // toggle placement mode
    add("SPAWN TOP");
    add("SPAWN BOT");
    add("SPAWN L");
    add("SPAWN R");
    add("RAY");     // spawn ray at cursor Y
    add("CLEAR");
    add("BEAM");
}

// draw a button rectangle (no text if GLUT missing). color changes when active
static void drawButtonRect(const Button& b) {
    glBegin(GL_QUADS);
    if (b.active) glColor3f(0.2f, 0.75f, 0.2f);
    else glColor3f(0.15f, 0.15f, 0.15f);
    glVertex2f(b.x0, b.y0);
    glVertex2f(b.x1, b.y0);
    glVertex2f(b.x1, b.y1);
    glVertex2f(b.x0, b.y1);
    glEnd();
    // Thin outline
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(b.x0, b.y0); glVertex2f(b.x1, b.y0); glVertex2f(b.x1, b.y1); glVertex2f(b.x0, b.y1);
    glEnd();
    // Note: we do not render text in this simple fallback. If you have GLUT, you can render label text inside.
}

// small HUD: draw the toolbox
static void drawToolboxHUD() {
    // Draw background translucent
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Buttons
    for (const auto& btn : gButtons) {
        drawButtonRect(btn);
    }
    glDisable(GL_BLEND);
}

// Draw primitives
static void drawFilledCircle(float radius, int segments, float rCenter, float gCenter, float bCenter, float rEdge, float gEdge, float bEdge) {
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(rCenter, gCenter, bCenter);
    Vec2 c = worldToNDC({ 0,0 });
    glVertex2f(c.x, c.y);
    for (int i = 0;i <= segments;i++) {
        float a = (float)i / segments * 2.0f * PI;
        float x = radius * std::cos(a); float y = radius * std::sin(a);
        Vec2 ndc = worldToNDC({ x,y });
        glColor3f(rEdge, gEdge, bEdge);
        glVertex2f(ndc.x, ndc.y);
    }
    glEnd();
}
static void drawCircle(float radius, int segments, float r, float g, float b) {
    glColor3f(r, g, b);
    glBegin(GL_LINE_LOOP);
    for (int i = 0;i < segments;i++) {
        float a = (float)i / segments * 2.0f * PI;
        Vec2 ndc = worldToNDC({ radius * std::cos(a), radius * std::sin(a) });
        glVertex2f(ndc.x, ndc.y);
    }
    glEnd();
}

// Accretion disk color helper (same as yours, clamped)
static void computeDiskColor(float r, float phi, float& R, float& G, float& B) {
    Vec2 tangent(-std::sin(phi), std::cos(phi));
    float cosTheta = tangent.x;
    float f = std::max(0.0f, 1.0f - Rs / std::max(r, Rs + 1e-4f));
    float grav = std::pow(f, 2.0f);
    float v = std::sqrt(std::max(0.0f, G * M / std::max(r, 1e-4f)));
    v = std::min(v, 0.99f);
    float gamma = 1.0f / std::sqrt(std::max(1e-4f, 1.0f - v * v));
    float delta = gamma * (1.0f - v * cosTheta);
    float beam = 1.0f / (delta * delta * delta + 1e-6f);
    beam *= (1.0f + 0.6f * cosTheta);
    beam = std::min(beam, 12.0f);
    float baseR = 1.0f, baseG = 0.6f, baseB = 0.2f;
    float blueShift = (cosTheta > 0 ? 1.15f : 0.85f);
    float greenShift = (cosTheta > 0 ? 1.05f : 0.9f);
    float redShift = (cosTheta > 0 ? 0.95f : 1.10f);
    R = baseR * redShift * grav * beam;
    G = baseG * greenShift * grav * beam;
    B = baseB * blueShift * grav * beam;
    R = std::min(R, 1.0f); G = std::min(G, 1.0f); B = std::min(B, 1.0f);
}

// stars
static std::vector<Vec2> gStars;
static void initStars() { if (!gStars.empty()) return; srand(42); int count = 200; for (int i = 0;i < count;i++) { float x = (frand() * 2.0f - 1.0f) * viewScale; float y = (frand() * 2.0f - 1.0f) * viewScale; if (lengthV({ x,y }) < DiskOuter * 1.1f) { i--; continue; } gStars.push_back({ x,y }); } }
static void drawStars() { glPointSize(2.0f); glBegin(GL_POINTS); glColor3f(1, 1, 1); for (const auto& s : gStars) { Vec2 ndc = worldToNDC(s); glVertex2f(ndc.x, ndc.y); } glEnd(); }

// Ray builder (same algorithm)
static void buildRays() {
    if (!gRays.empty()) return;
    srand(1337);
    std::vector<float> ySamples;
    ySamples.reserve((size_t)(rayCount * 2));
    for (int i = 0;i < rayCount;i++) ySamples.push_back(sampleY(i, rayCount));
    int extra = rayCount; float sigmaTight = 0.06f * Rs;
    for (int i = 0;i < extra / 2;i++) { float jitter = (frand() * 2.0f - 1.0f) * sigmaTight; ySamples.push_back(+PhotonSphereR + jitter); }
    for (int i = 0;i < extra / 2;i++) { float jitter = (frand() * 2.0f - 1.0f) * sigmaTight; ySamples.push_back(-PhotonSphereR + jitter); }

    gRays.resize(ySamples.size());
    for (size_t i = 0;i < ySamples.size(); ++i) {
        float y = ySamples[i];
        std::vector<Vec2> pts; float minR = 1e9f; float phiTravel = 0.0f; std::vector<unsigned char> inDisk;
        integrateRay(y, pts, minR, phiTravel, &inDisk);
        Ray r; r.pts = std::move(pts); r.inDisk = std::move(inDisk); r.progress = 0.0f; r.active = false;
        float band = 0.05f * Rs;
        r.graze = std::abs(minR - PhotonSphereR) < band;
        r.loops = (int)std::floor(phiTravel / (2.0f * PI) + 0.001f);
        float s = 0.03f * Rs; float d = (minR - PhotonSphereR);
        r.grazeWeight = std::exp(-(d * d) / (2.0f * s * s));
        gRays[i] = std::move(r);
    }
}

// spawn a ray at impact Y (useful for placement right-click)
static void spawnRayAtY(float y) {
    std::vector<Vec2> pts; float minR = 1e9f; float phiTravel = 0.0f; std::vector<unsigned char> inDisk;
    integrateRay(y, pts, minR, phiTravel, &inDisk);
    Ray r; r.pts = std::move(pts); r.inDisk = std::move(inDisk); r.progress = 0.0f; r.active = true;
    float band = 0.05f * Rs; r.graze = std::abs(minR - PhotonSphereR) < band;
    r.loops = (int)std::floor(phiTravel / (2.0f * PI) + 0.001f);
    float s = 0.03f * Rs; float d = (minR - PhotonSphereR); r.grazeWeight = std::exp(-(d * d) / (2.0f * s * s));
    gRays.push_back(std::move(r));
}

static void activateNextRay() {
    while (gNextRayIndex < gRays.size()) {
        if (!gRays[gNextRayIndex].active) { gRays[gNextRayIndex].active = true; break; }
        gNextRayIndex++;
    }
    gNextRayIndex++;
}
static void activateRayBurst(int N) { for (int i = 0;i < N;i++) activateNextRay(); }
static void clearRays() { for (auto& r : gRays) r.active = false; gNextRayIndex = 0; }

static void updateRays(float dt) { if (!animateRays) return; for (auto& r : gRays) { if (r.active) { r.progress += rayAnimSpeed * dt; if (r.progress > (float)r.pts.size() - 1) r.progress = (float)r.pts.size() - 1; } } }

static float grazingActivity() { float sum = 0; int act = 0; for (const auto& r : gRays) if (r.active) { sum += r.grazeWeight; act++; } return act ? sum / (float)act : 0.0f; }

static void drawActiveRays() {
    for (const auto& r : gRays) {
        if (!r.active) continue;
        size_t count = (size_t)std::floor(r.progress);
        if (count < 2) continue;
        if (beamMode) {
            glLineWidth(1.5f + 3.0f * r.grazeWeight);
            glBegin(GL_LINE_STRIP);
            for (size_t i = 0;i < count && i < r.pts.size(); ++i) {
                Vec2 p = r.pts[i]; Vec2 ndc = worldToNDC(p);
                float tcol = r.graze ? 1.0f : 0.7f;
                if (i < r.inDisk.size() && r.inDisk[i]) glColor3f(1.0f, 0.55f, 0.2f);
                else glColor3f(0.2f * tcol, 0.7f * tcol, 1.0f * tcol);
                glVertex2f(ndc.x, ndc.y);
            }
            glEnd();
            glLineWidth(1.0f);
        }
        else {
            glBegin(GL_LINE_STRIP);
            for (size_t i = 0;i < count && i < r.pts.size(); ++i) { Vec2 ndc = worldToNDC(r.pts[i]); glVertex2f(ndc.x, ndc.y); }
            glEnd();
        }
        if (r.loops >= 1) {
            glColor3f(1.0f, 0.7f, 0.3f);
            glBegin(GL_LINE_STRIP);
            for (size_t i = 0;i < count && i < r.pts.size(); ++i) { if (r.inDisk.empty() || !r.inDisk[i]) { Vec2 ndc = worldToNDC(r.pts[i]); glVertex2f(ndc.x, ndc.y); } }
            glEnd();
        }
    }
}

static void drawShadowRingGlow() {
    float activity = grazingActivity();
    float alphaScale = std::min(0.75f, 0.15f + 1.2f * activity);
    glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    float inner = ShadowR * 0.992f;
    float outer = ShadowR * 1.018f;
    int segs = 540;
    glBegin(GL_TRIANGLE_STRIP);
    for (int i = 0;i <= segs;i++) {
        float a = (float)i / segs * 2.0f * PI;
        float x0 = inner * std::cos(a), y0 = inner * std::sin(a);
        float x1 = outer * std::cos(a), y1 = outer * std::sin(a);
        Vec2 n0 = worldToNDC({ x0,y0 }), n1 = worldToNDC({ x1,y1 });
        glColor4f(0.95f, 0.98f, 1.0f, alphaScale); glVertex2f(n0.x, n0.y);
        glColor4f(0.95f, 0.98f, 1.0f, 0.0f); glVertex2f(n1.x, n1.y);
    }
    glEnd();
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

static void drawLensedDiskOnRays() {
    glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glPointSize(3.5f);
    glBegin(GL_POINTS);
    for (const auto& r : gRays) {
        if (!r.active) continue;
        size_t count = (size_t)std::floor(r.progress);
        for (size_t i = 0;i < count && i < r.pts.size(); ++i) {
            if (i >= r.inDisk.size() || !r.inDisk[i]) continue;
            const Vec2& p = r.pts[i];
            float rr = lengthV(p); if (rr < DiskInner || rr > DiskOuter) continue;
            float phi = std::atan2(p.y, p.x);
            float R, G, B; computeDiskColor(rr, phi, R, G, B);
            float boost = (p.x > 0.0f ? 1.6f : 0.8f);
            glColor4f(R * boost, G * boost, B * boost, 0.9f);
            Vec2 ndc = worldToNDC(p); glVertex2f(ndc.x, ndc.y);
        }
    }
    glEnd();
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

// ---------- Placement mode implementation ----------
static bool gPlaceBeamMode = false;
static float gLastPlacedY = 0.0f;
static float gLastPlacedX = 0.0f;

static void drawPlacementMarker() {
    if (!gPlaceBeamMode) return;
    glColor3f(0.9f, 0.9f, 0.2f);
    Vec2 a = worldToNDC({ gLastPlacedX, gLastPlacedY });
    glBegin(GL_LINES);
    float s = 0.03f;
    glVertex2f(a.x - s, a.y); glVertex2f(a.x + s, a.y);
    glVertex2f(a.x, a.y - s); glVertex2f(a.x, a.y + s);
    glEnd();
}

static void spawnRayAtXY(float x, float y) {
    spawnRayAtY(y); // integrateRay uses impact parameter y (x ignored)
}

static void spawnParticleAtXY(float x, float y) {
    Particle p;
    p.pos = Vec2(x, y);
    p.radius = lengthV(p.pos);
    p.phase = std::atan2(y, x);
    p.angularVel = 0.6f;
    p.r = 1.0f; p.g = 0.8f; p.b = 0.4f;
    p.inside = (p.radius < Rs);
    p.falling = false;
    p.life = 0.0f;
    gLastPlacedIndex = (int)gParticles.size();
    gParticles.push_back(p);
}

static float getWorldXFromCursor(GLFWwindow* w) {
    double mx, my; glfwGetCursorPos(w, &mx, &my);
    float ndcX = (float)((mx / gWidth) * 2.0 - 1.0);
    return ndcX * viewScale;
}
static float getWorldYFromCursor(GLFWwindow* w) {
    double mx, my; glfwGetCursorPos(w, &mx, &my);
    float ndcY = (float)(1.0 - (my / gHeight) * 2.0);
    return ndcY * viewScale;
}

// spawn particle at edges
static void spawnParticleEdge(int edge) {
    float margin = 0.7f;
    if (edge == 1) spawnParticleAtXY(0.0f, viewScale * 0.95f * margin);
    else if (edge == 2) spawnParticleAtXY(0.0f, -viewScale * 0.95f * margin);
    else if (edge == 3) spawnParticleAtXY(-viewScale * 0.95f * margin, 0.0f);
    else if (edge == 4) spawnParticleAtXY(viewScale * 0.95f * margin, 0.0f);
}

// ---- Toolbox button hit testing (pixel -> NDC -> find button) ----
static bool ndcPointInButton(float ndcX, float ndcY, const Button& b) {
    return ndcX >= b.x0 && ndcX <= b.x1 && ndcY >= b.y0 && ndcY <= b.y1;
}
static int hitTestToolbox(GLFWwindow* win) {
    // convert mouse to NDC
    double mx, my; glfwGetCursorPos(win, &mx, &my);
    float ndcX = (float)((mx / gWidth) * 2.0 - 1.0f);
    float ndcY = (float)(1.0 - (my / gHeight) * 2.0);
    for (size_t i = 0;i < gButtons.size();++i) if (ndcPointInButton(ndcX, ndcY, gButtons[i])) return (int)i;
    return -1;
}

static void drawParticles(double time) {
    glPointSize(4.0f);
    glBegin(GL_POINTS);
    for (const auto& p : gParticles) {
        float x = p.pos.x;
        float y = p.pos.y;
        float fade = p.inside ? 0.4f : (p.falling ? 0.8f : 1.0f);
        glColor3f(p.r * fade, p.g * fade, p.b * fade);
        Vec2 ndc = worldToNDC({ x,y });
        glVertex2f(ndc.x, ndc.y);
    }
    glEnd();
}

// Render everything
static void renderScene() {
    glClearColor(0, 0, 0, 1); glClear(GL_COLOR_BUFFER_BIT);
    initStars(); initParticles(); buildToolboxButtons(); if (gRays.empty()) buildRays();

    static double lastTime = glfwGetTime();
    double time = glfwGetTime();
    float dt = (float)(time - lastTime);
    lastTime = time;

    updateParticles(dt);
    updateRays(dt);

    drawStars();

    // draw disk
    int segments = 180; float r0 = DiskInner, r1 = DiskOuter;
    glBegin(GL_TRIANGLE_STRIP);
    for (int i = 0;i <= segments;i++) {
        float a = (float)i / segments * 2.0f * PI;
        for (int ring = 0; ring < 2; ++ring) {
            float r = (ring == 0) ? r0 : r1;
            float R, G, B; computeDiskColor(r, a, R, G, B);
            float fade = 1.0f - 0.6f * ((r - r0) / (r1 - r0));
            R *= fade; G *= fade; B *= fade;
            glColor3f(R, G, B);
            Vec2 ndc = worldToNDC({ r * std::cos(a), r * std::sin(a) });
            glVertex2f(ndc.x, ndc.y);
        }
    }
    glEnd();

    drawLensedDiskOnRays();
    drawShadowRingGlow();
    drawFilledCircle(Rs, 80, 0.4f, 0.0f, 0.0f, 0.9f, 0.1f, 0.0f);
    drawCircle(Rs, 64, 0.2f, 0.0f, 0.0f);
    drawParticles(time);
    drawActiveRays();
    drawPlacementMarker();
    drawToolboxHUD();

    static double lastPrint = 0.0;
    if (time - lastPrint > 1.0) {
        lastPrint = time;
        std::printf("bendingScale=%.2f | activeRays=%zu/%zu | particles=%zu | placeMode=%s | lifeToFall=%.1f\n",
            bendingScale, gNextRayIndex, gRays.size(), gParticles.size(), gPlaceBeamMode ? "ON" : "OFF", lifeToFall);
    }
}

// Input handling: keyboard nudges and toggles
static void processInput(GLFWwindow* win) {
    if (glfwGetKey(win, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(win, 1);
    if (glfwGetKey(win, GLFW_KEY_UP) == GLFW_PRESS) bendingScale *= 1.01f;
    if (glfwGetKey(win, GLFW_KEY_DOWN) == GLFW_PRESS) bendingScale /= 1.01f;
    bendingScale = std::max(0.01f, std::min(bendingScale, 50.0f));

    // nudges for last placed particle when placement ON
    if (gPlaceBeamMode && gLastPlacedIndex >= 0 && gLastPlacedIndex < (int)gParticles.size()) {
        Particle& p = gParticles[gLastPlacedIndex];
        float nudge = 0.02f;
        if (glfwGetKey(win, GLFW_KEY_LEFT) == GLFW_PRESS) p.pos.x -= nudge;
        if (glfwGetKey(win, GLFW_KEY_RIGHT) == GLFW_PRESS) p.pos.x += nudge;
        if (glfwGetKey(win, GLFW_KEY_W) == GLFW_PRESS || glfwGetKey(win, GLFW_KEY_UP) == GLFW_PRESS) p.pos.y += nudge;
        if (glfwGetKey(win, GLFW_KEY_S) == GLFW_PRESS || glfwGetKey(win, GLFW_KEY_DOWN) == GLFW_PRESS) p.pos.y -= nudge;
        p.radius = lengthV(p.pos);
        p.phase = std::atan2(p.pos.y, p.pos.x);
    }

    // toggle placement P (debounce)
    static int prevP = GLFW_RELEASE; int pState = glfwGetKey(win, GLFW_KEY_P);
    if (prevP == GLFW_RELEASE && pState == GLFW_PRESS) { gPlaceBeamMode = !gPlaceBeamMode; std::printf("Placement mode: %s\n", gPlaceBeamMode ? "ON" : "OFF"); }
    prevP = pState;
}

// Mouse callback handles toolbox clicks and placement actions
static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        // first check if clicked a toolbox button
        int hit = hitTestToolbox(window);
        if (hit >= 0) {
            const std::string& lbl = gButtons[hit].label;
            if (lbl == "PLACE") { gPlaceBeamMode = !gPlaceBeamMode; }
            else if (lbl == "SPAWN TOP") { spawnParticleEdge(1); }
            else if (lbl == "SPAWN BOT") { spawnParticleEdge(2); }
            else if (lbl == "SPAWN L") { spawnParticleEdge(3); }
            else if (lbl == "SPAWN R") { spawnParticleEdge(4); }
            else if (lbl == "RAY") { float y = getWorldYFromCursor(window); spawnRayAtY(y); }
            else if (lbl == "CLEAR") { clearRays(); }
            else if (lbl == "BEAM") { beamMode = !beamMode; }
            return;
        }

        // else toolbox miss -> normal placement/ray controls
        if (gPlaceBeamMode) {
            gLastPlacedX = getWorldXFromCursor(window);
            gLastPlacedY = getWorldYFromCursor(window);
            spawnParticleAtXY(gLastPlacedX, gLastPlacedY);
        }
        else {
            triggerNextFall();
            activateNextRay();
        }
    }
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        int hit = hitTestToolbox(window);
        if (hit >= 0) {
            // Right-click on toolbox: do nothing or could be alternate action
            return;
        }
        if (gPlaceBeamMode) {
            float x = getWorldXFromCursor(window);
            float y = getWorldYFromCursor(window);
            spawnRayAtXY(x, y);
        }
        else {
            for (int i = 0;i < 10;i++) activateNextRay();
        }
    }
}

// keyboard callback for other shortcuts
static void keyCallback(GLFWwindow* window, int key, int sc, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_SPACE) animateRays = !animateRays;
        else if (key == GLFW_KEY_C) clearRays();
        else if (key == GLFW_KEY_R) { gRays.clear(); buildRays(); }
        else if (key == GLFW_KEY_B) beamMode = !beamMode;
        else if (key == GLFW_KEY_LEFT_BRACKET) { rayCount = std::max(16, rayCount / 2); gRays.clear(); buildRays(); }
        else if (key == GLFW_KEY_RIGHT_BRACKET) { rayCount = std::min(2000, rayCount * 2); gRays.clear(); buildRays(); }
        else if (key == GLFW_KEY_N) { rayAnimSpeed *= 1.2f; }
        else if (key == GLFW_KEY_M) { rayAnimSpeed = std::max(10.0f, rayAnimSpeed / 1.2f); }
        else if (key == GLFW_KEY_1) spawnParticleEdge(1);
        else if (key == GLFW_KEY_2) spawnParticleEdge(2);
        else if (key == GLFW_KEY_3) spawnParticleEdge(3);
        else if (key == GLFW_KEY_4) spawnParticleEdge(4);
        else if (key == GLFW_KEY_P) { gPlaceBeamMode = !gPlaceBeamMode; std::printf("Placement mode: %s\n", gPlaceBeamMode ? "ON" : "OFF"); }
    }
}

static void framebufferSizeCallback(GLFWwindow* w, int width, int height) { gWidth = width; gHeight = height; glViewport(0, 0, width, height); }

int main() {
    std::printf("Interactive Black Hole 2D Lens (Option A) - placement toolbox active.\n P toggles placement; click toolbox buttons (top-right) or use hotkeys.\n");
    if (!glfwInit()) { std::cerr << "Failed to init GLFW\n"; return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2); glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    GLFWwindow* window = glfwCreateWindow(gWidth, gHeight, "Black Hole 2D - Option A (Placement + Auto-Fall)", nullptr, nullptr);
    if (!window) { std::cerr << "Failed to create window\n"; glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    const GLubyte* ver = glGetString(GL_VERSION);
    const GLubyte* renderer = glGetString(GL_RENDERER);
    std::cout << "GL_VERSION=" << (ver ? (const char*)ver : "?") << " | GL_RENDERER=" << (renderer ? (const char*)renderer : "?") << std::endl;

    glViewport(0, 0, gWidth, gHeight);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(4.0f);

    // prebuild toolbox layout
    buildToolboxButtons();

    while (!glfwWindowShouldClose(window)) {
        processInput(window);
        renderScene();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window); glfwTerminate(); return 0;
}

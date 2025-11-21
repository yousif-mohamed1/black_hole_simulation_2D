# Black Hole 2D Lens Simulation 

## Overview
Minimal interactive 2D visualization of a stylized black hole gravitational lens, accretion disk, photon sphere rays, and orbiting / falling particles. Written in C++14 using raw OpenGL (fixed pipeline) + GLFW for windowing/input.

## Main File
`empty project/2D(final).cpp`

## Core Features
- Pseudo–gravitational bending of light rays around a Schwarzschild black hole
- Photon sphere grazing detection and glow
- Accretion disk with Doppler beaming + gravitational dimming approximation
- Adaptive sampling of initial ray impact parameters
- Interactive activation and animation of rays (beam or simple line mode)
- Particle ring + inner swirl with auto-fall after lifetime or manual trigger
- Simple placement toolbox (clickable HUD buttons) + manual particle placement mode

## Controls
Toolbox buttons (top-right) duplicate hotkeys. Hotkeys:
- `P` : Toggle placement mode
- Left Mouse (placement OFF) : Activate next ray + trigger particle fall
- Left Mouse (placement ON)  : Place particle at cursor
- Right Mouse (placement OFF): Burst activate 10 rays
- Right Mouse (placement ON) : Spawn a single ray at cursor Y
- `1 2 3 4` : Spawn particle at Top / Bottom / Left / Right edges
- Arrow keys (placement ON) : Nudge last placed particle
- `UP / DOWN` : Increase / decrease `bendingScale`
- `[ / ]` : Halve / double `rayCount` (rebuild rays)
- `SPACE` : Toggle ray animation
- `C` : Clear (deactivate) all rays
- `R` : Rebuild rays
- `B` : Toggle beam visual mode
- `N / M` : Increase / decrease ray animation speed
- `Esc` : Exit

Console periodically prints current state (bending scale, rays, particle count, placement mode, life-to-fall).

## Dependencies
- C++14 compiler
- OpenGL 2.1+ (fixed pipeline calls: `glBegin`, etc.)
- GLFW 3.x

## Building (Windows - Visual Studio)
1. Install GLFW (prebuilt binaries or build from source).
2. Add GLFW include path and link library (`glfw3.lib`); ensure system links `opengl32.lib` (already pragma'd on `_WIN32`).
3. Compile `2D(final).cpp` as a Console / Win32 app (C++14).

## Building (Cross-Platform - Example CMake)
Create `CMakeLists.txt`:
```
cmake_minimum_required(VERSION 3.10)
project(black_hole_2d LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
add_executable(black_hole_2d "empty project/2D(final).cpp")
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
if (GLFW_FOUND)
    target_include_directories(black_hole_2d PRIVATE ${GLFW_INCLUDE_DIRS})
    target_link_libraries(black_hole_2d ${GLFW_LIBRARIES})
endif()
# On Linux you may also need: target_link_libraries(black_hole_2d GL X11 pthread dl)
```
Then:
```
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

## Running
Run the produced executable. A window opens with the simulation; interact using the controls above.

## Performance Notes
- `rayCount` and `maxSteps` strongly impact CPU cost; reduce for slower machines.
- `bendingScale` exaggerated for visual clarity; extreme values can cause ray clustering.

## Theory Behind Current Approximation
This MVP does not solve true general relativity geodesics. It applies a heuristic steering force pulling ray velocity toward the origin (black hole) to emulate deflection. Key elements:
- Deflection scale uses a term proportional to `2*Rs / b` (impact parameter b ? |y|) inspired by the weak-field light bending formula ? `4GM/(bc^2)` (constants simplified for visuals).
- A Gaussian proximity weight increases bending near closest approach; a ring factor amplifies deflection near the photon sphere radius (1.5 Rs) to create multiple loops.
- Loop count is estimated by accumulated absolute change in polar angle ?.
- Rays mark samples passing through the disk annulus to highlight lensed disk emission.
- Particles orbit by incrementing angular phase; falling particles reduce radius linearly (visual spiral, not physical accretion).

## Towards a Real Physics Engine
To upgrade from heuristic steering to physically correct light propagation:
1. Use Schwarzschild metric (start with equatorial plane) with units G = c = 1 so Rs = 2M.
2. Integrate null geodesics (photon paths) using differential equations derived from the metric's Christoffel symbols or effective potential.
3. Employ adaptive Runge–Kutta (RK45) for accurate step control near the photon sphere (critical impact parameter region).
4. Compute bending angle analytically for validation: weak-field limit ?? ? 4GM/(bc^2).
5. Replace artificial ringFactor with naturally emergent multiple orbits as impact parameter approaches critical value b_crit.
6. Map geodesic samples back to disk coordinates to compute emission + relativistic Doppler/gravitational shifts.
7. Separate integrators for photons (null) and massive particles (timelike) to simulate orbiting matter near ISCO (innermost stable circular orbit at 3 Rs/2 for light sphere, 3 Rs for massive particle orbit boundary).

### Geodesic Essentials (Equatorial Simplification)
For Schwarzschild in equatorial plane (? = ?/2), with affine parameter ?:
- Conserved quantities: Energy E, angular momentum L. Set E = 1 for scaling.
- Effective potential for photons: `V_eff(r) = (1 - Rs/r) * L^2 / r^2`.
- Radial equation: `(dr/d?)^2 = E^2 - V_eff(r)`.
- Angular evolution: `d?/d? = L / r^2`.
Near the photon sphere (r = 1.5 Rs) small changes in impact parameter produce large increases in ? travel (multiple loops).

### Engine Module Plan
- `math/` : Vector, small matrix, RK45 integrator.
- `metric/` : Schwarzschild metric provider (Christoffel symbol calculations if full 4D later).
- `integrators/` : GeodesicStepper (photon), ParticleStepper (massive).
- `engine/` : RayBuilder (samples impact parameters), DiskModel (emissivity & color), ParticleSystem.
- `render/` : GPU buffers (VBO) for rays, particles, disk, post-processing glow.

### Incremental Implementation Steps
1. Refactor current heuristic integrator behind an interface `IRayIntegrator`.
2. Add RK45 implementation and a Schwarzschild geodesic integrator returning polyline samples.
3. Side-by-side compare heuristic vs physical mode (`--mode=approx` / `--mode=phys`).
4. Add weak-field bending angle test harness (assert error < chosen tolerance for large b).
5. Introduce disk emission physically: emissivity ?(r) ~ r^{-q} (q ? 2–3), gravitational redshift factor ?(1 - Rs/r), Doppler factor using orbital velocity v = ?(GM/r).
6. Optimize: cache geodesic samples for reused impact parameters; multithread precomputation.
7. Replace immediate mode with modern OpenGL (VAO/VBO + GLSL) for performance.

### Key Data Structures (Conceptual)
```
struct GeodesicState { double r, phi, pr, pphi; } // minimal 2D form
struct RaySample { float x, y, r, phi; }          // stored for rendering
struct PhotonRay { std::vector<RaySample> samples; double minR; int loops; bool active; };
```

### Validation & Testing
- Compare numerical bending angle to analytic weak-field formula for large impact parameters.
- Ensure absorption when r < Rs (terminate ray).
- Observe divergence in loop count as b ? b_crit.

### Future Extensions
- Kerr (rotating) black hole metric (requires ? evolution & frame dragging terms).
- Shader-based color grading (tone mapping, bloom around photon sphere).
- Disk self-occlusion & secondary images (higher-order lensing passes).
- Recording geodesic parameters for educational overlay.

## Possible Extensions
- Replace fixed pipeline with modern OpenGL + shaders
- Add true GR lensing integration (geodesic solver)
- UI text rendering for button labels
- Save / load particle layouts
- Frame capture / recording

## License
Add a license of your choice (e.g., MIT) in a `LICENSE` file.

## Disclaimer
Visualization is artistic / approximate, not a physically rigorous simulation. Physics engine roadmap above describes future work, not current implementation.

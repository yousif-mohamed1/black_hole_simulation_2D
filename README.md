# Black Hole 2D Lens Simulation (MVP)

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

## Possible Extensions
- Replace fixed pipeline with modern OpenGL + shaders
- Add true GR lensing integration (geodesic solver)
- UI text rendering for button labels
- Save / load particle layouts
- Frame capture / recording

## License
Add a license of your choice (e.g., MIT) in a `LICENSE` file.

## Disclaimer
Visualization is artistic / approximate, not a physically rigorous simulation.

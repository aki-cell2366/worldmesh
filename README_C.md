# WorldMesh C Library

A complete C99 implementation of the World Grid Square Code (ISO 24108-1 Annex B) for calculating Cartesian grid squares compatible with JIS X 0410 standard.

## Overview

This library provides functions to:
1. **Convert coordinates to grid square codes** - Calculate grid square codes from latitude/longitude pairs at multiple resolution levels
2. **Convert grid square codes to coordinates** - Determine grid square boundaries from grid square codes
3. **Calculate geodesic distances** - Compute distances between points using Vincenty's formulae on WGS84 ellipsoid
4. **Calculate grid areas** - Determine the size and area of grid squares

## Features

- **Multiple resolution levels**: 8, 10, 12, 13, 14, 16-digit grid square codes
- **Extended high-resolution codes**: 100m, 50m, 10m, 1m precision variants
- **WGS84 ellipsoid support**: Accurate geodesic calculations using Vincenty's formulae
- **Robust error handling**: Input validation and boundary condition checking
- **C99 standard**: Portable and widely compatible with modern C compilers
- **No external dependencies**: Only uses standard C library (<math.h>, <string.h>, <stdio.h>)

## Building

### Compile as Static Library

```bash
gcc -c -Wall -std=c99 worldmesh.c -o worldmesh.o -lm
ar rcs libworldmesh.a worldmesh.o
```

### Compile as Shared Library

```bash
gcc -shared -fPIC -Wall -std=c99 -o libworldmesh.so worldmesh.c -lm
```

### With CMake

```bash
mkdir build && cd build
cmake ..
make
make install
```

## Usage

### Basic Example

```c
#include "worldmesh.h"
#include <stdio.h>

int main(void) {
    char meshcode[32];
    
    /* Convert coordinates to grid square mesh code (1km resolution, 10 digits) */
    cal_meshcode3(35.6762, 139.7674, meshcode);
    printf("World grid square code: %s\n", meshcode);
    
    /* Convert grid square code back to coordinates */
    WorldMesh_Grid grid;
    meshcode_to_latlong_grid(meshcode, 0, "", &grid);
    printf("Grid corners:\n");
    printf("  NW: %.6f, %.6f\n", grid.lat0, grid.long0);
    printf("  SE: %.6f, %.6f\n", grid.lat1, grid.long1);
    
    /* Calculate geodesic distance */
    double distance = Vincenty(35.6762, 139.7674, 34.6937, 135.5023);
    printf("Distance: %.2f m (%.2f km)\n", distance, distance / 1000.0);
    
    /* Calculate area of grid square */
    WorldMesh_Area area;
    cal_area_from_meshcode(meshcode, 0, "", &area);
    printf("Grid area: %.2f m²\n", area.A);
    
    return 0;
}
```

### Linking

```bash
gcc -o myapp myapp.c libworldmesh.a -lm
```

## API Reference

### Data Structures

#### WorldMesh_Grid
```c
typedef struct {
    double lat0;   /* Northern latitude (NW corner) */
    double long0;  /* Western longitude (NW corner) */
    double lat1;   /* Southern latitude (SE corner) */
    double long1;  /* Eastern longitude (SE corner) */
} WorldMesh_Grid;
```

#### WorldMesh_Area
```c
typedef struct {
    double W1;  /* Northern span (meters) */
    double W2;  /* Southern span (meters) */
    double H;   /* North-south span (meters) */
    double A;   /* Area (square meters) */
} WorldMesh_Area;
```

### Meshcode to Coordinates Functions

All these functions require the caller to allocate the output structure.

#### meshcode_to_latlong_grid()
```c
WorldMesh_Grid* meshcode_to_latlong_grid(const char *meshcode, int extension,
                                          const char *spec, WorldMesh_Grid *out);
```
Convert mesh code to grid square boundaries. Returns pointer to `out` on success, NULL on failure.

**Parameters:**
- `meshcode`: String of 8, 10, 11, 12, 13, 14, or 16 digits
- `extension`: 1 for extended high-resolution codes, 0 for standard codes
- `spec`: Specification for extended codes:
  - `"ex100m_13"` - 100m resolution (13 digits)
  - `"ex50m_14"` - 50m resolution (14 digits)
  - `"ex10m_14"` - 10m resolution (14 digits)
  - `""` - Default based on digit count and extension mode
- `out`: Pointer to output structure

#### Corner Variants
```c
WorldMesh_Grid* meshcode_to_latlong_NW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);
WorldMesh_Grid* meshcode_to_latlong_SW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);
WorldMesh_Grid* meshcode_to_latlong_NE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);
WorldMesh_Grid* meshcode_to_latlong_SE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);
```
Return a single corner coordinates in `lat0` and `long0` fields.

### Coordinates to Meshcode Functions

All these functions require the caller to provide output buffer (minimum 32 bytes).

#### Standard Resolution Codes

```c
char* cal_meshcode1(double latitude, double longitude, char *out);  /* 6 digits, 80km */
char* cal_meshcode2(double latitude, double longitude, char *out);  /* 8 digits, 10km */
char* cal_meshcode3(double latitude, double longitude, char *out);  /* 10 digits, 1km */
char* cal_meshcode4(double latitude, double longitude, char *out);  /* 11 digits, 500m */
char* cal_meshcode5(double latitude, double longitude, char *out);  /* 12 digits, 250m */
char* cal_meshcode6(double latitude, double longitude, char *out);  /* 13 digits, 125m */
```

All return pointer to `out` on success, or pointer to error code string "999..." on invalid input.

#### Extended High-Resolution Codes

```c
char* cal_meshcode_ex100m_12(double latitude, double longitude, char *out);
char* cal_meshcode_ex100m_13(double latitude, double longitude, char *out);
char* cal_meshcode_ex50m_13(double latitude, double longitude, char *out);
char* cal_meshcode_ex50m_14(double latitude, double longitude, char *out);
char* cal_meshcode_ex10m_14(double latitude, double longitude, char *out);
char* cal_meshcode_ex1m_16(double latitude, double longitude, char *out);
```

### Distance and Area Functions

#### Vincenty()
```c
double Vincenty(double latitude1, double longitude1,
                double latitude2, double longitude2);
```
Calculate geodesic distance between two points on WGS84 ellipsoid.

**Returns:** Distance in meters

#### cal_area_from_meshcode()
```c
WorldMesh_Area* cal_area_from_meshcode(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Area *out);
```
Calculate grid square size and area from mesh code. Returns pointer to `out` on success, NULL on failure.

#### cal_area_from_latlong()
```c
WorldMesh_Area* cal_area_from_latlong(const WorldMesh_Grid *grid,
                                       WorldMesh_Area *out);
```
Calculate grid square size and area from corner coordinates. Returns pointer to `out` on success, NULL on failure.

### Utility Functions

#### worldmesh_validate_coordinates()
```c
int worldmesh_validate_coordinates(double latitude, double longitude);
```
Validate that coordinates are within valid ranges (-90 to 90 latitude, -180 to 180 longitude).

**Returns:** 1 if valid, 0 if invalid

## Grid Square Code Structure

The grid square code is a string of digits encoding hierarchical grid squares:

```
ABBBBBCCDDEFG...
│     │ │ ││
│     │ │ │└─ Extended codes
│     │ │ └── 1km grid (2 digits)
│     │ └──── 10km grid (2 digits)
│     └─────── 80km grid (4 digits)
└────────────── Area code (1 digit: 1-8)
```

### Digit Counts and Resolutions

| Digits | Size | Resolution |
|--------|------|------------|
| 6 | 80km | 40 arc-min × 1 arc-deg |
| 8 | 10km | 5 arc-min × 7.5 arc-min |
| 10 | 1km | 30 arc-sec × 45 arc-sec |
| 11 | 500m | 15 arc-sec × 22.5 arc-sec |
| 12 | 250m | 7.5 arc-sec × 11.25 arc-sec |
| 13 | 125m | 3.75 arc-sec × 5.625 arc-sec |
| 12-16 | 100m-1m | Extended high-resolution variants |

## Testing

Compile and run the test suite:

```bash
gcc -std=c99 -Wall -Wextra worldmesh_test.c worldmesh.o -o worldmesh_test -lm
./worldmesh_test
```

The test suite covers:
- Input validation and boundary conditions
- Coordinate to meshcode conversion at all levels
- Meshcode to coordinate conversion
- Round-trip conversions (coordinate → code → coordinate)
- Area calculations
- Geodesic distance calculations
- Invalid input handling
- Extended high-resolution codes

## Performance Characteristics

- **Time complexity**: O(1) for all conversions
- **Space complexity**: O(1) for coordinate calculations, O(n) for string meshcodes
- **Precision**: IEEE 754 double-precision (64-bit) floating-point
- **Accuracy**: Sub-meter precision for geodesic calculations

## Compatibility

- **Standard**: C99 (ISO/IEC 9899:1999)
- **Platforms**: Linux, macOS, Windows (with MinGW)
- **Dependencies**: libc, libm (standard C math library)
- **Memory**: Stack-based, no dynamic allocation

## References

- **Vincenty's Formulae**: T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of nested equations", Survey Review XXIII, Vol 176 (1975)
- **JIS X 0410:1975**: Japanese Standard for Grid Square Systems
- **ISO 24108-1:2025**: Grid square statistics and their applications - Part 1: Fundamental principle of grid square statistics
- **WGS84**: World Geodetic System 1984 Ellipsoid Parameters

## License

Copyright (c) 2015-2026 Research Institute for World Grid Squares
Prof. Dr. Aki-Hiro Sato (Yokohama City University)

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Authors

- **Original R implementation**: Prof. Dr. Aki-Hiro Sato
- **C port**: 03 May 2026

## Contact

For issues, questions, or contributions, refer to the original research institute:
- Address: 22-2, Seto, Kanazawa-ku, Yokohama, Kanagawa 236-0027 Japan
- E-mail: ahsato@yokohama-cu.ac.jp

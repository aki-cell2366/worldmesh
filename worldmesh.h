/*
 * Copyright (c) 2015-2022 Research Institute for World Grid Squares
 * Prof. Dr. Aki-Hiro Sato (Yokohama City University)
 * All rights reserved.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * C implementation of world grid square code calculation.
 * The world grid square code computed by this library is compatible to JISX0410.
 *
 * Ported from R implementation (worldmesh.R v1.81)
 * C implementation version: 1.0 (May 2026)
 */

#ifndef WORLDMESH_H
#define WORLDMESH_H

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/**
 * WorldMesh_Grid: Represents the four corners of a grid square
 * lat0, long0: Northwest (top-left) corner
 * lat1, long1: Southeast (bottom-right) corner
 */
typedef struct {
    double lat0;   /* Northern latitude */
    double long0;  /* Western longitude */
    double lat1;   /* Southern latitude */
    double long1;  /* Eastern longitude */
} WorldMesh_Grid;

/**
 * WorldMesh_Area: Represents the size and area of a grid square
 * W1: Northern west-to-east span (meters)
 * W2: Southern west-to-east span (meters)
 * H:  North-to-south span (meters)
 * A:  Area approximated by trapezoid (square meters)
 */
typedef struct {
    double W1;  /* Northern span */
    double W2;  /* Southern span */
    double H;   /* North-south span */
    double A;   /* Area */
} WorldMesh_Area;

/* ============================================================================
 * 1. Meshcode to Coordinates Functions
 * ============================================================================ */

/**
 * meshcode_to_latlong_grid:
 * Calculate the four corners of a grid square from a meshcode string.
 *
 * Parameters:
 *   meshcode  - Meshcode string (6 to 16 digits)
 *   extension - If 1, interpret extended high-resolution codes
 *   spec      - Specification for extended codes: "ex100m_13", "ex50m_14", etc.
 *   out       - Output structure pointer (caller allocates)
 *
 * Returns:
 *   Pointer to out on success, NULL on failure
 */
WorldMesh_Grid* meshcode_to_latlong_grid(const char *meshcode, int extension,
                                          const char *spec, WorldMesh_Grid *out);

/**
 * meshcode_to_latlong_NW:
 * Calculate the northwest (top-left) corner of a grid square.
 */
WorldMesh_Grid* meshcode_to_latlong_NW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);

/**
 * meshcode_to_latlong_SW:
 * Calculate the southwest (bottom-left) corner of a grid square.
 */
WorldMesh_Grid* meshcode_to_latlong_SW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);

/**
 * meshcode_to_latlong_NE:
 * Calculate the northeast (top-right) corner of a grid square.
 */
WorldMesh_Grid* meshcode_to_latlong_NE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);

/**
 * meshcode_to_latlong_SE:
 * Calculate the southeast (bottom-right) corner of a grid square.
 */
WorldMesh_Grid* meshcode_to_latlong_SE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out);

/* ============================================================================
 * 2. Coordinates to Meshcode Functions
 * ============================================================================ */

/**
 * cal_meshcode1:
 * Calculate 80km grid square code (6 digits) from latitude, longitude.
 */
char* cal_meshcode1(double latitude, double longitude, char *out);

/**
 * cal_meshcode2:
 * Calculate 10km grid square code (8 digits) from latitude, longitude.
 */
char* cal_meshcode2(double latitude, double longitude, char *out);

/**
 * cal_meshcode3:
 * Calculate 1km grid square code (10 digits) from latitude, longitude.
 * (Default: cal_meshcode)
 */
char* cal_meshcode3(double latitude, double longitude, char *out);

/**
 * cal_meshcode4:
 * Calculate 500m grid square code (11 digits) from latitude, longitude.
 */
char* cal_meshcode4(double latitude, double longitude, char *out);

/**
 * cal_meshcode5:
 * Calculate 250m grid square code (12 digits) from latitude, longitude.
 */
char* cal_meshcode5(double latitude, double longitude, char *out);

/**
 * cal_meshcode6:
 * Calculate 125m grid square code (13 digits) from latitude, longitude.
 */
char* cal_meshcode6(double latitude, double longitude, char *out);

/* Extended high-resolution codes (not in JIS X0410 directly) */

/**
 * cal_meshcode_ex100m_12:
 * Calculate extended 100m grid code (12 digits) - 3 arc-sec x 4.5 arc-sec.
 */
char* cal_meshcode_ex100m_12(double latitude, double longitude, char *out);

/**
 * cal_meshcode_ex100m_13:
 * Calculate extended 100m grid code (13 digits) - 3 arc-sec x 4.5 arc-sec.
 */
char* cal_meshcode_ex100m_13(double latitude, double longitude, char *out);

/**
 * cal_meshcode_ex50m_13:
 * Calculate extended 50m grid code (13 digits) - 1.5 arc-sec x 2.25 arc-sec.
 */
char* cal_meshcode_ex50m_13(double latitude, double longitude, char *out);

/**
 * cal_meshcode_ex50m_14:
 * Calculate extended 50m grid code (14 digits) - 1.5 arc-sec x 2.25 arc-sec.
 */
char* cal_meshcode_ex50m_14(double latitude, double longitude, char *out);

/**
 * cal_meshcode_ex10m_14:
 * Calculate extended 10m grid code (14 digits) - 0.3 arc-sec x 0.45 arc-sec.
 */
char* cal_meshcode_ex10m_14(double latitude, double longitude, char *out);

/**
 * cal_meshcode_ex1m_16:
 * Calculate extended 1m grid code (16 digits) - 0.03 arc-sec x 0.045 arc-sec.
 */
char* cal_meshcode_ex1m_16(double latitude, double longitude, char *out);

/* ============================================================================
 * 3. Geodesic Distance and Area Functions
 * ============================================================================ */

/**
 * Vincenty:
 * Calculate geodesic distance between two points on WGS84 ellipsoid.
 *
 * Based on Vincenty's formulae (1975):
 *   T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
 *   with application of nested equations", Survey Review XXIII, Vol 176 (1975)
 *
 * Parameters:
 *   latitude1, longitude1  - First point (degrees)
 *   latitude2, longitude2  - Second point (degrees)
 *
 * Returns:
 *   Distance in meters
 */
double Vincenty(double latitude1, double longitude1,
                double latitude2, double longitude2);

/**
 * cal_area_from_meshcode:
 * Calculate size and area of a grid square from its meshcode.
 *
 * Parameters:
 *   meshcode  - Meshcode string
 *   extension - If 1, interpret extended codes
 *   spec      - Specification string for extended codes
 *   out       - Output structure pointer (caller allocates)
 *
 * Returns:
 *   Pointer to out on success, NULL on failure
 */
WorldMesh_Area* cal_area_from_meshcode(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Area *out);

/**
 * cal_area_from_latlong:
 * Calculate size and area of a trapezoid defined by four corners.
 *
 * Parameters:
 *   grid - Coordinates structure with lat0, long0, lat1, long1
 *   out  - Output area structure pointer (caller allocates)
 *
 * Returns:
 *   Pointer to out on success, NULL on failure
 */
WorldMesh_Area* cal_area_from_latlong(const WorldMesh_Grid *grid,
                                       WorldMesh_Area *out);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Validate latitude and longitude bounds.
 * Returns 1 if valid, 0 if invalid.
 */
int worldmesh_validate_coordinates(double latitude, double longitude);

#ifdef __cplusplus
}
#endif

#endif /* WORLDMESH_H */

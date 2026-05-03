/*
 * Copyright (c) 2015-2022 Research Institute for World Grid Squares
 * Prof. Dr. Aki-Hiro Sato (Yokohama City University)
 * All rights reserved.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * C implementation of world grid square code calculation.
 * Ported from R implementation (worldmesh.R v1.81)
 * The world grid square code computed by this library 
 * is compatible to Japanese grid square code (JIS X 0410). 
 *
 * This library for World grid square code complies ISO 24108-1:
 * https://www.iso.org/standard/87706.html
 *
 * Written by Prof. Dr. Aki-Hiro Sato
 * Graduate School of Data Science, Department of Data Science,
 * Yokohama City University
 *
 * Contact:
 * Address: 22-2, Seto, Kanazawa-ku, Yokohama, Kanagawa 236-0027 Japan
 * E-mail: ahsato@yokohama-cu.ac.jp
 * TEL: +81-45-787-2208
 *
 * Three types of functions are defined in this library.
 * 1. calculate representative geographical position(s) (latitude, longitude)
 * of a grid square from a grid square code
 * 2. calculate a grid square code from a geographical position (latitude, longitude)
 * 3. calculate geodesic distance and size of grid square (representative lengths and area)
 *
 * 1.
 *
 * meshcode_to_latlong(meshcode, extension=F, spec="")
 * : calculate northen western geographic position
 * of the grid (latitude, longitude) from meshcode
 * meshcode_to_latlong_NW(meshcode, extension=F, spec="")
 * : calculate northen western geographic position
 * of the grid (latitude, longitude) from meshcode
 * meshcode_to_latlong_SW(meshcode, extension=F, spec="")
 * : calculate sourthern western geographic position
 * of the grid (latitude, longitude) from meshcode
 * meshcode_to_latlong_NE(meshcode, extension=F, spec="")
 * : calculate northern eastern geographic position
 * of the grid (latitude, longitude) from meshcode
 * meshcode_to_latlong_SE(meshcode, extension=F, spec="")
 * : calculate sourthern eastern geographic position
 * of the grid (latitude, longitude) from meshcode
 * meshcode_to_latlong_grid(meshcode, extension=F, spec="")
 * : calculate northern western and sourthern eastern geographic positions
 * of the grid (latitude0, longitude0, latitude1, longitude1) from meshcode
 * The extension option determines high resolution gird square codes.
 * If extension=T, we can also specify several levels of grid squares
 * by using the spec option. If the spec option is ignored, then the default
 * level of grid square is selectedin accordance with the number
 * of digits of grid squre code:
 * 100m x 100m 3 arc-second x 4.5 arc second (12 digits) : default
 * 100m x 100m 3 arc-second x 4.5 arc second (13 digits) : spec = "ex100m_13"
 * 50m x 50m 1.5 arc-second x 2.25 arc second (13 digits) : default
 * 50m x 50m 1.5 arc-second x 2.25 arc second (14 digits) : sepc = "ex50m_14"
 * 10m x 10m 0.3 arc-second x 0.45 arc second (14 digits) : default
 * 1m x 1m 0.03 arc-second x 0.045 arc second (16 digits) : default
 *
 * 2.
 *
 * : calculate a basic (1km) grid square code (10 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode1(latitude,longitude)
 * : calculate an 80km grid square code (6 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode2(latitude,longitude)
 * : calculate a 10km grid square code (8 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode3(latitude,longitude)
 * : calculate a 1km grid square code (10 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode4(latitude,longitude)
 * : calculate a 500m grid square code (11 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode5(latitude,longitude)
 * : calculate a 250m grid square code (12 digits)
 * from a geographical position (latitude, longitude)
 * cal_meshcode6(latitude,longitude)
 * : calculate a 125m grid square code (13 digits)
 * from a geographical position (latitude, longitude)
 *
 * This grid square code set is not included in JIS X0410 directly but useful.
 *
 * cal_meshcode_ex100m_12(latitude,longitude)
 * : calculate an extended 100m (1km / 10) grid square code (12 digits)
 * from a geographical position (latitude, longitude) - 3 arc-second
 * for latitude and 4.5 arc-second for longitude
 *
 * cal_meshcode_ex100m_13(latitude,longitude)
 * : calculate an extended 100m (500m / 5) grid square code (13 digits)
 * from a geographical position (latitude, longitude) - 3 arc-second
 * for latitude and 4.5 arc-second for longitude 
 *
 * cal_meshcode_ex50m_13(latitude,longitude)
 * : calculate an extended 50m (100m / 2) grid square code (13 digits)
 * from a geographical position (latitude, longitude) - 1.5 arc-second
 * for latitude and 2.25 arc-second for longitude
 *
 * cal_meshcode_ex50m_14(latitude,longitude)
 * : calculate an extended 50m (100m / 2) grid square code (14 digits)
 * from a geographical position (latitude, longitude) - 1.5 arc-second
 * for latitude and 2.25 arc-second for longitude 
 *
 * cal_meshcode_ex10m_14(latitude,longitude)
 * : calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits)
 * from a geographical position (latitude, longitude) - 0.3 arc-second
 * for latitude and 0.45 arc-second for longitude
 *
 * cal_meshcode_ex1m_16(latitude,longitude)
 * : calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits)
 * from a geographical position (latitude, longitude) - 0.03 arc-second for
 * latitude and 0.045 arc-second for longitude
 *
 * Structure of the world grid square code with compatibility to JIS X 0410
 * A : area code (1 digit) A takes 1 to 8
 * ABBBBB : 80km grid square code (40 arc-minutes for latitude, 1 arc-degree for longitude) (6 digits)
 * ABBBBBCC : 10km grid square code (5 arc-minutes for latitude, 7.5 arc-minutes for longitude) (8 digits)
 * ABBBBBCCDD : 1km grid square code (30 arc-seconds for latitude, 45 arc-secondes for longitude) (10 digits)
 * ABBBBBCCDDE : 500m grid square code (15 arc-seconds for latitude, 22.5 arc-seconds for longitude) (11 digits)
 * ABBBBBCCDDEF : 250m grid square code (7.5 arc-seconds for latitude, 11.25 arc-seconds for longitude) (12 digits)
 * ABBBBBCCDDEFG : 125m grid square code (3.75 arc-seconds for latitude, 5.625 arc-seconds for longitude) (13 digits)
 * ABBBBBCCDDHH : Extended 100m grid square code with 12 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (12 digits)
 * ABBBBBCCDDEHH : Extended 100m grid square code with 13 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (13 digits)
 * ABBBBBCCDDHHE : Extended 50m grid square code with 13 digits (1.5 arc-seconds for latitude, 2.25 arc-seconds for longitude) (13 digits)
 * ABBBBBCCDDEHHI : Extended 50m grid square code with 14 digits (1.5 arc-seconds for latitude, 2.25 arc-seconds for longitude) (14 digits)
 * ABBBBBCCDDHHGG : Extended 10m grid square code with 14 digits (0.3 arc-seconds for latitude, 0.45 arc-seconds for longitude) (14 digits)
 * ABBBBBCCDDHHGGII : Extended 1m grid square code with 16 digits (0.03 arc-seconds for latitude, 0.045 arc-seconds for longitude) (16 digits)
 *
 * 3.
 * Calculate geodesic distance and size of world grid square 
 *
 * T. Vincenty, ``Direct and Inverse Solutions of Geodesics
 * on the Ellipsoid with application of nested equations'',
 * Survey Review XXIII, Vol 176 (1975) Vol. 88-93.
 *
 * Vincenty(latitude1, longitude1, latitude2, longitude)
 * : calculate geodesitc distance between two points (latitude1, longitude1)
 * and (latitude2, longitude2) placed on the WGS84 Earth ellipsoid
 * based on the Vincenty's formulae (1975)
 * cal_area_from_meshcode(meshcode,extension=F,spec="")
 * : calculate size (northern west-to-east span H1, sothern west-to-east span H2,
 * north-to-south span W, and area approximated by trapezoide A)
 * of world grid square indicated by meshcode
 * cal_area_from_latlong(latlong)
 * : calculate size (northern west-to-east span H1, sothern west-to-east span H2,
 * north-to-south span W, and area approximated by trapezoid A)
 * of a trapezoid on the WGS84 Earth ellipoid indicated
 * by (latlong$lat0, latlong$long0, latlong$lat1, latlong$long1)
 *
 */

#include "worldmesh.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* ============================================================================
 * Constants and Macros
 * ============================================================================ */

#define PI 3.14159265358979323846
#define DEG_TO_RAD(x) ((x) * PI / 180.0)
#define RAD_TO_DEG(x) ((x) * 180.0 / PI)

/* WGS84 Ellipsoid Parameters */
#define WGS84_A      6378137.0          /* Semi-major axis (meters) */
#define WGS84_B      6356752.314245     /* Semi-minor axis (meters) */
#define WGS84_F      (1.0/298.257223563) /* Flattening */

#define VINCENTY_TOLERANCE 1e-12

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Validate latitude and longitude bounds.
 */
int worldmesh_validate_coordinates(double latitude, double longitude) {
    return (latitude >= -90.0 && latitude <= 90.0 &&
            longitude >= -180.0 && longitude <= 180.0);
}

/**
 * Extract a substring from a string (0-indexed).
 * Similar to R's substring(str, start, end) but with 0-based indexing.
 * Returns NULL-terminated substring copied to dest.
 */
static char* extract_substring(const char *src, int start, int end, char *dest) {
    if (!src || !dest || start < 0 || end < start) return NULL;
    int len = end - start + 1;
    strncpy(dest, src + start, len);
    dest[len] = '\0';
    return dest;
}

/**
 * Safe string to long conversion.
 * Returns 0 on error, actual value on success.
 */
static long safe_strtol(const char *str) {
    char *endptr;
    long val = strtol(str, &endptr, 10);
    if (endptr == str) return 0; /* No conversion performed */
    return val;
}

/**
 * Format meshcode with zero-padding for code12 and code34.
 * R uses paste() which handles leading zeros; we replicate that behavior.
 */
static char* format_meshcode(int o, int code12, int code34, 
                              const char *rest, char *out) {
    if (code12 < 10) {
        if (code34 < 10) {
            snprintf(out, 256, "%d00%d0%d%s", o, code12, code34, rest);
        } else {
            snprintf(out, 256, "%d00%d%d%s", o, code12, code34, rest);
        }
    } else if (code12 < 100) {
        if (code34 < 10) {
            snprintf(out, 256, "%d0%d0%d%s", o, code12, code34, rest);
        } else {
            snprintf(out, 256, "%d0%d%d%s", o, code12, code34, rest);
        }
    } else {
        if (code34 < 10) {
            snprintf(out, 256, "%d%d0%d%s", o, code12, code34, rest);
        } else {
            snprintf(out, 256, "%d%d%d%s", o, code12, code34, rest);
        }
    }
    return out;
}

/* ============================================================================
 * Core Grid Calculation Function
 * ============================================================================ */

/**
 * meshcode_to_latlong_grid: The core function that interprets all meshcode formats.
 * Handles 6, 8, 10, 11, 12, 13, 14, and 16 digit codes.
 * With extension=1, supports extended high-resolution codes.
 */
WorldMesh_Grid* meshcode_to_latlong_grid(const char *meshcode, int extension,
                                          const char *spec, WorldMesh_Grid *out) {
    if (!meshcode || !out) return NULL;
    
    int len = strlen(meshcode);
    if (len < 6 || (len != 6 && len != 8 && len != 10 && len != 11 && 
                     len != 12 && len != 13 && len != 14 && len != 16)) {
        return NULL;
    }
    
    /* Extract area code (code0): 1-8 in R, becomes 0-7 internally */
    char buf[3];
    extract_substring(meshcode, 0, 0, buf);
    long code0_raw = safe_strtol(buf);
    if (code0_raw < 1 || code0_raw > 8) return NULL;
    long code0 = code0_raw - 1;  /* Transform 1-8 to 0-7 */
    
    /* Calculate x, y, z from area code */
    long z = code0 % 2;
    long y = ((code0 - z) / 2) % 2;
    long x = (code0 - 2*y - z) / 4;
    
    /* Extract code12 (2nd-4th digits, positions 1-3) */
    char code12_str[4] = {meshcode[1], meshcode[2], meshcode[3], '\0'};
    long code12 = safe_strtol(code12_str);
    
    /* Extract code34 (5th-6th digits, positions 4-5) */
    char code34_str[3] = {meshcode[4], meshcode[5], '\0'};
    long code34 = safe_strtol(code34_str);
    
    /* Variables for higher precision levels */
    long code5 = 0, code6 = 0;
    long code7 = 0, code8 = 0;
    long code9 = 0, code10 = 0, code11 = 0;
    long codeex9 = 0, codeex10 = 0, codeex11 = 0, codeex12 = 0;
    long codeex13 = 0, codeex14 = 0;
    
    if (len >= 8) {
        code5 = meshcode[6] - '0';  /* 7th position (0-indexed) */
        code6 = meshcode[7] - '0';  /* 8th position (0-indexed) */
    }
    
    if (len >= 10) {
        code7 = meshcode[8] - '0';   /* 9th position */
        code8 = meshcode[9] - '0';   /* 10th position */
    }
    
    if (!extension) {
        if (len >= 11) {
            code9 = meshcode[10] - '0'; /* 11th position */
        }
        if (len >= 12) {
            code10 = meshcode[11] - '0'; /* 12th position */
        }
        if (len >= 13) {
            code11 = meshcode[12] - '0'; /* 13th position */
        }
    } else {
        if (len >= 12) {
            codeex9 = meshcode[10] - '0';  /* 11th position */
            codeex10 = meshcode[11] - '0'; /* 12th position */
        }
        if (len >= 13) {
            code9 = meshcode[10] - '0';   /* 11th position */
            codeex10 = meshcode[11] - '0'; /* 12th position */
            codeex11 = meshcode[12] - '0'; /* 13th position */
        }
        if (len >= 14) {
            codeex11 = meshcode[12] - '0'; /* 13th position */
            codeex12 = meshcode[13] - '0'; /* 14th position */
        }
        if (len == 16) {
            codeex13 = meshcode[14] - '0'; /* 15th position */
            codeex14 = meshcode[15] - '0'; /* 16th position */
        }
    }
    
    double lat0, long0, dlat, dlong;
    
    /* 1st grid: 6 digits */
    if (len == 6) {
        lat0 = (code12 - x + 1) * 2.0 / 3.0;
        long0 = code34 + 100.0 * z;
        lat0 = (1.0 - 2.0*x) * lat0;
        long0 = (1.0 - 2.0*y) * long0;
        dlat = 2.0 / 3.0;
        dlong = 1.0;
    }
    /* 2nd grid: 8 digits */
    else if (len == 8) {
        lat0 = code12 * 2.0 / 3.0;
        long0 = code34 + 100.0 * z;
        lat0 = lat0 + ((code5 - x + 1) * 2.0 / 3.0) / 8.0;
        long0 = long0 + (code6 + y) / 8.0;
        lat0 = (1.0 - 2.0*x) * lat0;
        long0 = (1.0 - 2.0*y) * long0;
        dlat = 2.0 / 3.0 / 8.0;
        dlong = 1.0 / 8.0;
    }
    /* 3rd grid: 10 digits */
    else if (len == 10) {
        lat0 = code12 * 2.0 / 3.0;
        long0 = code34 + 100.0 * z;
        lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
        long0 = long0 + code6 / 8.0;
        lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
        long0 = long0 + (code8 + y) / 8.0 / 10.0;
        lat0 = (1.0 - 2.0*x) * lat0;
        long0 = (1.0 - 2.0*y) * long0;
        dlat = 2.0 / 3.0 / 8.0 / 10.0;
        dlong = 1.0 / 8.0 / 10.0;
    }
    /* 4th grid: 11 digits */
    else if (len == 11) {
        lat0 = code12 * 2.0 / 3.0;
        long0 = code34 + 100.0 * z;
        lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
        long0 = long0 + code6 / 8.0;
        lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
        long0 = long0 + (code8 + y) / 8.0 / 10.0;
        lat0 = lat0 + (floor((code9 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
        long0 = long0 + ((code9 - 1) % 2 - y) / 8.0 / 10.0 / 2.0;
        lat0 = (1.0 - 2.0*x) * lat0;
        long0 = (1.0 - 2.0*y) * long0;
        dlat = 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
        dlong = 1.0 / 8.0 / 10.0 / 2.0;
    }
    /* 5th grid: 12 digits */
    else if (len == 12) {
        if (!extension) {
            /* Standard 250m code */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + (code8 + y) / 8.0 / 10.0;
            lat0 = lat0 + (floor((code9 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
            long0 = long0 + ((code9 - 1) % 2 - y) / 8.0 / 10.0 / 2.0;
            lat0 = lat0 + (floor((code10 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0;
            long0 = long0 + ((code10 - 1) % 2 - y) / 8.0 / 10.0 / 2.0 / 2.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0;
            dlong = 1.0 / 8.0 / 10.0 / 2.0 / 2.0;
        } else {
            /* Extended 100m code */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + (code7 * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + code8 / 8.0 / 10.0;
            lat0 = lat0 + ((codeex9 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0;
            long0 = long0 + (codeex10 + y) / 8.0 / 10.0 / 10.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 10.0;
            dlong = 1.0 / 8.0 / 10.0 / 10.0;
        }
    }
    /* 6th grid: 13 digits */
    else if (len == 13) {
        if (!extension) {
            /* Standard 125m code */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + (code8 + y) / 8.0 / 10.0;
            lat0 = lat0 + (floor((code9 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
            long0 = long0 + ((code9 - 1) % 2 - y) / 8.0 / 10.0 / 2.0;
            lat0 = lat0 + (floor((code10 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0;
            long0 = long0 + ((code10 - 1) % 2 - y) / 8.0 / 10.0 / 2.0 / 2.0;
            lat0 = lat0 + (floor((code11 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0 / 2.0;
            long0 = long0 + ((code11 - 1) % 2 - y) / 8.0 / 10.0 / 2.0 / 2.0 / 2.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0 / 2.0;
            dlong = 1.0 / 8.0 / 10.0 / 2.0 / 2.0 / 2.0;
        } else if (spec && strcmp(spec, "ex100m_13") == 0) {
            /* Extended 100m code (13 digits) */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + (code8 + y) / 8.0 / 10.0;
            lat0 = lat0 + (floor((code9 - 1) / 2.0) + 2.0*x - 2.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
            long0 = long0 + ((code9 - 1) % 2 - 2.0*y) / 8.0 / 10.0 / 2.0;
            lat0 = lat0 + ((codeex10 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 2.0 / 5.0;
            long0 = long0 + (codeex11 + y) / 8.0 / 10.0 / 2.0 / 5.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 5.0;
            dlong = 1.0 / 8.0 / 10.0 / 2.0 / 5.0;
        } else {
            /* Extended 50m code (13 digits) */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + (code7 * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + code8 / 8.0 / 10.0;
            lat0 = lat0 + ((codeex9 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0;
            long0 = long0 + (codeex10 + y) / 8.0 / 10.0 / 10.0;
            lat0 = lat0 + (floor((codeex11 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 10.0 / 2.0;
            long0 = long0 + ((codeex11 - 1) % 2 - y) / 8.0 / 10.0 / 10.0 / 2.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 10.0 / 2.0;
            dlong = 1.0 / 8.0 / 10.0 / 10.0 / 2.0;
        }
    }
    /* 14 digits */
    else if (len == 14) {
        if (spec && strcmp(spec, "ex50m_14") == 0) {
            /* Extended 50m code (14 digits) */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + ((code7 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + (code8 + y) / 8.0 / 10.0;
            lat0 = lat0 + (floor((code9 - 1) / 2.0) + 2.0*x - 2.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0;
            long0 = long0 + ((code9 - 1) % 2 - 2.0*y) / 8.0 / 10.0 / 2.0;
            lat0 = lat0 + ((codeex10 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 2.0 / 5.0;
            long0 = long0 + (codeex11 + y) / 8.0 / 10.0 / 2.0 / 5.0;
            lat0 = lat0 + (floor((codeex12 - 1) / 2.0) + x - 1.0) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 5.0 / 2.0;
            long0 = long0 + ((codeex12 - 1) % 2 - y) / 8.0 / 10.0 / 2.0 / 5.0 / 2.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 5.0 / 2.0;
            dlong = 1.0 / 8.0 / 10.0 / 2.0 / 5.0 / 2.0;
        } else {
            /* Extended 10m code (14 digits) */
            lat0 = code12 * 2.0 / 3.0;
            long0 = code34 + 100.0 * z;
            lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
            long0 = long0 + code6 / 8.0;
            lat0 = lat0 + (code7 * 2.0 / 3.0) / 8.0 / 10.0;
            long0 = long0 + code8 / 8.0 / 10.0;
            lat0 = lat0 + (codeex9 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0;
            long0 = long0 + codeex10 / 8.0 / 10.0 / 10.0;
            lat0 = lat0 + ((codeex11 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0;
            long0 = long0 + (codeex12 + y) / 8.0 / 10.0 / 10.0 / 10.0;
            lat0 = (1.0 - 2.0*x) * lat0;
            long0 = (1.0 - 2.0*y) * long0;
            dlat = 2.0 / 3.0 / 8.0 / 10.0 / 10.0 / 10.0;
            dlong = 1.0 / 8.0 / 10.0 / 10.0 / 10.0;
        }
    }
    /* 16 digits */
    else if (len == 16) {
        /* Extended 1m code (16 digits) */
        lat0 = code12 * 2.0 / 3.0;
        long0 = code34 + 100.0 * z;
        lat0 = lat0 + (code5 * 2.0 / 3.0) / 8.0;
        long0 = long0 + code6 / 8.0;
        lat0 = lat0 + (code7 * 2.0 / 3.0) / 8.0 / 10.0;
        long0 = long0 + code8 / 8.0 / 10.0;
        lat0 = lat0 + (codeex9 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0;
        long0 = long0 + codeex10 / 8.0 / 10.0 / 10.0;
        lat0 = lat0 + (codeex11 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0;
        long0 = long0 + codeex12 / 8.0 / 10.0 / 10.0 / 10.0;
        lat0 = lat0 + ((codeex13 - x + 1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0 / 10.0;
        long0 = long0 + (codeex14 + y) / 8.0 / 10.0 / 10.0 / 10.0 / 10.0;
        lat0 = (1.0 - 2.0*x) * lat0;
        long0 = (1.0 - 2.0*y) * long0;
        dlat = 2.0 / 3.0 / 8.0 / 10.0 / 10.0 / 10.0 / 10.0;
        dlong = 1.0 / 8.0 / 10.0 / 10.0 / 10.0 / 10.0;
    }
    else {
        return NULL;
    }
    
    out->lat0 = lat0;
    out->long0 = long0;
    out->lat1 = lat0 - dlat;
    out->long1 = long0 + dlong;
    
    return out;
}

/**
 * meshcode_to_latlong_NW: Northwest (top-left) corner
 */
WorldMesh_Grid* meshcode_to_latlong_NW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out) {
    WorldMesh_Grid temp;
    WorldMesh_Grid *res = meshcode_to_latlong_grid(meshcode, extension, spec, &temp);
    if (!res) return NULL;
    out->lat0 = res->lat0;
    out->long0 = res->long0;
    out->lat1 = 0;
    out->long1 = 0;
    return out;
}

/**
 * meshcode_to_latlong_SW: Southwest (bottom-left) corner
 */
WorldMesh_Grid* meshcode_to_latlong_SW(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out) {
    WorldMesh_Grid temp;
    WorldMesh_Grid *res = meshcode_to_latlong_grid(meshcode, extension, spec, &temp);
    if (!res) return NULL;
    out->lat0 = res->lat1;
    out->long0 = res->long0;
    out->lat1 = 0;
    out->long1 = 0;
    return out;
}

/**
 * meshcode_to_latlong_NE: Northeast (top-right) corner
 */
WorldMesh_Grid* meshcode_to_latlong_NE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out) {
    WorldMesh_Grid temp;
    WorldMesh_Grid *res = meshcode_to_latlong_grid(meshcode, extension, spec, &temp);
    if (!res) return NULL;
    out->lat0 = res->lat0;
    out->long0 = res->long1;
    out->lat1 = 0;
    out->long1 = 0;
    return out;
}

/**
 * meshcode_to_latlong_SE: Southeast (bottom-right) corner
 */
WorldMesh_Grid* meshcode_to_latlong_SE(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Grid *out) {
    WorldMesh_Grid temp;
    WorldMesh_Grid *res = meshcode_to_latlong_grid(meshcode, extension, spec, &temp);
    if (!res) return NULL;
    out->lat0 = res->lat1;
    out->long0 = res->long1;
    out->lat1 = 0;
    out->long1 = 0;
    return out;
}

/* ============================================================================
 * Coordinates to Meshcode Functions
 * ============================================================================ */

/**
 * Helper function to determine area code from coordinates.
 * o: Output area code (1-8)
 * z, y, x: Internal parameters
 */
static void get_area_code(double latitude, double longitude, 
                          int *o, long *z, long *y, long *x) {
    *o = 0;
    if (latitude < 0.0) *o = 4;
    if (longitude < 0.0) *o += 2;
    if (fabs(longitude) >= 100.0) *o += 1;
    
    *z = *o % 2;
    *y = ((*o - *z) / 2) % 2;
    *x = (*o - 2 * (*y) - *z) / 4;
    
    *o = *o + 1;
}

/**
 * cal_meshcode1: 80km grid square code (6 digits)
 */
char* cal_meshcode1(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    long u = (long)floor(longitude - 100.0 * z);
    
    format_meshcode(o, (int)p, (int)u, "", out);
    return out;
}

/**
 * cal_meshcode2: 10km grid square code (8 digits)
 */
char* cal_meshcode2(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "99999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 / 40.0 - p;
    a = a * 40.0;
    long q = (long)floor(a / 5.0);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    
    char rest[10];
    sprintf(rest, "%ld%ld", q, v);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode3: 1km grid square code (10 digits)
 */
char* cal_meshcode3(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "9999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 / 40.0 - p;
    a = a * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a / 5.0 - q;
    b = b * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 / 7.5 - v;
    g = g * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld", q, v, r, w);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode4: 500m grid square code (11 digits)
 */
char* cal_meshcode4(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "99999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 / 40.0 - p;
    a = a * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a / 5.0 - q;
    b = b * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s2u = (long)floor(c / 15.0);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 / 7.5 - v;
    g = g * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long s2l = (long)floor(h / 22.5);
    long s2 = s2u * 2 + s2l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld", q, v, r, w, s2);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode5: 250m grid square code (12 digits)
 */
char* cal_meshcode5(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 / 40.0 - p;
    a = a * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a / 5.0 - q;
    b = b * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s2u = (long)floor(c / 15.0);
    double d = c - s2u * 15.0;
    long s4u = (long)floor(d / 7.5);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 / 7.5 - v;
    g = g * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long s2l = (long)floor(h / 22.5);
    double i = h - s2l * 22.5;
    long s4l = (long)floor(i / 11.25);
    long s2 = s2u * 2 + s2l + 1;
    long s4 = s4u * 2 + s4l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld", q, v, r, w, s2, s4);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode6: 125m grid square code (13 digits)
 */
char* cal_meshcode6(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "9999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s2u = (long)floor(c / 15.0);
    double d = c - s2u * 15.0;
    long s4u = (long)floor(d / 7.5);
    double e = d - s4u * 7.5;
    long s8u = (long)floor(e / 3.75);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long s2l = (long)floor(h / 22.5);
    double i = h - s2l * 22.5;
    long s4l = (long)floor(i / 11.25);
    double j = i - s4l * 11.25;
    long s8l = (long)floor(j / 5.625);
    long s2 = s2u * 2 + s2l + 1;
    long s4 = s4u * 2 + s4l + 1;
    long s8 = s8u * 2 + s8l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s2, s4, s8);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex100m_12: Extended 100m code (12 digits)
 */
char* cal_meshcode_ex100m_12(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s = (long)floor(c / 3.0);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long xx = (long)floor(h / 4.5);
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld", q, v, r, w, s, xx);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex100m_13: Extended 100m code (13 digits)
 */
char* cal_meshcode_ex100m_13(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "9999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s2u = (long)floor(c / 15.0);
    double d = c - s2u * 15.0;
    long et = (long)floor(d / 3.0);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long s2l = (long)floor(h / 22.5);
    double i = h - s2l * 22.5;
    long jt = (long)floor(i / 4.5);
    long s2 = s2u * 2 + s2l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s2, et, jt);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex50m_13: Extended 50m code (13 digits)
 */
char* cal_meshcode_ex50m_13(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "9999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s = (long)floor(c / 3.0);
    double d = c - s * 3.0;
    long t2u = (long)floor(d / 1.5);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long xx = (long)floor(h / 4.5);
    double i = h - xx * 4.5;
    long t2l = (long)floor(i / 2.25);
    long tt = t2u * 2 + t2l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s, xx, tt);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex50m_14: Extended 50m code (14 digits)
 */
char* cal_meshcode_ex50m_14(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "99999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s2u = (long)floor(c / 15.0);
    double d = c - s2u * 15.0;
    long et = (long)floor(d / 3.0);
    double e = d - et * 3.0;
    long t2u = (long)floor(e / 1.5);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long s2l = (long)floor(h / 22.5);
    double i = h - s2l * 22.5;
    long jt = (long)floor(i / 4.5);
    double j = i - jt * 4.5;
    long t2l = (long)floor(j / 2.25);
    long s2 = s2u * 2 + s2l + 1;
    long t2 = t2u * 2 + t2l + 1;
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s2, et, jt, t2);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex10m_14: Extended 10m code (14 digits)
 */
char* cal_meshcode_ex10m_14(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "99999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s = (long)floor(c / 3.0);
    double d = c - s * 3.0;
    long t = (long)floor(d / 0.3);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long xx = (long)floor(h / 4.5);
    double i = h - xx * 4.5;
    long yy = (long)floor(i / 0.45);
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s, xx, t, yy);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/**
 * cal_meshcode_ex1m_16: Extended 1m code (16 digits)
 */
char* cal_meshcode_ex1m_16(double latitude, double longitude, char *out) {
    if (!out) return NULL;
    if (!worldmesh_validate_coordinates(latitude, longitude)) {
        strcpy(out, "9999999999999999");
        return out;
    }
    
    int o;
    long z, y, x;
    get_area_code(latitude, longitude, &o, &z, &y, &x);
    
    latitude = (1.0 - 2.0*x) * latitude;
    longitude = (1.0 - 2.0*y) * longitude;
    
    long p = (long)floor(latitude * 60.0 / 40.0);
    double a = latitude * 60.0 - p * 40.0;
    long q = (long)floor(a / 5.0);
    double b = a - q * 5.0;
    long r = (long)floor(b * 60.0 / 30.0);
    double c = b * 60.0 - r * 30.0;
    long s = (long)floor(c / 3.0);
    double d = c - s * 3.0;
    long t = (long)floor(d / 0.3);
    double e = d - t * 0.3;
    long tt = (long)floor(e / 0.03);
    long u = (long)floor(longitude - 100.0 * z);
    double f = longitude - 100.0 * z - u;
    long v = (long)floor(f * 60.0 / 7.5);
    double g = f * 60.0 - v * 7.5;
    long w = (long)floor(g * 60.0 / 45.0);
    double h = g * 60.0 - w * 45.0;
    long xx = (long)floor(h / 4.5);
    double i = h - xx * 4.5;
    long yy = (long)floor(i / 0.45);
    double j = i - yy * 0.45;
    long zz = (long)floor(j / 0.045);
    
    char rest[256];
    sprintf(rest, "%ld%ld%ld%ld%ld%ld%ld%ld%ld%ld", q, v, r, w, s, xx, t, yy, tt, zz);
    format_meshcode(o, (int)p, (int)u, rest, out);
    return out;
}

/* ============================================================================
 * Vincenty's Formulae and Area Functions
 * ============================================================================ */

/**
 * Vincenty: Calculate geodesic distance between two points on WGS84 ellipsoid.
 */
double Vincenty(double latitude1, double longitude1,
                double latitude2, double longitude2) {
    if (fabs(latitude1 - latitude2) < 1e-18 && 
        fabs(longitude1 - longitude2) < 1e-18) {
        return 0.0;
    }
    
    double L = (longitude1 - longitude2) / 180.0 * PI;
    double U1 = atan((1.0 - WGS84_F) * tan(DEG_TO_RAD(latitude1)));
    double U2 = atan((1.0 - WGS84_F) * tan(DEG_TO_RAD(latitude2)));
    double lambda = L;
    double dlambda = 10.0;
    
    double cs, cscc, sinsigma, cossigma, sigma, sinalpha, cos2alpha;
    double cos2sigmam, C, lambda0;
    
    while (fabs(dlambda) > VINCENTY_TOLERANCE) {
        cs = cos(U2) * sin(lambda);
        cscc = cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lambda);
        sinsigma = sqrt(cs * cs + cscc * cscc);
        cossigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lambda);
        sigma = atan(sinsigma / cossigma);
        sinalpha = cos(U1) * cos(U2) * sin(lambda) / sinsigma;
        cos2alpha = 1.0 - sinalpha * sinalpha;
        
        if (cos2alpha == 0.0) {
            C = 0.0;
            lambda0 = L + WGS84_F * sinalpha * sigma;
        } else {
            cos2sigmam = cossigma - 2.0 * sin(U1) * sin(U2) / cos2alpha;
            C = WGS84_F / 16.0 * cos2alpha * (4.0 + WGS84_F * (4.0 - 3.0 * cos2alpha));
            lambda0 = L + (1.0 - C) * WGS84_F * sinalpha * 
                     (sigma + C * sinsigma * (cos2sigmam + C * cossigma * 
                     (-1.0 + 2.0 * cos2sigmam * cos2sigmam)));
        }
        
        dlambda = lambda0 - lambda;
        lambda = lambda0;
    }
    
    double A, dsigma, u2, B;
    if (C == 0.0) {
        A = 1.0;
        dsigma = 0.0;
    } else {
        u2 = cos2alpha * (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / 
             (WGS84_B * WGS84_B);
        A = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
        B = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
        dsigma = B * sinsigma * (cos2sigmam + 0.25 * B * 
                (cossigma * (-1.0 + 2.0 * cos2sigmam * cos2sigmam) - 
                 1.0 / 6.0 * B * cos2sigmam * (-3.0 + 4.0 * sinsigma * sinsigma) * 
                 (-3.0 + 4.0 * cos2sigmam * cos2sigmam)));
    }
    
    double s = WGS84_B * A * (sigma - dsigma);
    return s;
}

/**
 * cal_area_from_latlong: Calculate size and area from grid corners.
 */
WorldMesh_Area* cal_area_from_latlong(const WorldMesh_Grid *grid,
                                       WorldMesh_Area *out) {
    if (!grid || !out) return NULL;
    
    double W1 = Vincenty(grid->lat0, grid->long0, grid->lat0, grid->long1);
    double W2 = Vincenty(grid->lat1, grid->long0, grid->lat1, grid->long1);
    double H = Vincenty(grid->lat0, grid->long0, grid->lat1, grid->long0);
    double A = (W1 + W2) * H * 0.5;
    
    out->W1 = W1;
    out->W2 = W2;
    out->H = H;
    out->A = A;
    
    return out;
}

/**
 * cal_area_from_meshcode: Calculate area from meshcode.
 */
WorldMesh_Area* cal_area_from_meshcode(const char *meshcode, int extension,
                                        const char *spec, WorldMesh_Area *out) {
    if (!meshcode || !out) return NULL;
    
    WorldMesh_Grid grid;
    WorldMesh_Grid *res = meshcode_to_latlong_grid(meshcode, extension, spec, &grid);
    if (!res) return NULL;
    
    return cal_area_from_latlong(&grid, out);
}

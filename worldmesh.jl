#
# Copyright (c) 2015-2026 Research Institute for World Grid Squares 
# Prof. Dr. Aki-Hiro Sato (Yokohama City University)
# All rights reserved. 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Julia functions to calculate the world grid square code.
# The world grid square code computed by this library 
# is compatible to Japanese grid square code (JIS X 0410). 
#
# This library for World grid square code complies ISO 24108-1:
# https://www.iso.org/standard/87706.html
#
# Ported from R to Julia by Cline (03 May 2026)
# Original R version:
# Version 1.81 : Released on 12 December 2022
#
# Written by Prof. Dr. Aki-Hiro Sato
# Graduate School of Data Science, Department of Data Science,
# Yokohama City University
#
# Contact:
# Address: 22-2, Seto, Kanazawa-ku, Yokohama, Kanagawa 236-0027 Japan
# E-mail: ahsato@yokohama-cu.ac.jp
# TEL: +81-45-787-2208
#
# Three types of functions are defined in this library.
# 1. calculate representative geographical position(s) (latitude, longitude) of a grid square from a grid square code
# 2. calculate a grid square code from a geographical position (latitude, longitude)
# 3. calculate geodesic distance and size of grid square (representative lengths and area)
#
# 1.
#
# meshcode_to_latlong(meshcode, extension=F, spec="")
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NW(meshcode, extension=F, spec="")
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SW(meshcode, extension=F, spec="")
# : calculate sourthern western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NE(meshcode, extension=F, spec="")
# : calculate northern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SE(meshcode, extension=F, spec="")
# : calculate sourthern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_grid(meshcode, extension=F, spec="")
# : calculate northern western and sourthern eastern geographic positions of the grid (latitude0, longitude0, latitude1, longitude1) from meshcode
# The extension option determines high resolution gird square codes. If extension=T, we can also specify several levels 
# of grid squares by using the spec option. If the spec option is ignored, then the default level of grid square is selected
# in accordance with the number of digits of grid squre code:
# 100m x 100m 3 arc-second x 4.5 arc second (12 digits) : default
# 100m x 100m 3 arc-second x 4.5 arc second (13 digits) : spec = "ex100m_13"
# 50m x 50m 1.5 arc-second x 2.25 arc second (13 digits) : default
# 50m x 50m 1.5 arc-second x 2.25 arc second (14 digits) : sepc = "ex50m_14"
# 10m x 10m 0.3 arc-second x 0.45 arc second (14 digits) : default
# 1m x 1m 0.03 arc-second x 0.045 arc second (16 digits) : default
#
# 2.
#
# : calculate a basic (1km) grid square code (10 digits) from a geographical position (latitude, longitude)
# cal_meshcode1(latitude,longitude)
# : calculate an 80km grid square code (6 digits) from a geographical position (latitude, longitude)
# cal_meshcode2(latitude,longitude)
# : calculate a 10km grid square code (8 digits) from a geographical position (latitude, longitude)
# cal_meshcode3(latitude,longitude)
# : calculate a 1km grid square code (10 digits) from a geographical position (latitude, longitude)
# cal_meshcode4(latitude,longitude)
# : calculate a 500m grid square code (11 digits) from a geographical position (latitude, longitude)
# cal_meshcode5(latitude,longitude)
# : calculate a 250m grid square code (12 digits) from a geographical position (latitude, longitude)
# cal_meshcode6(latitude,longitude)
# : calculate a 125m grid square code (13 digits) from a geographical position (latitude, longitude)
#
# This grid square code set is not included in JIS X0410 directly but useful.
#
# cal_meshcode_ex100m_12(latitude,longitude)
# : calculate an extended 100m (1km / 10) grid square code (12 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude
#
# cal_meshcode_ex100m_13(latitude,longitude)
# : calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude 
#
# cal_meshcode_ex50m_13(latitude,longitude)
# : calculate an extended 50m (100m / 2) grid square code (13 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude
#
# cal_meshcode_ex50m_14(latitude,longitude)
# : calculate an extended 50m (100m / 2) grid square code (14 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude 
#
# cal_meshcode_ex10m_14(latitude,longitude)
# : calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) from a geographical position (latitude, longitude) - 0.3 arc-second for latitude and 0.45 arc-second for longitude
#
# cal_meshcode_ex1m_16(latitude,longitude)
# : calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) from a geographical position (latitude, longitude) - 0.03 arc-second for latitude and 0.045 arc-second for longitude
#
# Structure of the world grid square code with compatibility to JIS X0410
# A : area code (1 digit) A takes 1 to 8
# ABBBBB : 80km grid square code (40 arc-minutes for latitude, 1 arc-degree for longitude) (6 digits)
# ABBBBBCC : 10km grid square code (5 arc-minutes for latitude, 7.5 arc-minutes for longitude) (8 digits)
# ABBBBBCCDD : 1km grid square code (30 arc-seconds for latitude, 45 arc-secondes for longitude) (10 digits)
# ABBBBBCCDDE : 500m grid square code (15 arc-seconds for latitude, 22.5 arc-seconds for longitude) (11 digits)
# ABBBBBCCDDEF : 250m grid square code (7.5 arc-seconds for latitude, 11.25 arc-seconds for longitude) (12 digits)
# ABBBBBCCDDEFG : 125m grid square code (3.75 arc-seconds for latitude, 5.625 arc-seconds for longitude) (13 digits)
# ABBBBBCCDDHH : Extended 100m grid square code with 12 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (12 digits)
# ABBBBBCCDDEHH : Extended 100m grid square code with 13 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (13 digits)
# ABBBBBCCDDHHE : Extended 50m grid square code with 13 digits (1.5 arc-seconds for latitude, 2.25 arc-seconds for longitude) (13 digits)
# ABBBBBCCDDEHHI : Extended 50m grid square code with 14 digits (1.5 arc-seconds for latitude, 2.25 arc-seconds for longitude) (14 digits)
# ABBBBBCCDDHHGG : Extended 10m grid square code with 14 digits (0.3 arc-seconds for latitude, 0.45 arc-seconds for longitude) (14 digits)
# ABBBBBCCDDHHGGII : Extended 1m grid square code with 16 digits (0.03 arc-seconds for latitude, 0.045 arc-seconds for longitude) (16 digits)
#
# 3.
#
# Calculate geodesic distance and size of world grid square 
#
# T. Vincenty, ``Direct and Inverse Solutions of Geodesics
# on the Ellipsoid with application of nested equations'',
# Survey Review XXIII, Vol 176 (1975) Vol. 88-93.
#
# Vincenty(latitude1, longitude1, latitude2, longitude)
# : calculate geodesitc distance between two points (latitude1, longitude1) and (latitude2, longitude2) placed on the WGS84 Earth ellipsoid based on the Vincenty's formulae (1975)
# cal_area_from_meshcode(meshcode,extension=F,spec="")
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoide A) of world grid square indicated by meshcode
# cal_area_from_latlong(latlong)
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoid A) of a trapezoid on the WGS84 Earth ellipoid indicated by (latlong$lat0, latlong$long0, latlong$lat1, latlong$long1)
#


###############################################################################
# meshcode_to_latlong(meshcode; extension=false, spec="")
#
# Calculate northern western geographic position (latitude, longitude)
# of the grid from meshcode.Returns a named tuple with :lat and :long keys.
###############################################################################

function meshcode_to_latlong(meshcode; extension=false, spec="")
    res = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    return (lat=res.lat0, long=res.long0)
end

###############################################################################
#    meshcode_to_latlong_NW(meshcode; extension=false, spec="")
#
# Calculate northern western geographic position (latitude, longitude)
# of the grid from meshcode.
###############################################################################

function meshcode_to_latlong_NW(meshcode; extension=false, spec="")
    res = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    return (lat=res.lat0, long=res.long0)
end

###############################################################################
#    meshcode_to_latlong_SW(meshcode; extension=false, spec="")
#
# Calculate southern western geographic position (latitude, longitude)
# of the grid from meshcode.
###############################################################################

function meshcode_to_latlong_SW(meshcode; extension=false, spec="")
    res = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    return (lat=res.lat1, long=res.long0)
end

###############################################################################
#    meshcode_to_latlong_NE(meshcode; extension=false, spec="")
#
# Calculate northern eastern geographic position (latitude, longitude)
# of the grid from meshcode.
###############################################################################

function meshcode_to_latlong_NE(meshcode; extension=false, spec="")
    res = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    return (lat=res.lat0, long=res.long1)
end

###############################################################################
#    meshcode_to_latlong_SE(meshcode; extension=false, spec="")
#
# Calculate southern eastern geographic position (latitude, longitude)
# of the grid from meshcode.
###############################################################################

function meshcode_to_latlong_SE(meshcode; extension=false, spec="")
    res = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    return (lat=res.lat1, long=res.long1)
end

###############################################################################
#    meshcode_to_latlong_grid(meshcode; extension=false, spec="")
#
# Calculate northern western and southern eastern geographic positions of the grid.
# Returns a named tuple with :lat0, :long0, :lat1, :long1 keys.
###############################################################################

function meshcode_to_latlong_grid(meshcode; extension=false, spec="")
    code = string(meshcode)
    
    # Check if code has at least 6 digits
    if !occursin(r"^[0-9]{6}", code)
        return nothing
    end
    
    code0 = parse(Int, code[1:1]) - 1  # transforming code0 from 0 to 7
    code12_str = code[2:4]
    
    if startswith(code12_str, "00")
        code12 = parse(Int, code[4:4])
    elseif startswith(code12_str, "0")
        code12 = parse(Int, code[3:4])
    else
        code12 = parse(Int, code[2:4])
    end
    
    if code[5:5] == "0"
        code34 = parse(Int, code[6:6])
    else
        code34 = parse(Int, code[5:6])
    end
    
    # Extract codes for higher resolution grids
    if length(code) >= 8 && occursin(r"^[0-9]{8}", code)
        code5 = parse(Int, code[7:7])
        code6 = parse(Int, code[8:8])
    end
    
    if length(code) >= 10 && occursin(r"^[0-9]{10}", code)
        code7 = parse(Int, code[9:9])
        code8 = parse(Int, code[10:10])
    end
    
    # Handle extension and standard codes
    if !extension
        if length(code) >= 11 && occursin(r"^[0-9]{11}", code)
            code9 = parse(Int, code[11:11])
        end
        if length(code) >= 12 && occursin(r"^[0-9]{12}", code)
            code10 = parse(Int, code[12:12])
        end
        if length(code) >= 13 && occursin(r"^[0-9]{13}", code)
            code11 = parse(Int, code[13:13])
        end
    else
        if length(code) >= 12 && occursin(r"^[0-9]{12}", code)
            codeex9 = parse(Int, code[11:11])
            codeex10 = parse(Int, code[12:12])
        end
        if length(code) >= 13 && occursin(r"^[0-9]{13}", code)
            code9 = parse(Int, code[11:11])
            codeex10 = parse(Int, code[12:12])
            codeex11 = parse(Int, code[13:13])
        end
        if length(code) >= 14 && occursin(r"^[0-9]{14}", code)
            codeex11 = parse(Int, code[13:13])
            codeex12 = parse(Int, code[14:14])
        end
        if length(code) >= 16 && occursin(r"^[0-9]{16}", code)
            codeex13 = parse(Int, code[15:15])
            codeex14 = parse(Int, code[16:16])
        end
    end
    
    # Calculate x, y, z coordinates
    z = code0 % 2
    y = div(code0 - z, 2) % 2
    x = div(code0 - 2*y - z, 4)
    
    # Process based on code length
    if length(code) == 6
        # 1st grid (80km)
        lat0 = (code12 - x + 1) * 2 / 3
        long0 = code34 + y + 100*z
        lat0 = (1 - 2*x) * lat0
        long0 = (1 - 2*y) * long0
        dlat = 2/3
        dlong = 1
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 8
        # 2nd grid (10km)
        lat0 = code12 * 2 / 3
        long0 = code34 + 100*z
        lat0 = lat0 + ((code5 - x + 1) * 2 / 3) / 8
        long0 = long0 + (code6 + y) / 8
        lat0 = (1 - 2*x) * lat0
        long0 = (1 - 2*y) * long0
        dlat = 2/3/8
        dlong = 1/8
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 10
        # 3rd grid (1km)
        lat0 = code12 * 2 / 3
        long0 = code34 + 100*z
        lat0 = lat0 + (code5 * 2 / 3) / 8
        long0 = long0 + code6 / 8
        lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
        long0 = long0 + (code8 + y) / 8 / 10
        lat0 = (1 - 2*x) * lat0
        long0 = (1 - 2*y) * long0
        dlat = 2/3/8/10
        dlong = 1/8/10
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 11
        # 4th grid (500m)
        lat0 = code12 * 2 / 3
        long0 = code34 + 100*z
        lat0 = lat0 + (code5 * 2 / 3) / 8
        long0 = long0 + code6 / 8
        lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
        long0 = long0 + (code8 + y) / 8 / 10
        lat0 = lat0 + (floor(Int, (code9 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2
        long0 = long0 + ((code9 - 1) % 2 - y) / 8 / 10 / 2
        lat0 = (1 - 2*x) * lat0
        long0 = (1 - 2*y) * long0
        dlat = 2/3/8/10/2
        dlong = 1/8/10/2
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 12
        # 5th grid (250m) or Extended 100m
        if !extension
            lat0 = code12 * 2 / 3
            long0 = code34 + 100*z
            lat0 = lat0 + (code5 * 2 / 3) / 8
            long0 = long0 + code6 / 8
            lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
            long0 = long0 + (code8 + y) / 8 / 10
            lat0 = lat0 + (floor(Int, (code9 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2
            long0 = long0 + ((code9 - 1) % 2 - y) / 8 / 10 / 2
            lat0 = lat0 + (floor(Int, (code10 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2 / 2
            long0 = long0 + ((code10 - 1) % 2 - y) / 8 / 10 / 2 / 2
            lat0 = (1 - 2*x) * lat0
            long0 = (1 - 2*y) * long0
            dlat = 2/3/8/10/2/2
            dlong = 1/8/10/2/2
        else
            # Extended 100m grid square code (12 digits)
            lat0 = code12 * 2 / 3
            long0 = code34 + 100*z
            lat0 = lat0 + (code5 * 2 / 3) / 8
            long0 = long0 + code6 / 8
            lat0 = lat0 + (code7 * 2 / 3) / 8 / 10
            long0 = long0 + code8 / 8 / 10
            lat0 = lat0 + ((codeex9 - x + 1) * 2 / 3) / 8 / 10 / 10
            long0 = long0 + (codeex10 + y) / 8 / 10 / 10
            lat0 = (1 - 2*x) * lat0
            long0 = (1 - 2*y) * long0
            dlat = 2/3/8/10/10
            dlong = 1/8/10/10
        end
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 13
        # 6th grid (125m) or Extended variants
        if !extension
            lat0 = code12 * 2 / 3
            long0 = code34 + 100*z
            lat0 = lat0 + (code5 * 2 / 3) / 8
            long0 = long0 + code6 / 8
            lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
            long0 = long0 + (code8 + y) / 8 / 10
            lat0 = lat0 + (floor(Int, (code9 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2
            long0 = long0 + ((code9 - 1) % 2 - y) / 8 / 10 / 2
            lat0 = lat0 + (floor(Int, (code10 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2 / 2
            long0 = long0 + ((code10 - 1) % 2 - y) / 8 / 10 / 2 / 2
            lat0 = lat0 + (floor(Int, (code11 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2 / 2 / 2
            long0 = long0 + ((code11 - 1) % 2 - y) / 8 / 10 / 2 / 2 / 2
            lat0 = (1 - 2*x) * lat0
            long0 = (1 - 2*y) * long0
            dlat = 2/3/8/10/2/2/2
            dlong = 1/8/10/2/2/2
        else
            if spec == "ex100m_13"
                # Extended 100m grid square code with 13 digits
                lat0 = code12 * 2 / 3
                long0 = code34 + 100*z
                lat0 = lat0 + (code5 * 2 / 3) / 8
                long0 = long0 + code6 / 8
                lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
                long0 = long0 + (code8 + y) / 8 / 10
                lat0 = lat0 + (floor(Int, (code9 - 1)/2) + 2*x - 2) * 2 / 3 / 8 / 10 / 2
                long0 = long0 + ((code9 - 1) % 2 - 2*y) / 8 / 10 / 2
                lat0 = lat0 + ((codeex10 - x + 1) * 2 / 3) / 8 / 10 / 2 / 5
                long0 = long0 + (codeex11 + y) / 8 / 10 / 2 / 5
                lat0 = (1 - 2*x) * lat0
                long0 = (1 - 2*y) * long0
                dlat = 2/3/8/10/2/5
                dlong = 1/8/10/2/5
            else
                # Extended 50m grid square code with 13 digits
                lat0 = code12 * 2 / 3
                long0 = code34 + 100*z
                lat0 = lat0 + (code5 * 2 / 3) / 8
                long0 = long0 + code6 / 8
                lat0 = lat0 + (code7 * 2 / 3) / 8 / 10
                long0 = long0 + code8 / 8 / 10
                lat0 = lat0 + ((codeex9 - x + 1) * 2 / 3) / 8 / 10 / 10
                long0 = long0 + (codeex10 + y) / 8 / 10 / 10
                lat0 = lat0 + (floor(Int, (codeex11 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 10 / 2
                long0 = long0 + ((codeex11 - 1) % 2 - y) / 8 / 10 / 10 / 2
                lat0 = (1 - 2*x) * lat0
                long0 = (1 - 2*y) * long0
                dlat = 2/3/8/10/10/2
                dlong = 1/8/10/10/2
            end
        end
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 14
        # Extended 50m or 10m grid square code (14 digits)
        if !extension
            return nothing
        else
            if spec == "ex50m_14"
                # Extended 50m grid square code with 14 digits
                lat0 = code12 * 2 / 3
                long0 = code34 + 100*z
                lat0 = lat0 + (code5 * 2 / 3) / 8
                long0 = long0 + code6 / 8
                lat0 = lat0 + ((code7 - x + 1) * 2 / 3) / 8 / 10
                long0 = long0 + (code8 + y) / 8 / 10
                lat0 = lat0 + (floor(Int, (code9 - 1)/2) + 2*x - 2) * 2 / 3 / 8 / 10 / 2
                long0 = long0 + ((code9 - 1) % 2 - 2*y) / 8 / 10 / 2
                lat0 = lat0 + ((codeex10 - x + 1) * 2 / 3) / 8 / 10 / 2 / 5
                long0 = long0 + (codeex11 + y) / 8 / 10 / 2 / 5
                lat0 = lat0 + (floor(Int, (codeex12 - 1)/2) + x - 1) * 2 / 3 / 8 / 10 / 2 / 5 / 2
                long0 = long0 + ((codeex12 - 1) % 2 - y) / 8 / 10 / 2 / 5 / 2
                lat0 = (1 - 2*x) * lat0
                long0 = (1 - 2*y) * long0
                dlat = 2/3/8/10/2/5/2
                dlong = 1/8/10/2/5/2
            else
                # Extended 10m grid square code with 14 digits
                lat0 = code12 * 2 / 3
                long0 = code34 + 100*z
                lat0 = lat0 + (code5 * 2 / 3) / 8
                long0 = long0 + code6 / 8
                lat0 = lat0 + code7 * 2 / 3 / 8 / 10
                long0 = long0 + code8 / 8 / 10
                lat0 = lat0 + (codeex9 * 2 / 3) / 8 / 10 / 10
                long0 = long0 + codeex10 / 8 / 10 / 10
                lat0 = lat0 + ((codeex11 - x + 1) * 2 / 3) / 8 / 10 / 10 / 10
                long0 = long0 + (codeex12 + y) / 8 / 10 / 10 / 10
                lat0 = (1 - 2*x) * lat0
                long0 = (1 - 2*y) * long0
                dlat = 2/3/8/10/10/10
                dlong = 1/8/10/10/10
            end
        end
        lat1 = lat0 - dlat
        long1 = long0 + dlong
        return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        
    elseif length(code) == 16
        # Extended 1m grid square code (16 digits)
        if !extension
            return nothing
        else
            lat0 = code12 * 2 / 3
            long0 = code34 + 100*z
            lat0 = lat0 + (code5 * 2 / 3) / 8
            long0 = long0 + code6 / 8
            lat0 = lat0 + (code7 * 2 / 3) / 8 / 10
            long0 = long0 + code8 / 8 / 10
            lat0 = lat0 + (codeex9 * 2 / 3) / 8 / 10 / 10
            long0 = long0 + codeex10 / 8 / 10 / 10
            lat0 = lat0 + (codeex11 * 2 / 3) / 8 / 10 / 10 / 10
            long0 = long0 + codeex12 / 8 / 10 / 10 / 10
            lat0 = lat0 + ((codeex13 - x + 1) * 2 / 3) / 8 / 10 / 10 / 10 / 10
            long0 = long0 + (codeex14 + y) / 8 / 10 / 10 / 10 / 10
            lat0 = (1 - 2*x) * lat0
            long0 = (1 - 2*y) * long0
            dlat = 2/3/8/10/10/10/10
            dlong = 1/8/10/10/10/10
            lat1 = lat0 - dlat
            long1 = long0 + dlong
            return (lat0=lat0, long0=long0, lat1=lat1, long1=long1)
        end
    end
    
    return nothing
end

# ============================================================================
# Calculate meshcode from latitude and longitude
# ============================================================================

###############################################################################
#    cal_meshcode(latitude, longitude)
#
# Calculate a basic (1km) grid square code (10 digits) from a geographical position.
# Default is cal_meshcode3.
###############################################################################

function cal_meshcode(latitude, longitude)
    return cal_meshcode3(latitude, longitude)
end

###############################################################################
#    cal_meshcode1(latitude, longitude)
#
# Calculate an 80km grid square code (6 digits) from a geographical position.
###############################################################################

function cal_meshcode1(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    u = floor(Int, longitude - 100*z)
    
    return format_meshcode(o, p, u, "")
end

###############################################################################
#    cal_meshcode2(latitude, longitude)
#
# Calculate a 10km grid square code (8 digits) from a geographical position.
###############################################################################

function cal_meshcode2(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "99999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = (latitude * 60 / 40 - p) * 40
    q = floor(Int, a / 5)
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    
    return format_meshcode(o, p, u, "$(q)$(v)")
end

###############################################################################
#    cal_meshcode3(latitude, longitude)
#
# Calculate a 1km grid square code (10 digits) from a geographical position.
###############################################################################

function cal_meshcode3(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "9999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = (latitude * 60 / 40 - p) * 40
    q = floor(Int, a / 5)
    b = (a / 5 - q) * 5
    r = floor(Int, b * 60 / 30)
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = (f * 60 / 7.5 - v) * 7.5
    w = floor(Int, g * 60 / 45)
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)")
end

###############################################################################
#    cal_meshcode4(latitude, longitude)
#
# Calculate a 1km grid square code (11 digits) from a geographical position
# (500m resolution).
###############################################################################

function cal_meshcode4(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "99999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = (latitude * 60 / 40 - p) * 40
    q = floor(Int, a / 5)
    b = (a / 5 - q) * 5
    r = floor(Int, b * 60 / 30)
    c = (b * 60 / 30 - r) * 30
    s2u = floor(Int, c / 15)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = (f * 60 / 7.5 - v) * 7.5
    w = floor(Int, g * 60 / 45)
    h = (g * 60 / 45 - w) * 45
    s2l = floor(Int, h / 22.5)
    s2 = s2u * 2 + s2l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s2)")
end

###############################################################################
#    cal_meshcode5(latitude, longitude)
#
# Calculate a 500m grid square code (12 digits) from a geographical position.
###############################################################################

function cal_meshcode5(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = (latitude * 60 / 40 - p) * 40
    q = floor(Int, a / 5)
    b = (a / 5 - q) * 5
    r = floor(Int, b * 60 / 30)
    c = (b * 60 / 30 - r) * 30
    s2u = floor(Int, c / 15)
    d = (c / 15 - s2u) * 15
    s4u = floor(Int, d / 7.5)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = (f * 60 / 7.5 - v) * 7.5
    w = floor(Int, g * 60 / 45)
    h = (g * 60 / 45 - w) * 45
    s2l = floor(Int, h / 22.5)
    i = (h / 22.5 - s2l) * 22.5
    s4l = floor(Int, i / 11.25)
    s2 = s2u * 2 + s2l + 1
    s4 = s4u * 2 + s4l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s2)$(s4)")
end

###############################################################################
#    cal_meshcode6(latitude, longitude)
#
# Calculate a 250m grid square code (13 digits) from a geographical position.
###############################################################################

function cal_meshcode6(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "9999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s2u = floor(Int, c / 15)
    d = c - s2u * 15
    s4u = floor(Int, d / 7.5)
    e = d - s4u * 7.5
    s8u = floor(Int, e / 3.75)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    s2l = floor(Int, h / 22.5)
    i = h - s2l * 22.5
    s4l = floor(Int, i / 11.25)
    j = i - s4l * 11.25
    s8l = floor(Int, j / 5.625)
    s2 = s2u * 2 + s2l + 1
    s4 = s4u * 2 + s4l + 1
    s8 = s8u * 2 + s8l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s2)$(s4)$(s8)")
end

###############################################################################
#    cal_meshcode_ex100m_12(latitude, longitude)
#
# Calculate an extended 100m grid square code (12 digits).
# 3 arc-second for latitude, 4.5 arc-second for longitude.
###############################################################################

function cal_meshcode_ex100m_12(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s = floor(Int, c / 3)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    xx = floor(Int, h / 4.5)
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s)$(xx)")
end

###############################################################################
#    cal_meshcode_ex100m_13(latitude, longitude)
#
# Calculate an extended 100m grid square code (13 digits).
# 3 arc-second for latitude, 4.5 arc-second for longitude.
###############################################################################

function cal_meshcode_ex100m_13(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "9999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s2u = floor(Int, c / 15)
    d = c - s2u * 15
    et = floor(Int, d / 3)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    s2l = floor(Int, h / 22.5)
    i = h - s2l * 22.5
    jt = floor(Int, i / 4.5)
    s2 = s2u * 2 + s2l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s2)$(et)$(jt)")
end

###############################################################################
#    cal_meshcode_ex50m_13(latitude, longitude)
#
# Calculate an extended 50m grid square code (13 digits).
# 1.5 arc-second for latitude, 2.25 arc-second for longitude.
###############################################################################

function cal_meshcode_ex50m_13(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "9999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s = floor(Int, c / 3)
    d = c - s * 3
    t2u = floor(Int, d / 1.5)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    xx = floor(Int, h / 4.5)
    i = h - xx * 4.5
    t2l = floor(Int, i / 2.25)
    tt = t2u * 2 + t2l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s)$(xx)$(tt)")
end

###############################################################################
#    cal_meshcode_ex50m_14(latitude, longitude)
#
# Calculate an extended 50m grid square code (14 digits).
# 1.5 arc-second for latitude, 2.25 arc-second for longitude.
###############################################################################

function cal_meshcode_ex50m_14(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "99999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s2u = floor(Int, c / 15)
    d = c - s2u * 15
    et = floor(Int, d / 3)
    e = d - et * 3
    t2u = floor(Int, e / 1.5)
    
    u = floor(Int, longitude - 100*z)
    f_temp = longitude - 100*z - u
    v = floor(Int, f_temp * 60 / 7.5)
    g = f_temp * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    s2l = floor(Int, h / 22.5)
    i = h - s2l * 22.5
    jt = floor(Int, i / 4.5)
    j = i - jt * 4.5
    t2l = floor(Int, j / 2.25)
    s2 = s2u * 2 + s2l + 1
    t2 = t2u * 2 + t2l + 1
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s2)$(et)$(jt)$(t2)")
end

###############################################################################
#    cal_meshcode_ex10m_14(latitude, longitude)
#
# Calculate an extended 10m grid square code (14 digits).
# 0.3 arc-second for latitude, 0.45 arc-second for longitude.
###############################################################################

function cal_meshcode_ex10m_14(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "99999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s = floor(Int, c / 3)
    d = c - s * 3
    t = floor(Int, d / 0.3)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    xx = floor(Int, h / 4.5)
    i = h - xx * 4.5
    yy = floor(Int, i / 0.45)
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s)$(xx)$(t)$(yy)")
end

###############################################################################
#    cal_meshcode_ex1m_16(latitude, longitude)
#
# Calculate an extended 1m grid square code (16 digits).
# 0.03 arc-second for latitude, 0.045 arc-second for longitude.
###############################################################################

function cal_meshcode_ex1m_16(latitude, longitude)
    if latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180
        return "9999999999999999"
    end
    
    if latitude < 0
        o = 4
    else
        o = 0
    end
    
    if longitude < 0
        o = o + 2
    end
    
    if abs(longitude) >= 100
        o = o + 1
    end
    
    z = o % 2
    y = div(o - z, 2) % 2
    x = div(o - 2*y - z, 4)
    o = o + 1
    
    latitude = (1 - 2*x) * latitude
    longitude = (1 - 2*y) * longitude
    
    p = floor(Int, latitude * 60 / 40)
    a = latitude * 60 - p * 40
    q = floor(Int, a / 5)
    b = a - q * 5
    r = floor(Int, b * 60 / 30)
    c = b * 60 - r * 30
    s = floor(Int, c / 3)
    d = c - s * 3
    t = floor(Int, d / 0.3)
    e = d - t * 0.3
    tt = floor(Int, e / 0.03)
    
    u = floor(Int, longitude - 100*z)
    f = longitude - 100*z - u
    v = floor(Int, f * 60 / 7.5)
    g = f * 60 - v * 7.5
    w = floor(Int, g * 60 / 45)
    h = g * 60 - w * 45
    xx = floor(Int, h / 4.5)
    i = h - xx * 4.5
    yy = floor(Int, i / 0.45)
    j = i - yy * 0.45
    zz = floor(Int, j / 0.045)
    
    return format_meshcode(o, p, u, "$(q)$(v)$(r)$(w)$(s)$(xx)$(t)$(yy)$(tt)$(zz)")
end

# ============================================================================
# Helper function for formatting meshcode
# ============================================================================

function format_meshcode(o::Int, p::Int, u::Int, extra::String)
    """Helper function to format meshcode with proper zero padding."""
    if u < 10
        if p < 10
            return "$(o)00$(p)0$(u)$(extra)"
        elseif p < 100
            return "$(o)0$(p)0$(u)$(extra)"
        else
            return "$(o)$(p)0$(u)$(extra)"
        end
    else
        if p < 10
            return "$(o)00$(p)$(u)$(extra)"
        elseif p < 100
            return "$(o)0$(p)$(u)$(extra)"
        else
            return "$(o)$(p)$(u)$(extra)"
        end
    end
end

# ============================================================================
# Calculate geodesic distance and area
# ============================================================================

###############################################################################
#    Vincenty(latitude1, longitude1, latitude2, longitude2)
# 
# Calculate geodesic distance between two points on the WGS84 Earth ellipsoid
# using Vincenty's formulae (1975).
#
# Returns distance in meters.
###############################################################################

function Vincenty(latitude1::Float64, longitude1::Float64, latitude2::Float64, longitude2::Float64)
    # WGS84 parameters
    f = 1.0 / 298.257223563
    a = 6378137.0
    b = 6356752.314245
    
    # Check if two places are equivalent
    if abs(latitude1 - latitude2) < 1e-18 && abs(longitude1 - longitude2) < 1e-18
        return 0.0
    end
    
    # Convert degrees to radians
    lat1_rad = latitude1 * π / 180.0
    lat2_rad = latitude2 * π / 180.0
    lon1_rad = longitude1 * π / 180.0
    lon2_rad = longitude2 * π / 180.0
    
    # Initial values
    L = lon1_rad - lon2_rad
    U1 = atan((1.0 - f) * tan(lat1_rad))
    U2 = atan((1.0 - f) * tan(lat2_rad))
    
    lambda = L
    dlambda = 10.0
    lambda_prev = 0.0
    
    # Initialize variables for use after loop
    sigma = 0.0
    sinsigma = 0.0
    cossigma = 0.0
    sinalpha = 0.0
    cos2alpha = 0.0
    cos2sigmam = 0.0
    C = 0.0
    
    # Iterative calculation (Vincenty's forward calculation)
    max_iterations = 1000
    iteration = 0
    while abs(dlambda) > 1e-12 && iteration < max_iterations
        iteration += 1
        
        cs = cos(U2) * sin(lambda)
        cscc = cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lambda)
        sinsigma = sqrt(cs^2 + cscc^2)
        
        # Handle the case where sinsigma is zero (antipodal points)
        if sinsigma < 1e-18
            return 0.0
        end
        
        cossigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lambda)
        sigma = atan(sinsigma / cossigma)
        sinalpha = cos(U1) * cos(U2) * sin(lambda) / sinsigma
        cos2alpha = 1.0 - sinalpha * sinalpha
        
        if abs(cos2alpha) < 1e-18
            # Equatorial line
            C = 0.0
            cos2sigmam = 0.0
            lambda_prev = lambda
            lambda = L + f * sinalpha * sigma
        else
            cos2sigmam = cossigma - 2.0 * sin(U1) * sin(U2) / cos2alpha
            C = f / 16.0 * cos2alpha * (4.0 + f * (4.0 - 3.0 * cos2alpha))
            lambda_prev = lambda
            lambda = L + (1.0 - C) * f * sinalpha * (sigma + C * sinsigma * (cos2sigmam + C * cossigma * (-1.0 + 2.0 * cos2sigmam * cos2sigmam)))
        end
        
        dlambda = lambda - lambda_prev
    end
    
    # Calculate distance
    if abs(cos2alpha) < 1e-18
        A = 1.0
        dsigma = 0.0
    else
        u2 = cos2alpha * (a * a - b * b) / (b * b)
        A = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)))
        B = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)))
        dsigma = B * sinsigma * (cos2sigmam + 1.0/4.0 * B * (cossigma * (-1.0 + 2.0 * cos2sigmam * cos2sigmam) - 1.0/6.0 * B * cos2sigmam * (-3.0 + 4.0 * sinsigma * sinsigma) * (-3.0 + 4.0 * cos2sigmam * cos2sigmam)))
    end
    
    s = b * A * (sigma - dsigma)
    return s
end

###############################################################################
#    cal_area_from_meshcode(meshcode; extension=false, spec="")
#
# Calculate size of world grid square indicated by meshcode.
# Returns a named tuple with :W1, :W2, :H, :A (area) keys.
###############################################################################

function cal_area_from_meshcode(meshcode; extension=false, spec="")
    latlong = meshcode_to_latlong_grid(meshcode, extension=extension, spec=spec)
    if latlong === nothing
        return nothing
    end
    return cal_area_from_latlong(latlong)
end

###############################################################################
#    cal_area_from_latlong(latlong)
#
# Calculate size (northern and southern west-to-east span, north-to-south span, and area)
# of a trapezoid on the WGS84 Earth ellipsoid.
#
# Input: latlong named tuple with :lat0, :long0, :lat1, :long1 keys
# Returns: named tuple with :W1, :W2, :H, :A keys
###############################################################################

function cal_area_from_latlong(latlong)
    W1 = Vincenty(latlong.lat0, latlong.long0, latlong.lat0, latlong.long1)
    W2 = Vincenty(latlong.lat1, latlong.long0, latlong.lat1, latlong.long1)
    H = Vincenty(latlong.lat0, latlong.long0, latlong.lat1, latlong.long0)
    A = (W1 + W2) * H * 0.5
    
    return (W1=W1, W2=W2, H=H, A=A)
end
nothing

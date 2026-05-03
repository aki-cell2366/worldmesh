#
# Copyright (c) 2015-2022 Research Institute for World Grid Squares 
# Prof. Dr. Aki-Hiro Sato
# All rights reserved. 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Python functions to calculate the world grid square code.
# The world grid square code computed by this library
# is compatible to JISX0410.
#
# Release notes
# Version 1.0 : Released on 4 April 2017
# Version 1.1 : Released on 7 September 2017
# Version 1.2 : Released on 10 December 2018
# Version 1.4 : Released on 14 January 2020
# Version 1.6 : Released on 09 October 2020
# Version 1.7 : Released on 20 June 2021
# Version 1.75 : Released on 13 April 2022
# Version 1.76 : Released on 12 December 2022
# Version 1.8 : Released on 12 December 2022
#
# Written by Prof. Dr. Aki-Hiro Sato
# Department of Data Science, Graduate School of Data Science
# Yokohama City University
#
# Contact:
# Address: 22-2, Seto, Kanazawa-ku, Yokohama, Kanagawa 236-0027 Japan
# E-mail: ahsato@yokohama-cu.ac.jp
# TEL: +81-45-787-2208
#
# Difference from Version 1.76
# Adding two functions for exteded grid square codes with high spatial resolution such as
# cal_meshcode_ex50m_13() and cal_meshcode_ex50m_14().
# Difference from Version 1.75
# Debugging a code for cal_meshcode_ex1m16().
# Difference from Version 1.7
# Refinement of floating point calculation for cal_meshcode6(), cal_meshcode_ex100m_12(), cal_meshcode_ex100m_13(), cal_meshcode_ex10m_14(), and cal_meshcode_ex1m_16().
# Difference from Version 1.6
# Adding four functions for extended grid square codes with high spatial resolutions such as
# cal_meshcode_ex100m_12(), cal_meshcode_ex100m_13(), cal_meshcode_ex10m_14(), and cal_meshcode_ex1m_16(). 
# Difference from Version 1.2
# Three functions related to geodesic distance and representative size of world grid squares Vincenty(),
# cal_area_from_meshcode(), and cal_area_from_latlong() are added. 
# Difference from Version 1.1
# Compatibility of numeric calcuation for both Python2 and Python3
# Difference from Version 1.0
# Debugging for cal_meshcode5() and cal_meshcode6()
#
# Three types of functions are defined in this library.
# 1. calculate representative geographical position(s) (latitude, longitude) of a grid square from a grid square code
# 2. calculate a grid square code from a geographical position (latitude, longitude)
#
# 1.
#
# meshcode_to_latlong(meshcode, extension=False, spec="")
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NW(meshcode, extension=False, spec="")
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SW(meshcode, extension=False, spec="")
# : calculate sourthern western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NE(meshcode, extension=False, spec="")
# : calculate northern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SE(meshcode, extension=False, spec="")
# : calculate sourthern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_grid(meshcode, extension=False, spec="")
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
# : calculate an extended 100m (1km / 10) grid square code (12 digits) from a geographical position (latitude, longitude)
# - 3 arc-second for latitude and 4.5 arc-second for longitude
#
# cal_meshcode_ex100m_13(latitude,longitude)
# : calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude)
# - 3 arc-second for latitude and 4.5 arc-second for longitude 
#
# cal_meshcode_ex50m_13(latitude,longitude)
# : calculate an extended 50m (100m / 2) grid square code (13 digits) from a geographical position (latitude, longitude)
# - 1.5 arc-second for latitude and 2.25 arc-second for longitude 
#
# cal_meshcode_ex50m_14(latitude,longitude)
# : calculate an extended 50m (100m / 2) grid square code (14 digits) from a geographical position (latitude, longitude)
# - 1.5 arc-second for latitude and 2.25 arc-second for longitude 
#
# cal_meshcode_ex10m_14(latitude,longitude)
# : calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) from a geographical position (latitude, longitude)
# - 0.3 arc-second for latitude and 0.45 arc-second for longitude
#
# cal_meshcode_ex1m_16(latitude,longitude)
# : calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) from a geographical position (latitude, longitude)
# - 0.03 arc-second for latitude and 0.045 arc-second for longitude
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
# cal_area_from_meshcode(meshcode, extension=False, spec="")
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoide A) of world grid square indicated by meshcode
# cal_area_from_latlong(latlong)
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoid A) of a trapezoid on the WGS84 Earth ellipoid indicated by (latlong["lat0"], latlong["long0"], latlong["lat1"], latlong["long1"])
#

import math

def meshcode_to_latlong(meshcode, extension=False, spec=""):
    res=meshcode_to_latlong_grid(meshcode, extension, spec)
    xx={"lat":res["lat0"],"long":res["long0"]}
    return xx

def meshcode_to_latlong_NW(meshcode, extension=False, spec=""):
    res=meshcode_to_latlong_grid(meshcode, extension, spec)
    xx={"lat":res["lat0"],"long":res["long0"]}
    return xx

def meshcode_to_latlong_SW(meshcode, extension=False, spec=""):
  res=meshcode_to_latlong_grid(meshcode, extension, spec)
  xx={"lat":res["lat1"],"long":res["long0"]}
  return(xx)

def meshcode_to_latlong_NE(meshcode, extension=False, spec=""):
  res=meshcode_to_latlong_grid(meshcode, extension, spec)
  xx={"lat":res["lat0"],"long":res["long1"]}
  return xx

def meshcode_to_latlong_SE(meshcode, extension=False, spec=""):
  res=meshcode_to_latlong_grid(meshcode, extension, spec)
  xx={"lat":res["lat1"],"long":res["long1"]}
  return xx

def meshcode_to_latlong_grid(meshcode, extension=False, spec=""):
    # The extension option determines high resolution gird square codes if extension = T
    # 100m x 100m 3 arc-second x 4.5 arc second (12 digits) : default
    # 100m x 100m 3 arc-second x 4.5 arc second (13 digits) : spec = "ex100m_13"
    # 50m x 50m 1.5 arc-second x 2.25 arc second (13 digits) : default
    # 50m x 50m 1.5 arc-second x 2.25 arc second (14 digits) : sepc = "ex50m_14"
    # 10m x 10m 0.3 arc-second x 0.45 arc second (14 digits) : default
    # 1m x 1m 0.03 arc-second x 0.045 arc second (16 digits) : default
    code=str(meshcode)
    # more than 1st grid
    if len(code)>=6 :
        code0=int(code[0:1])-1
        code12=int(code[1:4])
        code34=int(code[4:6])
    else: return None

    # more than 2nd grid
    if len(code)>=8:
        code5=int(code[6:7])
        code6=int(code[7:8])
    # more than 3rd grid
    if len(code)>=10:
        code7=int(code[8:9])
        code8=int(code[9:10])
    if(not extension):
        # more than 4th grid
        if len(code)>=11:
            code9=int(code[10:11])
        # morethan 5th grid
        if len(code)>=12:
            code10=int(code[11:12])
        # more than 6th grid        
        if len(code)>=13:
            code11=int(code[12:13])
    else:
        # Extended 100m grid square code (12 digits)
        if len(code)>=12:
            codeex9=int(code[10:11])
            codeex10=int(code[11:12])
        # Extended 100m or 50m grid square code (13 digits)
        if len(code)>=13:
            code9=int(code[10:11])            
            codeex10=int(code[11:12])
            codeex11=int(code[12:13])
        # Extended 50m or 10m grid square code (14 digits)
        if len(code)>=14:
            codeex11=int(code[12:13])
            codeex12=int(code[13:14])
        # Extended 1m grid square code (16 digits)
        if len(code)>=16:
            codeex13=int(code[14:15])
            codeex14=int(code[15:16])
            
    # 0'th grid
    z = code0 % 2
    y = ((code0-z)/2) % 2
    x = (code0-2*y-z)/4

    # 1st grid
    if len(code)==6:
        lat0  = (code12-x+1) * 2.0 / 3.0
        long0 = (code34+y) + 100*z
        lat0  = (1-2*x)*lat0
        long0 = (1-2*y)*long0
        dlat = 2.0/3.0
        dlong = 1.0
        lat1  = "%.8f" % (lat0-dlat)
        long1 = "%.8f" % (long0+dlong)
        lat0  = "%.8f" % lat0
        long0 = "%.8f" % long0
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
    # 2nd grid
    elif len(code)==8:
        #    Code5
        #      7********
        #      6********    
        #      5******** 
        #      4********
        #      3********
        #      2********
        #      1********
        #      0********
        #       01234567 Code6
        lat0  = code12 * 2.0 / 3.0
        long0 = code34 + 100*z
        lat0  = lat0  + ((code5-x+1) * 2.0 / 3.0) / 8.0
        long0 = long0 +  (code6+y) / 8.0
        lat0 = (1-2*x) * lat0
        long0 = (1-2*y) * long0
        dlat = 2.0/3.0/8.0
        dlong = 1.0/8.0
        lat1  = "%.8f" % (lat0-dlat)
        long1 = "%.8f" % (long0+dlong)
        lat0  = "%.8f" % lat0
        long0 = "%.8f" % long0
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
     # 3rd grid
    elif len(code)==10:
        #    Code7
        #      9**********    
        #      8**********    
        #      7**********
        #      6**********    
        #      5********** 
        #      4**********
        #      3**********
        #      2**********
        #      1**********
        #      0**********
        #       0123456789  Code8
        lat0  = code12 * 2.0 / 3.0
        long0 = code34 + 100*z
        lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
        long0 = long0 +  code6 / 8.0
        lat0  = lat0  + ((code7-x+1) * 2.0 / 3.0) / 8.0 / 10.0
        long0 = long0 +  (code8+y) / 8.0 / 10.0
        lat0 = (1-2*x)*lat0
        long0 = (1-2*y)*long0
        dlat = 2.0/3.0/8.0/10.0
        dlong = 1/8.0/10.0
        lat1  = "%.8f" % (lat0-dlat)
        long1 = "%.8f" % (long0+dlong)
        lat0  = "%.8f" % (lat0)
        long0 = "%.8f" % (long0)
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
    # code 9
    #     N
    #   3 | 4
    # W - + - E
    #   1 | 2
    #     S
    # 4th frid (11 digits)
    elif len(code)==11:
         lat0  = code12 * 2.0 / 3.0
         long0 = code34 + 100*z
         lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
         long0 = long0 +  code6 / 8.0
         lat0  = lat0  + ((code7-x+1) * 2.0 / 3.0) / 8.0 / 10.0
         long0 = long0 +  (code8+y) / 8.0 / 10.0
         lat0  = lat0  + (math.floor((code9-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0
         long0 = long0 + ((code9-1)%2-y) / 8.0 / 10.0 / 2.0
         lat0 = (1-2*x)*lat0
         long0 = (1-2*y)*long0
         dlat = 2.0/3.0/8.0/10.0/2.0
         dlong = 1.0/8.0/10.0/2.0
         lat1  = "%.8f" % (lat0-dlat)
         long1 = "%.8f" % (long0+dlong)
         lat0  = "%.8f" % lat0
         long0 = "%.8f" % long0
         xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
         return xx
    elif len(code)==12:
        if(not extension):
            # code 10
            #     N
            #   3 | 4
            # W - + - E
            #   1 | 2
            #     S
            # 5th grid (12 digits)
            lat0  = code12 * 2.0 / 3.0
            long0 = code34 + 100*z
            lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
            long0 = long0 +  code6 / 8.0
            lat0  = lat0  + ((code7-x+1) * 2.0 / 3.0) / 8.0 / 10.0
            long0 = long0 +  (code8+y) / 8.0 / 10.0
            lat0  = lat0  + (math.floor((code9-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0
            long0 = long0 + ((code9-1)%2-y) / 8.0 / 10.0 / 2.0
            lat0  = lat0  + (math.floor((code10-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0
            long0 = long0 + ((code10-1)%2-y) / 8.0 / 10.0 / 2.0 / 2.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/2.0/2.0
            dlong = 1.0/8.0/10.0/2.0/2.0
        else:
            # Extended 100m grid square code (12 digits)
            #   Code9
            #      9**********
            #      8**********
            #      7**********
            #      6**********
            #      5**********
            #      4**********
            #      3**********
            #      2**********
            #      1**********
            #      0**********
            #      0123456789 Code10
            lat0  = code12 * 2.0 / 3.0
            long0 = code34 + 100*z
            lat0  = lat0 + code5 * 2.0 / 3.0 / 8.0
            long0 = long0 + code6 / 8.0
            lat0  = lat0 + code7 * 2.0 / 3.0 / 8.0 / 10.0
            long0 = long0 + code8 / 8.0 / 10.0
            lat0  = lat0 + ((codeex9-x+1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0
            long0 = long0 + (codeex10+y) / 8.0 / 10.0 / 10.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/10.0
            dlong = 1.0/8.0/10.0/10.0
            
        lat1  = "%.10f" % (lat0-dlat)
        long1 = "%.10f" % (long0+dlong)
        lat0  = "%.10f" % (lat0)
        long0 = "%.10f" % (long0)
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
    
    elif len(code)==13:
        if(not extension):
            # 6rd grid (13 digits)
            # code 11
            #     N
            #   3 | 4
            # W - + - E
            #   1 | 2
            #     S
            # 6rd grid (13 digits)        
            lat0  = code12 * 2.0 / 3.0
            long0 = code34 + 100*z
            lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
            long0 = long0 + code6 / 8.0
            lat0  = lat0  + ((code7-x+1) * 2.0 / 3.0) / 8.0 / 10.0
            long0 = long0 + (code8+y) / 8.0 / 10.0
            lat0  = lat0  + (math.floor((code9-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0
            long0 = long0 + ((code9-1)%2-y) / 8.0 / 10.0 / 2.0
            lat0  = lat0  + (math.floor((code10-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0
            long0 = long0 + ((code10-1)%2-y) / 8.0 / 10.0 / 2.0 / 2.0
            lat0  = lat0  + (math.floor((code11-1)/2)+x-1) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 2.0 / 2.0
            long0 = long0 + ((code11-1)%2-y) / 8.0 / 10.0 / 2.0 / 2.0 / 2.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/2.0/2.0/2.0
            dlong = 1.0/8.0/10.0/2.0/2.0/2.0
        elif spec == "ex100m_13":
            # Extended 100m grid square code (13 digits)
            #    Code10
            #      4*****
            #      3*****
            #      2*****
            #      1*****
            #      0*****
            #       01234 Code11
            lat0 = code12*2.0 / 3.0 
            long0 = code34 + 100.0*z
            lat0 = lat0 + (code5*2.0/3.0) / 8.0
            long0 = long0 + code6 / 8.0
            lat0 = lat0 + ((code7-x+1)*2.0/3.0) / 8.0 / 10.0
            long0 = long0 + (code8+y) / 8.0 / 10.0
            lat0 = lat0 + (math.floor((code9-1)/2)+2*x-2) * 2.0 / 3.0 / 8.0 / 10.0 / 2.0
            long0 = long0 + ((code9-1)%2-2*y) / 8.0 / 10.0 / 2.0
            lat0 = lat0 + (codeex10-x+1)*2.0 / 3.0 / 8.0 / 10.0 / 2.0 / 5.0
            long0 = long0 + (codeex11+y) / 8.0 / 10.0 / 2.0 / 5.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/2.0/5.0
            dlong = 1.0/8.0/10.0/2.0/5.0
        else:
            # Extended 50m grid square code with 13 digits
            # code 11
            #     N
            #   3 | 4
            # W - + - E
            #   1 | 2
            #     S
            lat0 = code12 * 2 / 3  
            long0 = code34 + 100*z
            lat0 = lat0  + (code5 * 2 / 3) / 8
            long0 = long0 +  code6 / 8
            lat0 = lat0  + (code7 * 2 / 3) / 8 / 10
            long0 = long0 + code8 / 8 / 10
            lat0  = lat0  + ((codeex9-x+1) * 2 / 3) / 8 / 10 / 10
            long0 = long0 + (codeex10+y) / 8 / 10 / 10
            lat0  = lat0  + (math.floor((codeex11-1)/2)+x-1) * 2 / 3 / 8 / 10 / 10 / 2
            long0 = long0 + ((codeex11-1)%2-y) / 8 / 10 / 10 / 2
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2/3/8/10/10/2
            dlong = 1/8/10/10/2

        lat1  = "%.10f" % (lat0-dlat)
        long1 = "%.10f" % (long0+dlong)
        lat0  = "%.10f" % lat0
        long0 = "%.10f" % long0
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
    
    elif len(code)==14:
        if(not extension):
            pass
        elif spec == "ex50m_14":
            # Extended 50m grid square code with 14 digits
            # code 12
            #     N
            #   3 | 4
            # W - + - E
            #   1 | 2
            #     S
            lat0 = code12 * 2 / 3  
            long0 = code34 + 100*z
            lat0 = lat0  + (code5 * 2 / 3) / 8 
            long0 = long0 + code6 / 8
            lat0 = lat0  + ((code7-x+1) * 2 / 3) / 8 / 10
            long0 = long0 + (code8+y) / 8 / 10
            lat0 = lat0  + (math.floor((code9-1)/2)+2*x-2) * 2 / 3 / 8 / 10 / 2
            long0 = long0 + ((code9-1)%2-2*y) / 8 / 10 / 2
            lat0 = lat0  + ((codeex10-x+1) * 2 / 3) / 8 / 10 / 2 / 5
            long0 = long0 + (codeex11+y) / 8 / 10 / 2 / 5
            lat0 = lat0  + (math.floor((codeex12-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 5 / 2
            long0 = long0 + ((codeex12-1)%2-y) / 8 / 10 / 2 / 5 / 2
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2/3/8/10/2/5/2
            dlong = 1/8/10/2/5/2
        else:
            # Extended 10m grid square code (14 digits)
            #   Code11
            #      9**********
            #      8**********
            #      7**********
            #      6**********
            #      5**********
            #      4**********
            #      3**********
            #      2**********
            #      1**********
            #      0**********
            #      0123456789 Code12
            lat0  = code12 * 2.0 / 3.0  
            long0 = code34 + 100*z
            lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
            long0 = long0 +  code6 / 8.0
            lat0 = lat0  + code7 * 2.0 / 3.0 / 8.0 / 10.0
            long0 = long0 + code8 / 8.0 / 10.0
            lat0  = lat0  + (codeex9 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0
            long0 = long0 + codeex10 / 8.0 / 10.0 / 10.0
            lat0  = lat0  + ((codeex11-x+1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0
            long0 = long0 + (codeex12+y) / 8.0 / 10.0 / 10.0 / 10.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/10.0/10.0
            dlong = 1.0/8.0/10.0/10.0/10.0
        
        lat1  = "%.12f" % (lat0-dlat)
        long1 = "%.12f" % (long0+dlong)
        lat0  = "%.12f" % lat0
        long0 = "%.12f" % long0
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
 
    elif len(code)==16:
        if(not extension):
            pass
        else:
            # Extended 1m grid square code (16 digits)
            #   Code13
            #      9**********
            #      8**********
            #      7**********
            #      6**********
            #      5**********
            #      4**********
            #      3**********
            #      2**********
            #      1**********
            #      0**********
            #      0123456789 Code14
            lat0  = code12 * 2.0 / 3.0
            long0 = code34 + 100*z
            lat0  = lat0  + (code5 * 2.0 / 3.0) / 8.0
            long0 = long0 +  code6 / 8.0
            lat0  = lat0  + (code7 * 2.0 / 3.0) / 8.0 / 10.0
            long0 = long0 + code8 / 8.0 / 10.0
            lat0  = lat0  + (codeex9 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0
            long0 = long0 + codeex10 / 8.0 / 10.0 / 10.0
            lat0 = lat0  + (codeex11 * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0
            long0 = long0 + codeex12 / 8.0 / 10.0 / 10.0 / 10.0
            lat0 = lat0  + ((codeex13-x+1) * 2.0 / 3.0) / 8.0 / 10.0 / 10.0 / 10.0 / 10.0
            long0 = long0 + (codeex14+y) / 8.0 / 10.0 / 10.0 / 10.0 / 10.0
            lat0 = (1-2*x)*lat0
            long0 = (1-2*y)*long0
            dlat = 2.0/3.0/8.0/10.0/10.0/10.0/10.0
            dlong = 1.0/8.0/10.0/10.0/10.0/10.0

        lat1  = "%.14f" % (lat0-dlat)
        long1 = "%.14f" % (long0+dlong)
        lat0  = "%.14f" % lat0
        long0 = "%.14f" % long0
        xx = {"lat0":float(lat0), "long0":float(long0), "lat1":float(lat1), "long1":float(long1)}
        return xx
            
    xx = {"lat0":int("99999"), "long0":int("99999"), "lat1":int("99999"), "long1":int("99999")}
    return xx

# calculate 3rd mesh code
def cal_meshcode(latitude, longitude):
  return cal_meshcode3(latitude,longitude)

# calculate 1st mesh code
def cal_meshcode1(latitude, longitude):
  if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
    return "999999"
  if latitude < 0:
    o = 4
  else:
    o = 0
  if longitude < 0:
    o = o + 2
  if abs(longitude) >= 100:
    o = o + 1
  z = o % 2
  y = ((o - z)/2) % 2
  x = (o - 2*y - z)/4
#
  o = o + 1
#
  latitude = (1-2*x)*latitude
  longitude = (1-2*y)*longitude
  #
  p = math.floor(latitude*60/40)
  u = math.floor(longitude-100*z)
  #
  o = int(o)
  p = int(p)
  u = int(u)
  #
  if u < 10:
    if p < 10:
      mesh = str(o)+"00"+str(p)+"0"+str(u)
    else:
      if p < 100:
       mesh = str(o)+"0"+str(p)+"0"+str(u)
      else:
       mesh = str(o)+str(p)+"0"+str(u)
  else:
    if p < 10:
      mesh = str(o)+"00"+str(p)+str(u)
    else:
      if p < 100:
          mesh = str(o)+"0"+str(p)+str(u)
      else:
        mesh = str(o)+str(p)+str(u)
  return mesh

# calculate 2nd mesh code
def cal_meshcode2(latitude, longitude):
 if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
   return "99999999"
 if latitude < 0:
    o = 4
 else:
    o = 0
 if longitude < 0:
    o = o + 2
 if abs(longitude) >= 100: o = o + 1
 z = o % 2
 y = ((o - z)/2) % 2
 x = (o - 2*y - z)/4
#
 o = o + 1
#
 latitude = (1-2*x)*latitude
 longitude = (1-2*y)*longitude
#
 p = math.floor(latitude*60/40)
 a = (latitude*60/40-p)*40
 q = math.floor(a/5)
 u = math.floor(longitude-100*z)
 f = longitude-100*z-u
 v = math.floor(f*60/7.5)
 #
 o = int(o)
 p = int(p)
 u = int(u)
 q = int(q)
 v = int(v)
 #
 if u < 10:
    if p < 10:
      mesh = str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)
      else:
        mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)
 else:
    if p < 10:
      mesh =str(o)+"00"+str(p)+str(u)+str(q)+str(v)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)
      else:
        mesh = str(o)+str(p)+str(u)+str(q)+str(v)
 return mesh

# calculate 3rd mesh code
def cal_meshcode3(latitude, longitude):
  if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
    return "9999999999"
  if latitude < 0:
    o = 4
  else:
    o = 0
  if longitude < 0:
    o = o + 2
  if abs(longitude) >= 100:
      o = o + 1
  z = o % 2
  y = ((o - z)/2) % 2
  x = (o - 2*y - z)/4
#
  o = o + 1
#
  latitude = (1-2*x)*latitude
  longitude = (1-2*y)*longitude
  #
  p = math.floor(latitude*60/40)
  a = (latitude*60/40-p)*40
  q = math.floor(a/5)
  b = (a/5-q)*5
  r = math.floor(b*60/30)
  c = (b*60/30-r)*30
  u = math.floor(longitude-100*z)
  f = longitude-100*z-u
  v = math.floor(f*60/7.5)
  g = (f*60/7.5-v)*7.5
  w = math.floor(g*60/45)
  h = (g*60/45-w)*45
  #
  o = int(o)
  p = int(p)
  u = int(u)
  q = int(q)
  v = int(v)
  r = int(r)
  w = int(w)
  #
  if u < 10:
    if p < 10:
      mesh = str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)
      else:
        mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)
  else:
    if p < 10:
      mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)
    else:
      if(p < 100):
        mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)
      else:
        mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)
  return mesh

# calculate 4th mesh code
def cal_meshcode4(latitude, longitude):
  if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
    return "99999999999"
  if latitude < 0:
    o = 4
  else:
    o = 0
  if longitude < 0:
    o = o + 2
  if abs(longitude) >= 100:
    o = o + 1
  z = o % 2
  y = ((o - z)/2) % 2
  x = (o - 2*y - z)/4
#
  o = o + 1
#
  latitude = (1-2*x)*latitude
  longitude = (1-2*y)*longitude
  #
  p = math.floor(latitude*60/40)
  a = (latitude*60/40-p)*40
  q = math.floor(a/5)
  b = (a/5-q)*5
  r = math.floor(b*60/30)
  c = (b*60/30-r)*30
  s2u = math.floor(c/15)
  u = math.floor(longitude-100*z)
  f = longitude-100*z-u
  v = math.floor(f*60/7.5)
  g = (f*60/7.5-v)*7.5
  w = math.floor(g*60/45)
  h = (g*60/45-w)*45
  s2l = math.floor(h/22.5)
  s2 = s2u*2+s2l+1
  #
  o = int(o)
  p = int(p)
  u = int(u)
  q = int(q)
  v = int(v)
  r = int(r)
  w = int(w)
  s2 = int(s2)
  #
  if u < 10:
    if p < 10:
      mesh = str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
      else:
        mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
  else:
    if p < 10:
      mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
      else:
        mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)
  return mesh

# calculate 5rd mesh code
def cal_meshcode5(latitude, longitude):
  if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
    return "999999999999"
  if latitude < 0:
    o = 4
  else:
    o = 0
  if longitude < 0:
    o = o + 2
  if abs(longitude) >= 100:
    o = o + 1
  z = o % 2
  y = ((o - z)/2) % 2
  x = (o - 2*y - z)/4
#
  o = o + 1
#
  latitude = (1-2*x)*latitude
  longitude = (1-2*y)*longitude
  #
  p = math.floor(latitude*60/40)
  a = (latitude*60/40-p)*40
  q = math.floor(a/5)
  b = (a/5-q)*5
  r = math.floor(b*60/30)
  c = (b*60/30-r)*30
  s2u = math.floor(c/15)
  d = (c/15-s2u)*15
  s4u = math.floor(d/7.5)
  u = math.floor(longitude-100*z)
  f = longitude-100*z-u
  v = math.floor(f*60/7.5)
  g = (f*60/7.5-v)*7.5
  w = math.floor(g*60/45)
  h = (g*60/45-w)*45
  s2l =math.floor(h/22.5)
  i = (h/22.5-s2l)*22.5
  s4l = math.floor(i/11.25)
  s2 = s2u*2+s2l+1
  s4 = s4u*2+s4l+1
  #
  o = int(o)
  p = int(p)
  u = int(u)
  q = int(q)
  v = int(v)
  r = int(r)
  w = int(w)
  s2 = int(s2)
  s4 = int(s4)
  #
  if u < 10:
    if p < 10:
      mesh = str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
      else:
        mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
  else:
    if p < 10:
      mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
    else:
      if p < 100:
        mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
      else:
        mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)
  return mesh

# calculate 6rd mesh code
def cal_meshcode6(latitude, longitude):
  if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
      return "9999999999999"
  if latitude < 0:
    o = 4
  else:
    o = 0
  if longitude < 0:
    o = o + 2
  if abs(longitude) >= 100:
    o = o + 1
  z = o % 2
  y = ((o - z)/2) % 2
  x = (o - 2*y - z)/4
#
  o = o + 1
#
  latitude = (1-2*x)*latitude
  longitude = (1-2*y)*longitude
  #
  p = math.floor(latitude*60/40)
  a = latitude*60-p*40
  q = math.floor(math.floor(a/5))
  b = a-q*5
  r = math.floor(b*60/30)
  c = b*60-r*30
  s2u = math.floor(c/15)
  d = c-s2u*15
  s4u = math.floor(d/7.5)
  e = d-s4u*7.5
  s8u = math.floor(e/3.75)
  u = math.floor(longitude-100*z)
  f = longitude-100*z-u
  v = math.floor(f*60/7.5)
  g = f*60-v*7.5
  w = math.floor(g*60/45)
  h = g*60-w*45
  s2l = math.floor(h/22.5)
  i = h-s2l*22.5
  s4l = math.floor(i/11.25)
  j = i-s4l*11.25
  s8l = math.floor(j/5.625)
  s2 = s2u*2+s2l+1
  s4 = s4u*2+s4l+1
  s8 = s8u*2+s8l+1
  #
  o = int(o)
  p = int(p)
  u = int(u)
  q = int(q)
  v = int(v)
  r = int(r)
  w = int(w)
  s2 = int(s2)
  s4 = int(s4)
  s8 = int(s8)
  #
  if u < 10:
   if p < 10:
     mesh = str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
   else:
     if p < 100:
       mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
     else:
       mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
  else:
   if p < 10:
     mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
   else:
     if p < 100:
       mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
     else:
       mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(s4)+str(s8)
  return int(mesh)

# calculate an extended 100m (1km / 10) grid square code (12 digits) from a geographical position (latitude, longitude)
# - 3 arc-second for latitude and 4.5 arc-second for longitude
def cal_meshcode_ex100m_12(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s = math.floor(c/3)
    d = c-s*3
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    xx = math.floor(h/4.5)
    i = h-xx*4.5
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)
    return(mesh)

# calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude)
# - 3 arc-second for latitude and 4.5 arc-second for longitude
def cal_meshcode_ex100m_13(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "9999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s2u = math.floor(c/15)
    d = c-s2u*15
    et = math.floor(d/3)
    e = d-et*3
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    s2l = math.floor(h/22.5)
    i = h-s2l*22.5  
    jt = math.floor(i/4.5)
    j = i*0.5-jt*4.5
    s2 = s2u*2+s2l+1
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)
    return(mesh)

# calculate an extended 50m (100m / 2) grid square code (13 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude
def cal_meshcode_ex50m_13(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "9999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s = math.floor(c/3)
    d = c-s*3
    t2u = math.floor(d/1.5)
    e = d-t2u*1.5
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    xx = math.floor(h/4.5)
    i = h-xx*4.5
    t2l = math.floor(i/2.25)
    j = i-t2l*2.25
    tt = t2u*2+t2l+1
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(tt)
    return(mesh)

# calculate an extended 50m (100m / 2) grid square code (14 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude
def cal_meshcode_ex50m_14(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "9999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s2u = math.floor(c/15)
    d = c-s2u*15
    et = math.floor(d/3)
    e = d-et*3
    t2u = math.floor(e/1.5)
    f = e-t2u*1.5
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    s2l = math.floor(h/22.5)
    i = h-s2l*22.5  
    jt = math.floor(i/4.5)
    j = i-jt*4.5
    t2l = math.floor(j/2.25)
    k = j-t2l*2.25
    s2 = s2u*2+s2l+1
    t2 = t2u*2+t2l+1
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s2)+str(et)+str(jt)+str(t2)
    return(mesh)

# calculate an extended 10m (100m / 10) grid square code (14 digits) from a geographical position (latitude, longitude)
# - 0.3 arc-second for latitude and 0.45 arc-second for longitude
def cal_meshcode_ex10m_14(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "99999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s = math.floor(c/3)
    d = c-s*3
    t = math.floor(d/0.3)
    e = d-t*0.3
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    xx = math.floor(h/4.5)
    i = h-xx*4.5  
    yy = math.floor(i/0.45)
    j = i-yy*0.45
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)
    return(mesh)

# calculate an extended 1m (10m / 10) grid square code (16 digits) from a geographical position (latitude, longitude)
# - 0.03 arc-second for latitude and 0.045 arc-second for longitude
def cal_meshcode_ex1m_16(latitude, longitude):
    if latitude < -90 or latitude > 90 or longitude < -180 or longitude > 180:
        return "9999999999999999"
    if latitude < 0:
        o = 4
    else:
        o = 0
        
    if longitude < 0:
        o = o + 2
  
    if abs(longitude) >= 100:
        o = o + 1
      
    z = o % 2
    y = ((o - z)/2) % 2
    x = (o - 2*y - z)/4
    #
    o = o + 1
    #
    latitude = (1-2*x)*latitude
    longitude = (1-2*y)*longitude
    #
    p = math.floor(latitude*60/40)
    a = latitude*60-p*40
    q = math.floor(a/5)
    b = a-q*5
    r = math.floor(b*60/30)
    c = b*60-r*30
    s = math.floor(c/3)
    d = c-s*3
    t = math.floor(d/0.3)
    e = d-t*0.3
    tt = math.floor(e/0.03)
    ee = e-tt*0.03
    u = math.floor(longitude-100*z)
    f = longitude-100*z-u
    v = math.floor(f*60/7.5)
    g = f*60-v*7.5
    w = math.floor(g*60/45)
    h = g*60-w*45
    xx = math.floor(h/4.5)
    i = h-xx*4.5  
    yy = math.floor(i/0.45)
    j = i-yy*0.45
    zz = math.floor(j/0.045)
    k = j-zz*0.045
    #
    if u < 10:
        if p < 10:
            mesh =str(o)+"00"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
        else:
            if p < 100:
                mesh = str(o)+"0"+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
            else:
                mesh = str(o)+str(p)+"0"+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
    else:
        if p < 10:
            mesh = str(o)+"00"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
        else:
            if(p < 100):
                mesh = str(o)+"0"+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
            else:
                mesh = str(o)+str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(s)+str(xx)+str(t)+str(yy)+str(tt)+str(zz)
    return(mesh)


def Vincenty(latitude1,longitude1,latitude2,longitude2):
    # WGS84
    f = 1.0/298.257223563
    a = 6378137.0
    b = 6356752.314245
    #
    if (math.fabs(latitude1-latitude2)<1e-18) and (math.fabs(longitude1-longitude2)<1e-18):
        return 0.0 # two places are equivalent
    L = (longitude1 - longitude2)/180.0*math.pi
    U1 = math.atan((1.0-f)*math.tan(latitude1/180.0*math.pi))
    U2 = math.atan((1.0-f)*math.tan(latitude2/180.0*math.pi))
    vlambda = L
    dlambda = 10.0
    while abs(dlambda)>1e-12:
        cs = math.cos(U2)*math.sin(vlambda)
        cscc = math.cos(U1)*math.sin(U2)-math.sin(U1)*math.cos(U2)*math.cos(vlambda)
        sinsigma = math.sqrt(cs*cs + cscc*cscc)
        cossigma = math.sin(U1)*math.sin(U2)+math.cos(U1)*math.cos(U2)*math.cos(vlambda)
        sigma = math.atan(sinsigma/cossigma)
        sinalpha = math.cos(U1)*math.cos(U2)*math.sin(vlambda)/sinsigma
        cos2alpha = 1.0 - sinalpha*sinalpha
        if cos2alpha == 0.0:
            C = 0.0
            lambda0 = L + f*sinalpha*sigma
        else:
            cos2sigmam = cossigma - 2.0*math.sin(U1)*math.sin(U2)/cos2alpha
            C = f/16.0*cos2alpha*(4.0+f*(4-3*cos2alpha))
            lambda0 = L + (1.0-C)*f*sinalpha*(sigma + C*sinsigma*(cos2sigmam + C*cossigma*(-1+2*cos2sigmam*cos2sigmam)))
        dlambda = lambda0 - vlambda
        vlambda = lambda0
    if C == 0.0:
        A = 1.0
        dsigma = 0.0
    else:
        u2 = cos2alpha * (a*a-b*b)/(b*b)
        A = 1.0 + u2/16384.0*(4096.0 + u2 * (-768.0 + u2*(320.0-175.0*u2)))
        B = u2/1024.0*(256.0+u2*(-128.0+u2*(74.0-47.0*u2)))
        dsigma = B*sinsigma*(cos2sigmam + 1.0/4.0*B*(cossigma*(-1.0+2.0*cos2sigmam*cos2sigmam)-1.0/6.0*B*cos2sigmam*(-3.0+4.0*sinsigma*sinsigma)*(-3.0+4.0*cos2sigmam*cos2sigmam)))
    s = b*A*(sigma-dsigma)
    return s

def cal_area_from_meshcode(meshcode, extension=False,spec=""):
    latlong = meshcode_to_latlong_grid(meshcode, extension, spec)
    return cal_area_from_latlong(latlong)

def cal_area_from_latlong(latlong):
    W1 = Vincenty(latlong["lat0"],latlong["long0"],latlong["lat0"],latlong["long1"])
    W2 = Vincenty(latlong["lat1"],latlong["long0"],latlong["lat1"],latlong["long1"])
    H = Vincenty(latlong["lat0"],latlong["long0"],latlong["lat1"],latlong["long0"])
    A=(W1+W2)*H*0.5
    xx={"W1":W1,"W2":W2,"H":H,"A":A}
    return xx

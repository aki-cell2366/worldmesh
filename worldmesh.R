#
# Copyright (c) 2015-2022 Research Institute for World Grid Squares 
# Prof. Dr. Aki-Hiro Sato (Yokohama City University)
# All rights reserved. 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# R functions to calculate the world grid square code.
# The world grid square code computed by this library 
# is compatible to Japanese grid square code (JIS X 0410). 
#
# This library for World grid square code complies ISO 24108-1:
# https://www.iso.org/standard/87706.html
# 
# Release notes
# Version 1.0 : Released on 02 June 2015
# Version 1.1 : Released on 09 July 2015
# Version 1.2 : Released on 30 July 2015
# Version 1.3 : Released on 03 August 2015
# Version 1.3.1 : Released on 05 August 2015
# Version 1.4 : Released on 19 December 2015
# Version 1.5 : Released on 20 November 2019
# Version 1.6 : Released on 09 October 2020
# Version 1.7 : Released on 19 June 2021
# Version 1.75 : Released on 13 April 2022
# Version 1.8 : Released on 02 October 2022
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

meshcode_to_latlong <- function (meshcode, extension=F, spec=""){
  res<-meshcode_to_latlong_grid(meshcode, extension, spec)
  xx<-list(res$lat0,res$long0)
  names(xx)<-c("lat","long")
  return(xx)
}

meshcode_to_latlong_NW <- function (meshcode, extension=F, spec=""){
  res<-meshcode_to_latlong_grid(meshcode, extension, spec)
  xx<-list(res$lat0,res$long0)
  names(xx)<-c("lat","long")
  return(xx)
}

meshcode_to_latlong_SW <- function(meshcode, extension=F, spec=""){
  res<-meshcode_to_latlong_grid(meshcode, extension, spec)
  xx<-list(res$lat1,res$long0)
  names(xx)<-c("lat","long")
  return(xx)
}

meshcode_to_latlong_NE <- function(meshcode, extension=F, spec=""){
  res<-meshcode_to_latlong_grid(meshcode, extension, spec)
  xx<-list(res$lat0,res$long1)
  names(xx)<-c("lat","long")
  return(xx)
}

meshcode_to_latlong_SE <- function(meshcode, extension=F, spec=""){
  res<-meshcode_to_latlong_grid(meshcode, extension, spec)
  xx<-list(res$lat1,res$long1)
  names(xx)<-c("lat","long")
  return(xx)
}

meshcode_to_latlong_grid <- function (meshcode, extension=F, spec=""){
    # The extension option determines high resolution gird square codes if extension = T
    # 100m x 100m 3 arc-second x 4.5 arc second (12 digits) : default
    # 100m x 100m 3 arc-second x 4.5 arc second (13 digits) : spec = "ex100m_13"
    # 50m x 50m 1.5 arc-second x 2.25 arc second (13 digits) : default
    # 50m x 50m 1.5 arc-second x 2.25 arc second (14 digits) : sepc = "ex50m_14"
    # 10m x 10m 0.3 arc-second x 0.45 arc second (14 digits) : default
    # 1m x 1m 0.03 arc-second x 0.045 arc second (16 digits) : default
    code <- as.character(meshcode)
    if (length(grep("^[0-9]{6}", code)) == 1) { # more than 1st grid
        code0 <- as.numeric(substring(code, 1, 1)) # code0 : 1 to 8
        code0 <- code0 - 1 # transforming code0 from 0 to 7
          code12 <- substring(code, 2, 4)
        if(substring(code12,1,2)=="00"){
          code12 <- as.numeric(substring(code, 4, 4))
        }
        else{
          if(substring(code12,1,2)=="0"){
            code12 <- as.numeric(substring(code, 3, 4))
          }
          else{
            code12 <- as.numeric(substring(code, 2, 4))
          }
        }
        if(substring(code,5,5)=="0"){
          code34 <- as.numeric(substring(code, 6, 6))
        }
        else{
          code34 <- as.numeric(substring(code, 5, 6))
        }
    }
    else {
        return(NULL)
    }

    if (length(grep("^[0-9]{8}", code)) == 1) { # more than 2nd grid
        code5 <- as.numeric(substring(code, 7, 7))
        code6 <- as.numeric(substring(code, 8, 8))
    }

    if (length(grep("^[0-9]{10}", code)) == 1) { # more than 3rd grid
        code7 <- as.numeric(substring(code, 9, 9))
        code8 <- as.numeric(substring(code, 10, 10))
    }

    if(!extension){
      if (length(grep("^[0-9]{11}", code)) == 1) { # more than 4th grid
        code9 <- as.numeric(substring(code, 11, 11))
      }
      if (length(grep("^[0-9]{12}", code)) == 1) { # more than 5th grid
        code10 <- as.numeric(substring(code, 12, 12))
      }
      if (length(grep("^[0-9]{13}", code)) == 1) { # more than 6th grid
        code11 <- as.numeric(substring(code, 13, 13))
      }
    }else{
      if (length(grep("^[0-9]{12}", code)) == 1) { # Extended 100m grid square code (12digits)
          codeex9 <- as.numeric(substring(code, 11, 11))      
          codeex10 <- as.numeric(substring(code, 12, 12))
      }
      if (length(grep("^[0-9]{13}", code)) == 1) { # Extended 100m or 50m grid square code (13digits)
          code9 <- as.numeric(substring(code, 11, 11))
          codeex10 <- as.numeric(substring(code, 12, 12))      
          codeex11 <- as.numeric(substring(code, 13, 13))
      }
      if (length(grep("^[0-9]{14}", code)) == 1) { # Extended 50m or 10m grid square code (14digits)
          codeex11 <- as.numeric(substring(code, 13, 13))      
          codeex12 <- as.numeric(substring(code, 14, 14))
      }
      if (length(grep("^[0-9]{16}", code)) == 1) { # Extended 1m grid square code (16digits)
          codeex13 <- as.numeric(substring(code, 15, 15))      
          codeex14 <- as.numeric(substring(code, 16, 16))
      }
    }

    # 0'th grid
    z <- code0 %% 2
    y <- ((code0 - z)/2) %% 2
    x <- (code0 - 2*y - z)/4

    if(nchar(code)==6){  #  1st grid
      lat0  <- (code12-x+1) * 2 / 3;        
      long0 <- (code34+y) + 100*z;
      lat0  <- (1-2*x)*lat0;        
      long0 <- (1-2*y)*long0;
      dlat <- 2/3;
      dlong <- 1;
      lat1  <- sprintf("%.8f", lat0-dlat);  
      long1 <- sprintf("%.8f", long0+dlong);
      lat0  <- sprintf("%.8f", lat0);  
      long0 <- sprintf("%.8f", long0);
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1));
      names(xx) <- c("lat0","long0","lat1","long1");
      return(xx)
    }
    if(nchar(code)==8) {  # 2nd grid
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
      lat0  <- code12 * 2 / 3;  
      long0 <- code34 + 100*z;
      lat0  <- lat0  + ((code5-x+1) * 2 / 3) / 8; 
      long0 <- long0 +  (code6+y) / 8;
      lat0 <- (1-2*x) * lat0;
      long0 <- (1-2*y) * long0;
      dlat <- 2/3/8;
      dlong <- 1/8;
      lat1  <- sprintf("%.8f", lat0-dlat);  
      long1 <- sprintf("%.8f", long0+dlong);
      lat0  <- sprintf("%.8f", lat0);  
      long0 <- sprintf("%.8f", long0);
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    if(nchar(code)==10) {  # 3rd grid
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
      lat0  <- code12 * 2 / 3;  
      long0 <- code34 + 100*z;
      lat0  <- lat0  + (code5 * 2 / 3) / 8; 
      long0 <- long0 +  code6 / 8;
      lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
      long0 <- long0 +  (code8+y) / 8 / 10;
      lat0 <- (1-2*x)*lat0;
      long0 <- (1-2*y)*long0;
      dlat <- 2/3/8/10;
      dlong <- 1/8/10;
      lat1  <- sprintf("%.8f", lat0-dlat); 
      long1 <- sprintf("%.8f", long0+dlong);
      lat0  <- sprintf("%.8f", lat0); 
      long0 <- sprintf("%.8f", long0);
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    # code 9
    #     N
    #   3 | 4
    # W - + - E
    #   1 | 2
    #     S
    if(nchar(code)==11) {  # 4th grid (11 digits)
      lat0  <- code12 * 2 / 3  
      long0 <- code34 + 100*z
      lat0  <- lat0  + (code5 * 2 / 3) / 8 
      long0 <- long0 +  code6 / 8
      lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10
      long0 <- long0 +  (code8+y) / 8 / 10
      lat0  <- lat0  + (floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2
      long0 <- long0 + ((code9-1)%%2-y) / 8 / 10 / 2
      lat0 <- (1-2*x)*lat0
      long0 <- (1-2*y)*long0
      dlat <- 2/3/8/10/2
      dlong <- 1/8/10/2
      lat1  <- sprintf("%.8f", lat0-dlat) 
      long1 <- sprintf("%.8f", long0+dlong)
      lat0  <- sprintf("%.8f", lat0)
      long0 <- sprintf("%.8f", long0)
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    if(nchar(code)==12) {  # 5th grid (12 digits)

      if(!extension){
        # code 10
        #     N
        #   3 | 4
        # W - + - E
        #   1 | 2
        #     S
        lat0  <- code12 * 2 / 3;  
        long0 <- code34 + 100*z;
        lat0  <- lat0  + (code5 * 2 / 3) / 8; 
        long0 <- long0 +  code6 / 8;
        lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
        long0 <- long0 +  (code8+y) / 8 / 10;
        lat0  <- lat0  + (floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2;
        long0 <- long0 + ((code9-1)%%2-y) / 8 / 10 / 2;
        lat0  <- lat0  + (floor((code10-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2;
        long0 <- long0 + ((code10-1)%%2-y) / 8 / 10 / 2 / 2;
        lat0 <- (1-2*x)*lat0;
        long0 <- (1-2*y)*long0;
        dlat <- 2/3/8/10/2/2;
        dlong <- 1/8/10/2/2;
      }else{ # Extended 100m grid square code
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
        lat0  <- code12 * 2 / 3  
        long0 <- code34 + 100*z
        lat0  <- lat0  + (code5 * 2 / 3) / 8
        long0 <- long0 +  code6 / 8
        lat0  <- lat0  + (code7 * 2 / 3) / 8 / 10
        long0 <- long0 + code8 / 8 / 10
        lat0  <- lat0  + ((codeex9-x+1) * 2 / 3) / 8 / 10 / 10
        long0 <- long0 + (codeex10+y) / 8 / 10 / 10
        lat0 <- (1-2*x)*lat0;
        long0 <- (1-2*y)*long0;
        dlat <- 2/3/8/10/10;
        dlong <- 1/8/10/10;
      }
      lat1  <- sprintf("%.10f", lat0-dlat); 
      long1 <- sprintf("%.10f", long0+dlong);
      lat0  <- sprintf("%.10f", lat0); 
      long0 <- sprintf("%.10f", long0);
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    if(nchar(code)==13) {  # 6rd grid (13 digits)
      if(!extension){
        # code 11
        #     N
        #   3 | 4
        # W - + - E
        #   1 | 2
        #     S
        lat0  <- code12 * 2 / 3;  
        long0 <- code34 + 100*z;
        lat0  <- lat0  + (code5 * 2 / 3) / 8; 
        long0 <- long0 +  code6 / 8;
        lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
        long0 <- long0 + (code8+y) / 8 / 10;
        lat0  <- lat0  + (floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2;
        long0 <- long0 + ((code9-1)%%2-y) / 8 / 10 / 2;
        lat0  <- lat0  + (floor((code10-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2;
        long0 <- long0 + ((code10-1)%%2-y) / 8 / 10 / 2 / 2;
        lat0  <- lat0  + (floor((code11-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2 / 2;
        long0 <- long0 + ((code11-1)%%2-y) / 8 / 10 / 2 / 2 / 2;
        lat0 <- (1-2*x)*lat0;
        long0 <- (1-2*y)*long0;
        dlat <- 2/3/8/10/2/2/2;
        dlong <- 1/8/10/2/2/2;
      }else{
	if(spec == "ex100m_13"){
        # Extended 100m grid square code with 13 digits
        #   Code10
        #      4*****
        #      3*****
        #      2*****
        #      1*****
        #      0*****
        #      01234 Code11
           lat0  <- code12 * 2 / 3;  
           long0 <- code34 + 100*z;
           lat0  <- lat0  + (code5 * 2 / 3) / 8; 
           long0 <- long0 + code6 / 8;
           lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
           long0 <- long0 + (code8+y) / 8 / 10;
           lat0  <- lat0  + (floor((code9-1)/2)+2*x-2) * 2 / 3 / 8 / 10 / 2;
           long0 <- long0 + ((code9-1)%%2-2*y) / 8 / 10 / 2;
           lat0  <- lat0  + ((codeex10-x+1) * 2 / 3) / 8 / 10 / 2 / 5
           long0 <- long0 + (codeex11+y) / 8 / 10 / 2 / 5
           lat0 <- (1-2*x)*lat0;
           long0 <- (1-2*y)*long0;
           dlat <- 2/3/8/10/2/5;
           dlong <- 1/8/10/2/5;
	}else{
        # Extended 50m grid square code with 13 digits
        # code 11
        #     N
        #   3 | 4
        # W - + - E
        #   1 | 2
        #     S
          lat0  <- code12 * 2 / 3  
          long0 <- code34 + 100*z
          lat0  <- lat0  + (code5 * 2 / 3) / 8
          long0 <- long0 +  code6 / 8
          lat0  <- lat0  + (code7 * 2 / 3) / 8 / 10
          long0 <- long0 + code8 / 8 / 10
          lat0  <- lat0  + ((codeex9-x+1) * 2 / 3) / 8 / 10 / 10
          long0 <- long0 + (codeex10+y) / 8 / 10 / 10
          lat0  <- lat0  + (floor((codeex11-1)/2)+x-1) * 2 / 3 / 8 / 10 / 10 / 2
          long0 <- long0 + ((codeex11-1)%%2-y) / 8 / 10 / 10 / 2
          lat0 <- (1-2*x)*lat0
          long0 <- (1-2*y)*long0
          dlat <- 2/3/8/10/10/2
          dlong <- 1/8/10/10/2
	}
      }
      lat1  <- sprintf("%.10f", lat0-dlat); 
      long1 <- sprintf("%.10f", long0+dlong);
      lat0  <- sprintf("%.10f", lat0); 
      long0 <- sprintf("%.10f", long0);
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    if(nchar(code)==14) {  # (14 digits)
      if(!extension){
      }else{ 
        if(spec == "ex50m_14"){
          # Extended 50m grid square code with 14 digits
          # code 12
          #     N
          #   3 | 4
          # W - + - E
          #   1 | 2
          #     S
          lat0  <- code12 * 2 / 3  
          long0 <- code34 + 100*z
          lat0  <- lat0  + (code5 * 2 / 3) / 8 
          long0 <- long0 + code6 / 8
          lat0  <- lat0  + ((code7-x+1) * 2 / 3) / 8 / 10
          long0 <- long0 + (code8+y) / 8 / 10
          lat0  <- lat0  + (floor((code9-1)/2)+2*x-2) * 2 / 3 / 8 / 10 / 2
          long0 <- long0 + ((code9-1)%%2-2*y) / 8 / 10 / 2
          lat0  <- lat0  + ((codeex10-x+1) * 2 / 3) / 8 / 10 / 2 / 5
          long0 <- long0 + (codeex11+y) / 8 / 10 / 2 / 5
          lat0  <- lat0  + (floor((codeex12-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 5 / 2
          long0 <- long0 + ((codeex12-1)%%2-y) / 8 / 10 / 2 / 5 / 2
          lat0 <- (1-2*x)*lat0
          long0 <- (1-2*y)*long0
          dlat <- 2/3/8/10/2/5/2
          dlong <- 1/8/10/2/5/2
        }else{
        # Extended 10m grid square code with 14 digits
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
          lat0  <- code12 * 2 / 3  
          long0 <- code34 + 100*z
          lat0  <- lat0  + (code5 * 2 / 3) / 8 
          long0 <- long0 +  code6 / 8
          lat0  <- lat0  + code7 * 2 / 3 / 8 / 10
          long0 <- long0 + code8 / 8 / 10
          lat0  <- lat0  + (codeex9 * 2 / 3) / 8 / 10 / 10
          long0 <- long0 + codeex10 / 8 / 10 / 10
          lat0  <- lat0  + ((codeex11-x+1) * 2 / 3) / 8 / 10 / 10 / 10
          long0 <- long0 + (codeex12+y) / 8 / 10 / 10 / 10
          lat0 <- (1-2*x)*lat0
          long0 <- (1-2*y)*long0
          dlat <- 2/3/8/10/10/10
          dlong <- 1/8/10/10/10
	}
      }
      lat1  <- sprintf("%.12f", lat0-dlat) 
      long1 <- sprintf("%.12f", long0+dlong)
      lat0  <- sprintf("%.12f", lat0)
      long0 <- sprintf("%.12f", long0)
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    if(nchar(code)==16) {  # (16 digits)
      if(!extension){
      }else{ # Extended 1m grid square code
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
        lat0  <- code12 * 2 / 3
        long0 <- code34 + 100*z
        lat0  <- lat0  + (code5 * 2 / 3) / 8
        long0 <- long0 +  code6 / 8
        lat0  <- lat0  + (code7 * 2 / 3) / 8 / 10
        long0 <- long0 + code8 / 8 / 10
        lat0  <- lat0  + (codeex9 * 2 / 3) / 8 / 10 / 10
        long0 <- long0 + codeex10 / 8 / 10 / 10
        lat0  <- lat0  + (codeex11 * 2 / 3) / 8 / 10 / 10 / 10
        long0 <- long0 + codeex12 / 8 / 10 / 10 / 10
        lat0  <- lat0  + ((codeex13-x+1) * 2 / 3) / 8 / 10 / 10 / 10 / 10
        long0 <- long0 + (codeex14+y) / 8 / 10 / 10 / 10 / 10
        lat0 <- (1-2*x)*lat0
        long0 <- (1-2*y)*long0
        dlat <- 2/3/8/10/10/10/10
        dlong <- 1/8/10/10/10/10
      }
      lat1  <- sprintf("%.14f", lat0-dlat) 
      long1 <- sprintf("%.14f", long0+dlong)
      lat0  <- sprintf("%.14f", lat0)
      long0 <- sprintf("%.14f", long0)
      xx <- list(as.numeric(lat0), as.numeric(long0), as.numeric(lat1), as.numeric(long1))
      names(xx) <- c("lat0","long0","lat1","long1")
      return(xx)
    }
    xx <- list(as.numeric(99999), as.numeric(99999),as.numeric(99999), as.numeric(99999))
    names(xx) <- c("lat0","long0","lat1","long1")
    return(xx)
}

# calculate 3rd mesh code
cal_meshcode <- function(latitude, longitude){
  return(cal_meshcode3(latitude,longitude))
}

# calculate 1st mesh code
cal_meshcode1 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  u <- floor(longitude-100*z)
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,sep="")
      }
      else{
        mesh <- paste(o,p,u,sep="")
      }
    }
  }
  return(mesh)
}

# calculate 2nd mesh code
cal_meshcode2 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("99999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- (latitude*60/40-p)*40
  q <- floor(a/5)
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,sep="")
      }
    }
  }
  return(mesh)
}

# calculate 3rd mesh code
cal_meshcode3 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("9999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- (latitude*60/40-p)*40
  q <- floor(a/5)
  b <- (a/5-q)*5
  r <- floor(b*60/30)
  c <- (b*60/30-r)*30
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- (f*60/7.5-v)*7.5
  w <- floor(g*60/45)
  h <- (g*60/45-w)*45
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,sep="")
      }
    }
  }
  return(mesh)
}

# calculate 4th mesh code
cal_meshcode4 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("99999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- (latitude*60/40-p)*40
  q <- floor(a/5)
  b <- (a/5-q)*5
  r <- floor(b*60/30)
  c <- (b*60/30-r)*30
  s2u <- floor(c/15)
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- (f*60/7.5-v)*7.5
  w <- floor(g*60/45)
  h <- (g*60/45-w)*45
  s2l <- floor(h/22.5)
  s2 <- s2u*2+s2l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s2,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s2,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s2,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s2,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s2,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s2,sep="")
      }
    }
  }
  return(mesh)
}

# calculate 5rd mesh code
cal_meshcode5 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- (latitude*60/40-p)*40
  q <- floor(a/5)
  b <- (a/5-q)*5
  r <- floor(b*60/30)
  c <- (b*60/30-r)*30
  s2u <- floor(c/15)
  d <- (c/15-s2u)*15
  s4u <- floor(d/7.5)
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- (f*60/7.5-v)*7.5
  w <- floor(g*60/45)
  h <- (g*60/45-w)*45
  s2l <- floor(h/22.5)
  i <- (h/22.5-s2l)*22.5
  s4l <- floor(i/11.25)
  s2 <- s2u*2+s2l+1
  s4 <- s4u*2+s4l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s2,s4,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s2,s4,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s2,s4,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s2,s4,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s2,s4,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s2,s4,sep="")
      }
    }
  }
  return(mesh)
}

# calculate 6rd mesh code
cal_meshcode6 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("9999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s2u <- floor(c/15)
  d <- c-s2u*15
  s4u <- floor(d/7.5)
  e <- d-s4u*7.5
  s8u <- floor(e/3.75)
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  s2l <- floor(h/22.5)
  i <- h-s2l*22.5
  s4l <- floor(i/11.25)
  j <- i-s4l*11.25
  s8l <- floor(j/5.625)
  s2 <- s2u*2+s2l+1
  s4 <- s4u*2+s4l+1
  s8 <- s8u*2+s8l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s2,s4,s8,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s2,s4,s8,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s2,s4,s8,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s2,s4,s8,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s2,s4,s8,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s2,s4,s8,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 100m (1km / 10) grid square code (12 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude
cal_meshcode_ex100m_12 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s <- floor(c/3)
  d <- c-s*3
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  xx <- floor(h/4.5)
  i <- h-xx*4.5
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s,xx,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s,xx,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s,xx,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s,xx,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s,xx,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s,xx,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude
cal_meshcode_ex100m_13 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("9999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s2u <- floor(c/15)
  d <- c-s2u*15
  et <- floor(d/3)
  e <- d-et*3
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  s2l <- floor(h/22.5)
  i <- h-s2l*22.5
  jt <- floor(i/4.5)
  j <- i*0.5-jt*4.5
  s2 <- s2u*2+s2l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s2,et,jt,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s2,et,jt,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s2,et,jt,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s2,et,jt,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s2,et,jt,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s2,et,jt,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 50m (100m / 2) grid square code (13 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude
cal_meshcode_ex50m_13 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("9999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s <- floor(c/3)
  d <- c-s*3
  t2u <- floor(d/1.5)
  e <- d-t2u*1.5
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  xx <- floor(h/4.5)
  i <- h-xx*4.5
  t2l <- floor(i/2.25)
  j <- i-t2l*2.25
  tt <- t2u*2+t2l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s,xx,tt,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s,xx,tt,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s,xx,tt,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s,xx,tt,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s,xx,tt,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s,xx,tt,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 50m (100m / 2) grid square code (14 digits) from a geographical position (latitude, longitude) - 1.5 arc-second for latitude and 2.25 arc-second for longitude
cal_meshcode_ex50m_14 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("99999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s2u <- floor(c/15)
  d <- c-s2u*15
  et <- floor(d/3)
  e <- d-et*3
  t2u <- floor(e/1.5)
  f <- e-t2u*1.5
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  s2l <- floor(h/22.5)
  i <- h-s2l*22.5
  jt <- floor(i/4.5)
  j <- i-jt*4.5
  t2l <- floor(j/2.25)
  k <- j-t2l*2.25
  s2 <- s2u*2+s2l+1
  t2 <- t2u*2+t2l+1
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s2,et,jt,t2,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s2,et,jt,t2,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s2,et,jt,t2,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s2,et,jt,t2,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s2,et,jt,t2,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s2,et,jt,t2,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) from a geographical position (latitude, longitude) - 0.3 arc-second for latitude and 0.45 arc-second for longitude
cal_meshcode_ex10m_14 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("99999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s <- floor(c/3)
  d <- c-s*3
  t <- floor(d/0.3)
  e <- d-t*0.3
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  xx <- floor(h/4.5)
  i <- h-xx*4.5
  yy <- floor(i/0.45)
  j <- i-yy*0.45
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s,xx,t,yy,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s,xx,t,yy,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s,xx,t,yy,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s,xx,t,yy,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s,xx,t,yy,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s,xx,t,yy,sep="")
      }
    }
  }
  return(mesh)
}
# calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) from a geographical position (latitude, longitude) - 0.03 arc-second for latitude and 0.045 arc-second for longitude
cal_meshcode_ex1m_16 <- function(latitude, longitude){
  if(latitude < -90 || latitude > 90 || longitude < -180 || longitude > 180){
    return("9999999999999999")
  }
  if(latitude < 0){
    o <- 4
  }
  else{
    o <- 0
  }
  if(longitude < 0){
    o <- o + 2
  }
  if(abs(longitude) >= 100) o <- o + 1
  z <- o %% 2
  y <- ((o - z)/2) %% 2
  x <- (o - 2*y - z)/4
#
  o <- o + 1
#
  latitude <- (1-2*x)*latitude
  longitude <- (1-2*y)*longitude
  #
  p <- floor(latitude*60/40)
  a <- latitude*60-p*40
  q <- floor(a/5)
  b <- a-q*5
  r <- floor(b*60/30)
  c <- b*60-r*30
  s <- floor(c/3)
  d <- c-s*3
  t <- floor(d/0.3)
  e <- d-t*0.3
  tt <- floor(e/0.03)
  ee <- e-tt*0.03
  u <- floor(longitude-100*z)
  f <- longitude-100*z-u
  v <- floor(f*60/7.5)
  g <- f*60-v*7.5
  w <- floor(g*60/45)
  h <- g*60-w*45
  xx <- floor(h/4.5)
  i <- h-xx*4.5
  yy <- floor(i/0.45)
  j <- i-yy*0.45
  zz <- floor(j/0.045)
  k <- j-zz*0.045
  #
  if(u < 10){
    if(p < 10){
      mesh <- paste(o,"00",p,"0",u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,"0",u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
      }
      else{
        mesh <- paste(o,p,"0",u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
      }
    }
  }
  else{
    if(p < 10){
      mesh <- paste(o,"00",p,u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
    }else{
      if(p < 100){
        mesh <- paste(o,"0",p,u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
      }
      else{
        mesh <- paste(o,p,u,q,v,r,w,s,xx,t,yy,tt,zz,sep="")
      }
    }
  }
  return(mesh)
}

Vincenty <- function(latitude1,longitude1,latitude2,longitude2){
   # WGS84
   f = 1/298.257223563
   a = 6378137.0
   b = 6356752.314245
   if(abs(latitude1-latitude2)<1e-18 && abs(longitude1-longitude2)<1e-18){ # two places are equivalent
     return(0.0)
   }
   #
   L = (longitude1 - longitude2)/180*pi
   U1 = atan((1-f)*tan(latitude1/180*pi))
   U2 = atan((1-f)*tan(latitude2/180*pi))
   lambda = L
   dlambda = 10
   while(abs(dlambda) > 1e-12){
      cs = cos(U2)*sin(lambda)
      cscc = cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)
      sinsigma = sqrt(cs*cs + cscc*cscc)
      cossigma = sin(U1)*sin(U2)+cos(U1)*cos(U2)*cos(lambda)
      sigma = atan(sinsigma/cossigma)
      sinalpha = cos(U1)*cos(U2)*sin(lambda)/sinsigma
      cos2alpha = 1 - sinalpha*sinalpha
      if(cos2alpha == 0.0){
         C = 0.0
         lambda0 = L + f*sinalpha*sigma
      }else{
         cos2sigmam = cossigma - 2*sin(U1)*sin(U2)/cos2alpha
         C = f/16*cos2alpha*(4+f*(4-3*cos2alpha))
         lambda0 = L + (1-C)*f*sinalpha*(sigma + C*sinsigma*(cos2sigmam + C*cossigma*(-1+2*cos2sigmam*cos2sigmam)))
      }
      dlambda = lambda0 - lambda
      lambda = lambda0
   }
   if(C == 0.0){
      A = 1.0
      dsigma = 0.0
   }else{
      u2 = cos2alpha * (a*a-b*b)/(b*b)
      A = 1 + u2/16384*(4096 + u2 * (-768 + u2*(320-175*u2)))
      B = u2/1024*(256+u2*(-128+u2*(74-47*u2)))
      dsigma = B*sinsigma*(cos2sigmam + 1/4*B*(cossigma*(-1+2*cos2sigmam*cos2sigmam)-1/6*B*cos2sigmam*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmam*cos2sigmam)))
   }
   s = b*A*(sigma-dsigma)
   return(s)
}

cal_area_from_meshcode <- function(meshcode,extension=F,spec=""){
   latlong <- meshcode_to_latlong_grid(meshcode,extension,spec)
   return(cal_area_from_latlong(latlong))
}

cal_area_from_latlong <- function(latlong){
   W1 = Vincenty(latlong$lat0,latlong$long0,latlong$lat0,latlong$long1)
   W2 = Vincenty(latlong$lat1,latlong$long0,latlong$lat1,latlong$long1)
   H = Vincenty(latlong$lat0,latlong$long0,latlong$lat1,latlong$long0)
   A=(W1+W2)*H*0.5
   xx<-list(W1,W2,H,A)
   names(xx)<-c("W1","W2","H","A")
   return(xx)
}

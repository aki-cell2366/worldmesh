<?php
#
# Copyright (c) 2015-2022 Research Institute for World Grid Squares 
# Prof. Dr. Aki-Hiro Sato (Yokohama City University)
# All rights reserved. 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# PHP fnctions to calculate the world grid square code.
# The world grid square code computed by this library 
# is compatible to Japanese grid square code (JIS X 0410). 
#
# This library for World grid square code complies ISO 24108-1:
# https://www.iso.org/standard/87706.html
#
# Version 1.0: Released on 07 February 2017
# Version 1.6: Released on 09 October 2020
# Version 1.7: Released on 10 October 2021
# Version 1.75: Released on 13 April 2022
# Version 1.76: Released on 08 February 2025
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
# meshcode_to_latlong($meshcode, $extension=false)
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NW($meshcode, $extension=false)
# : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SW($meshcode, $extension=false)
# : calculate sourthern western geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_NE($meshcode, $extension=false)
# : calculate northern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_SE($meshcode, $extension=false)
# : calculate sourthern eastern geographic position of the grid (latitude, longitude) from meshcode
# meshcode_to_latlong_grid($meshcode, $extension=false)
# : calculate northern western and sourthern eastern geographic positions of the grid (latitude0, longitude0, latitude1, longitude1) from meshcode
# cal_meshcode(latitude,longitude)
#
# 2.
#
# : calculate a basic (1km) grid square code (10 digits) from a geographical position (latitude, longitude)
# cal_meshcode1($latitude,$longitude)
# : calculate an 80km grid square code (6 digits) from a geographical position (latitude, longitude)
# cal_meshcode2($latitude,$longitude)
# : calculate a 10km grid square code (8 digits) from a geographical position (latitude, longitude)
# cal_meshcode3($latitude,$longitude)
# : calculate a 1km grid square code (10 digits) from a geographical position (latitude, longitude)
# cal_meshcode4($latitude,$longitude)
# : calculate a 500m grid square code (11 digits) from a geographical position (latitude, longitude)
# cal_meshcode5($latitude,$longitude)
# : calculate a 250m grid square code (12 digits) from a geographical position (latitude, longitude)
# cal_meshcode6($latitude,$longitude)
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
# Vincenty($latitude1, $longitude1, $latitude2, $longitude)
# : calculate geodesitc distance between two points ($latitude1, $longitude1) and ($latitude2, $longitude2) placed on the WGS84 Earth ellipsoid based on the Vincenty's formulae (1975)
# cal_area_from_meshcode($meshcode,$extension=false)
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoide A) of world grid square indicated by meshcode
# cal_area_from_latlong($latlong)
# : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoid A) of a trapezoid on the WGS84 Earth ellipoid indicated by ($latlong['lat0'], $latlong['long0'], $latlong['lat1'], $latlong['long1'])
#

function meshcode_to_latlong_grid($meshcode,$extension=false){
    # The extension option determines high resolution gird square code for
    # 12 digits (100m x 100m grid squares 3 arc-seconds for latitude, 4.5 arc-seconds for longitude)
    # 13 digits (100m x 100m grid squares 3 arc-seconds for latitude, 4.5 arc-seconds for longitude)
    # 14 digits (10m x 10m grid squares 0.3 arc-seconds for latitude, 0.45 arc-seconds for longitude)
    # 16 digits (1m x 1m grid squares 0.03 arc-seconds for latitude, 0.045 arc-seconds for longitude)    
    #
   $code = strval($meshcode);
   $ncode = strlen($code);
   if ($ncode >= 6) { # more than 1st grid
       $code0 = intval(substr($code, 0, 1)); # code0 : 1 to 8
       $code0 = $code0 - 1; # transforming code0 from 0 to 7
       $code12 = substr($code, 1, 3);
       if(strcmp(substr($code12,0,2),"00")==0){
         $code12 = intval(substr($code, 3, 1));
       }
       else{
         if(strcmp(substr($code12,0,1),"0")==0){
           $code12 = intval(substr($code, 2, 2));
         }
         else{
           $code12 = intval(substr($code, 1, 3));
         }
       }
       if(substr($code,4,1)=="0"){
         $code34 = intval(substr($code, 5, 1));
       }
       else{
         $code34 = intval(substr($code, 4, 2));
       }
   }
   else {
       return(NULL);
   }

   if ($ncode >= 8) { # more than 2nd grid
       $code5 = intval(substr($code, 6, 1));
       $code6 = intval(substr($code, 7, 1));
   }

   if ($ncode >= 10) { # more than 3rd grid
       $code7 = intval(substr($code, 8, 1));
       $code8 = intval(substr($code, 9, 1));
   }
   if(!$extension){
     if ($ncode >= 11) { # more than 4th grid
         $code9 = intval(substr($code, 10, 1));
     }
     if ($ncode >= 12) { # more than 5th grid
         $code10 = intval(substr($code, 11, 1));
     }
     if ($ncode >= 13) { # more than 6th grid
         $code11 = intval(substr($code, 12, 1));
     }
   }else{
     if ($ncode == 12) { # Extended 100m grid square code (12 digits)
          $codeex9 = intval(substr($code, 10, 1));
          $codeex10 = intval(substr($code, 11, 1));
     }
     if ($ncode == 13) { # Extended 100m grid square code (13 digits)
          $code9 = intval(substr($code, 10, 1));     
          $codeex10 = intval(substr($code, 11, 1));
          $codeex11 = intval(substr($code, 12, 1));
     }
     if ($ncode == 14) { # Extended 10m grid square code (14 digits)
          $codeex9 = intval(substr($code, 10, 1));
          $codeex10 = intval(substr($code, 11, 1));
          $codeex11 = intval(substr($code, 12, 1));
          $codeex12 = intval(substr($code, 13, 1));
     }
     if ($ncode == 16) { # Extended 1m grid square code (16 digits)
          $codeex9 = intval(substr($code, 10, 1));
          $codeex10 = intval(substr($code, 11, 1));
          $codeex11 = intval(substr($code, 12, 1));
          $codeex12 = intval(substr($code, 13, 1));
          $codeex13 = intval(substr($code, 14, 1));
          $codeex14 = intval(substr($code, 15, 1));
     }
   }

   # 0'th grid
   $z = $code0 % 2;
   $y = (($code0 - $z)/2) % 2;
   $x = ($code0 - 2*$y - $z)/4;

   if($ncode==6){ # 1st grid
     $lat0 = ($code12-$x+1) * 2 / 3;
     $long0 = ($code34+$y) + 100*$z;
     $lat0 = (1-2*$x)*$lat0;
     $long0 = (1-2*$y)*$long0;
     $dlat = 2/3;
     $dlong = 1;
     $lat1 = sprintf("%.8f", $lat0-$dlat);
     $long1 = sprintf("%.8f", $long0+$dlong);
     $lat0 = sprintf("%.8f", $lat0);
     $long0 = sprintf("%.8f", $long0);
     $xx = array("lat0" => floatval($lat0), "long0" => floatval($long0), "lat1" => floatval($lat1), "long1" => floatval($long1));
     return($xx);
   }
   if($ncode==8) { # 2nd grid
     $lat0 = $code12 * 2 / 3;
     $long0 = $code34 + 100*$z;
     $lat0 = $lat0  + (($code5-$x+1) * 2 / 3) / 8;
     $long0 = $long0 +  ($code6+$y) / 8;
     $lat0 = (1-2*$x) * $lat0;
     $long0 = (1-2*$y) * $long0;
     $dlat = 2/3/8;
     $dlong = 1/8;
     $lat1 = sprintf("%.8f", $lat0-$dlat);
     $long1 = sprintf("%.8f", $long0+$dlong);
     $lat0 = sprintf("%.8f", $lat0);
     $long0 = sprintf("%.8f", $long0);
     $xx = array("lat0" => floatval($lat0), "long0" => floatval($long0), "lat1" => floatval($lat1), "long1" => floatval($long1));
     return($xx);
   }
   if($ncode==10) { # 3rd grid
#    echo $code12;
#    echo $code34;
     $lat0 = $code12 * 2 / 3;
     $long0 = $code34 + 100*$z;
     $lat0 = $lat0  + ($code5 * 2 / 3) / 8;
     $long0 = $long0 +  $code6 / 8;
     $lat0 = $lat0  + (($code7-$x+1) * 2 / 3) / 8 / 10;
     $long0 = $long0 + ($code8+$y) / 8 / 10;
     $lat0 = (1-2*$x)*$lat0;
     $long0 = (1-2*$y)*$long0;
     $dlat = 2/3/8/10;
     $dlong = 1/8/10;
     $lat1 = sprintf("%.8f", $lat0-$dlat);
     $long1 = sprintf("%.8f", $long0+$dlong);
     $lat0 = sprintf("%.8f", $lat0);
     $long0 = sprintf("%.8f", $long0);
     $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
     return($xx);
   }
   # code 9
   #     N
   #   3 | 4
   # W - + - E
   #   1 | 2
   #     S
   if($ncode==11) { # 4th grid (11 digits)
     $lat0 = $code12 * 2 / 3;
     $long0 = $code34 + 100*$z;
     $lat0 = $lat0  + ($code5 * 2 / 3) / 8;
     $long0 = $long0 +  $code6 / 8;
     $lat0 = $lat0  + (($code7-$x+1) * 2 / 3) / 8 / 10;
     $long0 = $long0 + ($code8+$y) / 8 / 10;
     $lat0 = $lat0  + (floor(($code9-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2;
     $long0 = $long0 + (($code9-1)%2-$y) / 8 / 10 / 2;
     $lat0 = (1-2*$x)*$lat0;
     $long0 = (1-2*$y)*$long0;
     $dlat = 2/3/8/10/2;
     $dlong = 1/8/10/2;
     $lat1 = sprintf("%.8f", $lat0-$dlat);
     $long1 = sprintf("%.8f", $long0+$dlong);
     $lat0 = sprintf("%.8f", $lat0);
     $long0 = sprintf("%.8f", $long0);
     $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
     return($xx);
   }
   # code 10
   #     N
   #   3 | 4
   # W - + - E
   #   1 | 2
   #     S
   if($ncode==12) { # 5th grid (12 digits)
     if(!$extension){
       # code 10
       #     N
       #   3 | 4
       # W - + - E
       #   1 | 2
       #     S
       $lat0 = $code12 * 2 / 3;
       $long0 = $code34 + 100*$z;
       $lat0 = $lat0  + ($code5 * 2 / 3) / 8;
       $long0 = $long0 + $code6 / 8;
       $lat0 = $lat0  + (($code7-$x+1) * 2 / 3) / 8 / 10;
       $long0 = $long0 +  ($code8+$y) / 8 / 10;
       $lat0 = $lat0  + (floor(($code9-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2;
       $long0 = $long0 + (($code9-1)%2-$y) / 8 / 10 / 2;
       $lat0 = $lat0  + (floor(($code10-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2 / 2;
       $long0 = $long0 + (($code10-1)%2-$y) / 8 / 10 / 2 / 2;
       $lat0 = (1-2*$x)*$lat0;
       $long0 = (1-2*$y)*$long0;
       $dlat = 2/3/8/10/2/2;
       $dlong = 1/8/10/2/2;
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
        $lat0 = $code12 * 2 / 3;
        $long0 = $code34 + 100*$z;
        $lat0 = $lat0 + ($code5 * 2 / 3) / 8;
        $long0 = $long0 + $code6 / 8;
        $lat0 = $lat0 + ($code7 * 2 / 3) / 8 / 10;
        $long0 = $long0 + $code8 / 8 / 10;
        $lat0 = $lat0 + (($codeex9-$x+1) * 2 / 3) / 8 / 10 / 10;
        $long0 = $long0 + ($codeex10+$y) / 8 / 10 / 10;
        $lat0 = (1-2*$x)*$lat0;
        $long0 = (1-2*$y)*$long0;
        $dlat = 2/3/8/10/10;
        $dlong = 1/8/10/10;
     }
     $lat1 = sprintf("%.10f", $lat0-$dlat);
     $long1 = sprintf("%.10f", $long0+$dlong);
     $lat0 = sprintf("%.10f", $lat0);
     $long0 = sprintf("%.10f", $long0);
     $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
     return($xx);
   }
   if($ncode==13) { # 6rd grid (13 digits)
      if(!$extension){
        # code 11
        #     N
        #   3 | 4
        # W - + - E
        #   1 | 2
        #     S
        $lat0 = $code12 * 2 / 3;
        $long0 = $code34 + 100*$z;
        $lat0 = $lat0  + ($code5 * 2 / 3) / 8;
        $long0 = $long0 + $code6 / 8;
        $lat0 = $lat0  + (($code7-$x+1) * 2 / 3) / 8 / 10;
        $long0 = $long0 + ($code8+$y) / 8 / 10;
        $lat0 = $lat0  + (floor(($code9-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2;
        $long0 = $long0 + (($code9-1)%2-$y) / 8 / 10 / 2;
        $lat0 = $lat0  + (floor(($code10-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2 / 2;
        $long0 = $long0 + (($code10-1)%2-$y) / 8 / 10 / 2 / 2;
        $lat0 = $lat0  + (floor(($code11-1)/2)+$x-1) * 2 / 3 / 8 / 10 / 2 / 2 / 2;
        $long0 = $long0 + (($code11-1)%2-$y) / 8 / 10 / 2 / 2 / 2;
        $lat0 = (1-2*$x)*$lat0;
        $long0 = (1-2*$y)*$long0;
        $dlat = 2/3/8/10/2/2/2;
        $dlong = 1/8/10/2/2/2;
      }else{ # Extended 100m grid square
      #    Code10
#      4*****
#      3*****
#      2*****
#      1*****
#      0*****
#       01234 Code11
        $lat0  = $code12 * 2 / 3;  
        $long0 = $code34 + 100*$z;
        $lat0 = $lat0  + ($code5 * 2 / 3) / 8; 
        $long0 = $long0 +  $code6 / 8;
        $lat0 = $lat0  + (($code7-$x+1) * 2 / 3) / 8 / 10;
        $long0 = $long0 +  ($code8+$y) / 8 / 10;
        $lat0 = $lat0  + (floor(($code9-1)/2)+2*$x-2) * 2 / 3 / 8 / 10 / 2;
        $long0 = $long0 + (($code9-1)%2-2*$y) / 8 / 10 / 2;
        $lat0 = $lat0 + ($codeex10-$x+1) * 2 / 3 / 8 / 10 / 2 / 5;
        $long0 = $long0 + ($codeex11+$y) / 8 / 10 / 2 / 5;
        $lat0 = (1-2*$x)*$lat0;
        $long0 = (1-2*$y)*$long0;
        $dlat = 2/3/8/10/2/5;
        $dlong = 1/8/10/2/5;
      }
      $lat1 = sprintf("%.8f", $lat0-$dlat);
      $long1 = sprintf("%.8f", $long0+$dlong);
      $lat0 = sprintf("%.8f", $lat0);
      $long0 = sprintf("%.8f", $long0);
      $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
      return($xx);
   }
   if($ncode==14) {  # (14 digits)
      if(!$extension){
      }else{ # Extended 10m grid square code
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
        $lat0 = $code12 * 2 / 3;
        $long0 = $code34 + 100*$z;
        $lat0 = $lat0  + ($code5 * 2 / 3) / 8;
        $long0 = $long0 + $code6 / 8;
        $lat0 = $lat0 + $code7 * 2 / 3 / 8 / 10;
        $long0 = $long0 + $code8 / 8 / 10;
        $lat0 = $lat0 + ($codeex9 * 2 / 3) / 8 / 10 / 10;
        $long0 = $long0 + $codeex10 / 8 / 10 / 10;
        $lat0 = $lat0 + (($codeex11-$x+1) * 2 / 3) / 8 / 10 / 10 / 10;
        $long0 = $long0 + ($codeex12+$y) / 8 / 10 / 10 / 10;
        $lat0 = (1-2*$x)*$lat0;
        $long0 = (1-2*$y)*$long0;
        $dlat = 2/3/8/10/10/10;
        $dlong = 1/8/10/10/10;
      }
      $lat1 = sprintf("%.12f", $lat0-$dlat); 
      $long1 = sprintf("%.12f", $long0+$dlong);
      $lat0 = sprintf("%.12f", $lat0); 
      $long0 = sprintf("%.12f", $long0);
      $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
      return($xx);
    }
    if($ncode==16) {  # (16 digits)
      if(!$extension){
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
        $lat0 = $code12 * 2 / 3;
        $long0 = $code34 + 100*$z;
        $lat0 = $lat0 + ($code5 * 2 / 3) / 8;
        $long0 = $long0 + $code6 / 8;
        $lat0 = $lat0  + ($code7 * 2 / 3) / 8 / 10;
        $long0 = $long0 + $code8 / 8 / 10;
        $lat0 = $lat0  + ($codeex9 * 2 / 3) / 8 / 10 / 10;
        $long0 = $long0 + $codeex10 / 8 / 10 / 10;
        $lat0 = $lat0  + ($codeex11 * 2 / 3) / 8 / 10 / 10 / 10;
        $long0 = $long0 + $codeex12 / 8 / 10 / 10 / 10;
        $lat0 = $lat0  + (($codeex13-$x+1) * 2 / 3) / 8 / 10 / 10 / 10 / 10;
        $long0 = $long0 + ($codeex14+$y) / 8 / 10 / 10 / 10 / 10;
        $lat0 = (1-2*$x)*$lat0;
        $long0 = (1-2*$y)*$long0;
        $dlat = 2/3/8/10/10/10/10;
        $dlong = 1/8/10/10/10/10;
      }
      $lat1 = sprintf("%.14f", $lat0-$dlat);
      $long1 = sprintf("%.14f", $long0+$dlong);
      $lat0 = sprintf("%.14f", $lat0);
      $long0 = sprintf("%.14f", $long0);
      $xx = array("lat0"=>floatval($lat0), "long0"=>floatval($long0), "lat1"=>floatval($lat1), "long1"=>floatval($long1));
      return($xx);      
   }
   $xx = array("lat0"=>floatval(99999), "long0"=>floatval(99999), "lat1"=>floatval(99999), "long1"=>floatval(99999));
   return($xx);
}

function meshcode_to_latlong($meshcode,$extension=false){
    $res = meshcode_to_latlong_grid($meshcode,$extension);
    $xx = array();
    $xx["lat"] = $res["lat0"];
    $xx["long"] = $res["long0"];
    return($xx);
}

function meshcode_to_latlong_NW($meshcode,$extension=false){
    $res = meshcode_to_latlong_grid($meshcode,$extension);
    $xx = array();
    $xx["lat"] = $res["lat0"];
    $xx["long"] = $res["long0"];
    return($xx);
}

function meshcode_to_latlong_SW($meshcode,$extension=false){
    $res = meshcode_to_latlong_grid($meshcode,$extension);
    $xx = array();
    $xx["lat"] = $res["lat1"];
    $xx["long"] = $res["long0"];
    return($xx);
}

function meshcode_to_latlong_NE($meshcode,$extension=false){
    $res = meshcode_to_latlong_grid($meshcode,$extension);
    $xx = array();
    $xx["lat"] = $res["lat0"];
    $xx["long"] = $res["long1"];
    return($xx);
}

function meshcode_to_latlong_SE($meshcode,$extension=false){
    $res = meshcode_to_latlong_grid($meshcode,$extension);
    $xx = array();
    $xx["lat"] = $res["lat1"];
    $xx["long"] = $res["long1"];
    return($xx);
}

function cal_meshcode6($latitude,$longitude){
   if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
    return("999999999999");
   }
   if($latitude < 0){
     $o = 4;
   } else {
     $o = 0;
   }
   if($longitude < 0){
     $o = $o + 2;
   }
   if(abs($longitude) >= 100){
     $o = $o + 1;
   }
   $z = $o % 2;
   $y = (($o - $z)/2) % 2;
   $x = ($o - 2*$y - $z)/4;
   #
   $o = $o + 1;
   #
   $latitude = (1-2*$x)*$latitude;
   $longitude = (1-2*$y)*$longitude;
   #
   $p = floor($latitude*60/40);
   $a = $latitude*60-$p*40;
   $q = floor($a/5);
   $b = $a-$q*5;
   $r = floor($b*60/30);
   $c = $b*60-$r*30;
   $s2u = floor($c/15);
   $d = $c-$s2u*15;
   $s4u = floor($d/7.5);
   $e = $d-$s4u*7.5;
   $s8u = floor($e/3.75);
   $u = floor($longitude-100*$z);
   $f = $longitude-100*$z-$u;
   $v = floor($f*60/7.5);
   $g = $f*60-$v*7.5;
   $w = floor($g*60/45);
   $h = $g*60-$w*45;
   $s2l = floor($h/22.5);
   $i = $h-$s2l*22.5;
   $s4l = floor($i/11.25);
   $j = $i-$s4l*11.25;
   $s8l = floor($j/5.625);
   $s2 = $s2u*2+$s2l+1;
   $s4 = $s4u*2+$s4l+1;
   $s8 = $s8u*2+$s8l+1;
   #
   if($u < 10){
     if($p < 10){
       $mesh = $o."00".$p."0".$u.$q.$v.$r.$w.$s2.$s4.$s8;
     }else{
       if($p < 100){
          $mesh = $o."0".$p."0".$u.$q.$v.$r.$w.$s2.$s4.$s8;
       }
       else{
         $mesh = $o.$p."0".$u.$q.$v.$r.$w.$s2.$s4.$s8;
       }
     }
   }
   else{
     if($p < 10){
       $mesh = $o."00".$p.$u.$q.$v.$r.$w.$s2.$s4.$s8;
     }else{
       if($p < 100){
         $mesh = $o."0".$p.$u.$q.$v.$r.$w.$s2.$s4.$s8;
       }
       else{
         $mesh = $o.$p.$u.$q.$v.$r.$w.$s2.$s4.$s8;
       }
     }
   }
   return($mesh);
}
function cal_meshcode($latitude,$longitude){
    return(cal_meshcode3($latitude,$longitude));
}
function cal_meshcode1($latitude,$longitude){
    if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
      return("999999");
    }
    return(substr(cal_meshcode6($latitude,$longitude),0,6));
}
function cal_meshcode2($latitude,$longitude){
    if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
     return("99999999");
    }
    return(substr(cal_meshcode6($latitude,$longitude),0,8));
}
function cal_meshcode3($latitude,$longitude){
    if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
      return("9999999999");
    }
    return(substr(cal_meshcode6($latitude,$longitude),0,10));
}
function cal_meshcode4($latitude,$longitude){
    if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
      return("99999999999");
    }
    return(substr(cal_meshcode6($latitude,$longitude),0,11));
}
function cal_meshcode5($latitude,$longitude){
    if ($latitude < -90 or $latitude >90 or $longitude < -180 or $longitude >180){
      return("999999999999");
    }
    return(substr(cal_meshcode6($latitude,$longitude),0,12));
}

# calculate an extended 100m (1km / 10) grid square code (12 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude
function cal_meshcode_ex100m_12($latitude, $longitude){
  if($latitude < 0){
    $o = 4;
  }
  else{
    $o = 0;
  }
  if($longitude < 0){
    $o = $o + 2;
  }
  if(abs($longitude) >= 100) $o = $o + 1;
  $z = $o % 2;
  $y = (($o - $z)/2) % 2;
  $x = ($o - 2*$y - $z)/4;
#
  $o = $o + 1;
#
  $latitude = (1-2*$x)*$latitude;
  $longitude = (1-2*$y)*$longitude;
  #
  $p = floor($latitude*60/40);
  $a = $latitude*60-$p*40;
  $q = floor($a/5);
  $b = $a-$q*5;
  $r = floor($b*60/30);
  $c = $b*60-$r*30;
  $s = floor($c/3);
  $d = $c-$s*3;
  $u = floor($longitude-100*$z);
  $f = $longitude-100*$z-$u;
  $v = floor($f*60/7.5);
  $g = $f*60-$v*7.5;
  $w = floor($g*60/45);
  $h = $g*60-$w*45;
  $xx = floor($h/4.5);
  $i = $h-$xx*4.5;
  #
  if($u < 10){
    if($p < 10){
      $mesh = $o."00".$p."0".$u.$q.$v.$r.$w.$s.$xx;
    }else{
      if($p < 100){
        $mesh = $o."0".$p."0".$u.$q.$v.$r.$w.$s.$xx;
      }
      else{
        $mesh = $o.$p."0".$u.$q.$v.$r.$w.$s.$xx;
      }
    }
  }
  else{
    if($p < 10){
      $mesh = $o."00".$p.$u.$q.$v.$r.$w.$s.$xx;
    }else{
      if($p < 100){
        $mesh = $o."0".$p.$u.$q.$v.$r.$w.$s.$xx;
      }
      else{
        $mesh = $o.$p.$u.$q.$v.$r.$w.$s.$xx;
      }
    }
  }
  return($mesh);
}
# calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude
function cal_meshcode_ex100m_13($latitude, $longitude){
  if($latitude < 0){
    $o = 4;
  }
  else{
    $o = 0;
  }
  if($longitude < 0){
    $o = $o + 2;
  }
  if(abs($longitude) >= 100) $o = $o + 1;
  $z = $o % 2;
  $y = (($o - $z)/2) % 2;
  $x = ($o - 2*$y - $z)/4;
#
  $o = $o + 1;
#
  $latitude = (1-2*$x)*$latitude;
  $longitude = (1-2*$y)*$longitude;
  #
  $p = floor($latitude*60/40);
  $a = $latitude*60-$p*40;
  $q = floor($a/5);
  $b = $a-$q*5;
  $r = floor($b*60/30);
  $c = $b*60-$r*30;
  $s2u = floor($c/15);
  $d = $c-$s2u*15;
  $et = floor($d/3);
  $e = $d-$et*3;
  $u = floor($longitude-100*$z);
  $f = $longitude-100*$z-$u;
  $v = floor($f*60/7.5);
  $g = $f*60-$v*7.5;
  $w = floor($g*60/45);
  $h = $g*60-$w*45;
  $s2l = floor($h/22.5);
  $i = $h-$s2l*22.5;
  $jt = floor($i/4.5);
  $j = $i*0.5-$jt*4.5;
  $s2 = $s2u*2+$s2l+1;
  #
  if($u < 10){
    if($p < 10){
      $mesh = $o."00".$p."0".$u.$q.$v.$r.$w.$s2.$et.$jt;
    }else{
      if($p < 100){
        $mesh = $o."0".$p."0".$u.$q.$v.$r.$w.$s2.$et.$jt;
      }
      else{
        $mesh = $o.$p."0".$u.$q.$v.$r.$w.$s2.$et.$jt;
      }
    }
  }
  else{
    if($p < 10){
      $mesh = $o."00".$p.$u.$q.$v.$r.$w.$s2.$et.$jt;
    }else{
      if($p < 100){
        $mesh = $o."0".$p.$u.$q.$v.$r.$w.$s2.$et.$jt;
      }
      else{
        $mesh = $o.$p.$u.$q.$v.$r.$w.$s2.$et.$jt;
      }
    }
  }
  return($mesh);
}
# calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) from a geographical position (latitude, longitude) - 0.3 arc-second for latitude and 0.45 arc-second for longitude
function cal_meshcode_ex10m_14($latitude, $longitude){
  if($latitude < 0){
    $o = 4;
  }
  else{
    $o = 0;
  }
  if($longitude < 0){
    $o = $o + 2;
  }
  if(abs($longitude) >= 100) $o = $o + 1;
  $z = $o % 2;
  $y = (($o - $z)/2) % 2;
  $x = ($o - 2*$y - $z)/4;
#
  $o = $o + 1;
#
  $latitude = (1-2*$x)*latitude;
  $longitude = (1-2*$y)*longitude;
  #
  $p = floor($latitude*60/40);
  $a = $latitude*60-$p*40;
  $q = floor($a/5);
  $b = $a-$q*5;
  $r = floor($b*60/30);
  $c = $b*60-$r*30;
  $s = floor($c/3);
  $d = $c-$s*3;
  $t = floor($d/0.3);
  $e = $d-$t*0.3;
  $u = floor($longitude-100*$z);
  $f = $longitude-100*$z-$u;
  $v = floor($f*60/7.5);
  $g = $f*60-$v*7.5;
  $w = floor($g*60/45);
  $h = $g*60-$w*45;
  $xx = floor($h/4.5);
  $i = $h-$xx*4.5;
  $yy = floor($i/0.45);
  $j = $i-$yy*0.45;
  #
  if($u < 10){
    if($p < 10){
      $mesh = $o."00".$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
    }else{
      if($p < 100){
        $mesh = $o."0".$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
      }
      else{
        $mesh = $o.$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
      }
    }
  }
  else{
    if($p < 10){
      $mesh = $o."00".$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
    }else{
      if($p < 100){
        $mesh = $o."0".$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
      }
      else{
        $mesh = $o.$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy;
      }
    }
  }
  return($mesh);
}
# calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) from a geographical position (latitude, longitude) - 0.03 arc-second for latitude and 0.045 arc-second for longitude
function cal_meshcode_ex1m_16($latitude, $longitude){
  if($latitude < 0){
    $o = 4;
  }
  else{
    $o = 0;
  }
  if($longitude < 0){
    $o <- $o + 2;
  }
  if(abs($longitude) >= 100) $o = $o + 1;
  $z = $o % 2;
  $y = (($o - $z)/2) % 2;
  $x = ($o - 2*$y - $z)/4;
#
  $o = $o + 1;
#
  $latitude = (1-2*$x)*$latitude;
  $longitude = (1-2*$y)*$longitude;
  #
  $p = floor($latitude*60/40);
  $a = $latitude*60-$p*40;
  $q = floor($a/5);
  $b = $a-$q*5;
  $r = floor($b*60/30);
  $c = $b*60-$r*30;
  $s = floor($c/3);
  $d = $c-$s*3;
  $t = floor($d/0.3);
  $e = $d-$t*0.3;
  $tt = floor($e/0.03);
  $ee = $e-$tt*0.03;
  $u = floor($longitude-100*$z);
  $f = $longitude-100*$z-$u;
  $v = floor($f*60/7.5);
  $g = $f*60-$v*7.5;
  $w = floor($g*60/45);
  $h = $g*60-$w*45;
  $xx = floor($h/4.5);
  $i = $h-$xx*4.5;
  $yy = floor($i/0.45);
  $j = $i-$yy*0.45;
  $zz = floor($j/0.045);
  $k = $j-$zz*0.045;
  #
  if($u < 10){
    if($p < 10){
      $mesh = $o."00".$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
    }else{
      if($p < 100){
        $mesh = $o."0".$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
      }
      else{
        $mesh = $o.$p."0".$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
      }
    }
  }
  else{
    if($p < 10){
      $mesh = $o."00".$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
    }else{
      if($p < 100){
        $mesh = $o."0".$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
      }
      else{
        $mesh = $o.$p.$u.$q.$v.$r.$w.$s.$xx.$t.$yy.$tt.$zz;
      }
    }
  }
  return($mesh);
}
#
function Vincenty ($latitude1,$longitude1,$latitude2,$longitude2){
   if(abs($latitude1-$latitude2)+abs($longitude1-$longitude2) < 1e-5) return(0.0);
   # WGS84
   $f = 1/298.257223563;
   $a = 6378137.0;
   $b = 6356752.314245;
   #
   $L = ($longitude1 - $longitude2)/180*pi();
   $U1 = atan((1-$f)*tan($latitude1/180*pi()));
   $U2 = atan((1-$f)*tan($latitude2/180*pi()));
   $lambda = $L;
   $dlambda = 10;
   while(abs($dlambda) > 1e-12){
      $cs = cos($U2)*sin($lambda);
      $cscc = cos($U1)*sin($U2)-sin($U1)*cos($U2)*cos($lambda);
      $sinsigma = sqrt($cs*$cs + $cscc*$cscc);
      $cossigma = sin($U1)*sin($U2)+cos($U1)*cos($U2)*cos($lambda);
      $sigma = atan($sinsigma/$cossigma);
      $sinalpha = cos($U1)*cos($U2)*sin($lambda)/$sinsigma;
      $cos2alpha = 1 - $sinalpha*$sinalpha;
      if($cos2alpha == 0.0){
         $C = 0.0;
         $lambda0 = $L + $f*$sinalpha*$sigma;
      }else{
         $cos2sigmam = $cossigma - 2*sin($U1)*sin($U2)/$cos2alpha;
         $C = $f/16*$cos2alpha*(4+$f*(4-3*$cos2alpha));
         $lambda0 = $L + (1-$C)*$f*$sinalpha*($sigma + $C*$sinsigma*($cos2sigmam + $C*$cossigma*(-1+2*$cos2sigmam*$cos2sigmam)));
      }
      $dlambda = $lambda0 - $lambda;
      $lambda = $lambda0;
   }
   if($C == 0.0){
      $A = 1.0;
      $dsigma = 0.0;
   }else{
      $u2 = $cos2alpha * ($a*$a-$b*$b)/($b*$b);
      $A = 1 + $u2/16384*(4096 + $u2 * (-768 + $u2*(320-175*$u2)));
      $B = $u2/1024*(256+$u2*(-128+$u2*(74-47*$u2)));
      $dsigma = $B*$sinsigma*($cos2sigmam + 1/4*$B*($cossigma*(-1+2*$cos2sigmam*$cos2sigmam)-1/6*$B*$cos2sigmam*(-3+4*$sinsigma*$sinsigma)*(-3+4*$cos2sigmam*$cos2sigmam)));
   }
   $s = $b*$A*($sigma-$dsigma);
   return($s);
}

function cal_area_from_meshcode($meshcode,$extension=false){
   $latlong = meshcode_to_latlong_grid($meshcode,$extension);
   return(cal_area_from_latlong($latlong));
}

function cal_area_from_latlong($latlong){
   $W1 = Vincenty($latlong['lat0'],$latlong['long0'],$latlong['lat0'],$latlong['long1']);
   $W2 = Vincenty($latlong['lat1'],$latlong['long0'],$latlong['lat1'],$latlong['long1']);
   $H = Vincenty($latlong['lat0'],$latlong['long0'],$latlong['lat1'],$latlong['long0']);
   $A=($W1+$W2)*$H*0.5;
   $xx=array();
   $xx['W1'] = $W1;
   $xx['W2'] = $W2;
   $xx['H'] = $H;
   $xx['A'] = $A;
   return($xx);
}
?>

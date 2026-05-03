//
// Copyright (c) 2015-2022 Research Institute for World Grid Squares 
// Prof. Dr. Aki-Hiro Sato (Yokohama City University)
// All rights reserved. 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// javascript fnctions to calculate the world grid square code.
// The world grid square code computed by this library 
// is compatible to Japanese grid square code (JIS X 0410). 
//
// This library for World grid square code complies ISO 24108-1:
// https://www.iso.org/standard/87706.html
//
// Version 1.0: Released on 2 February 2017
// Version 1.01: Released on 6 February 2017
// Version 1.1: Released on 21 November 2019
// Version 1.6: Released on 09 October 2020
// Version 1.7: Released on 05 March 2022
// Version 1.75: Released on 13 April 2022
//
// Written by Prof. Dr. Aki-Hiro Sato    
// Graduate School of Data Science, Department of Data Science,
// Yokohama City University
//
// Contact:
// Address: 22-2, Seto, Kanazawa-ku, Yokohama, Kanagawa 236-0027 Japan
// E-mail: ahsato@yokohama-cu.ac.jp
// TEL: +81-45-787-2208
//
// Three types of functions are defined in this library.
// 1. calculate representative geographical position(s) (latitude, longitude) of a grid square from a grid square code
// 2. calculate a grid square code from a geographical position (latitude, longitude)
// 3. calculate geodesic distance and size of grid square (representative lengths and area)
//
// 1.
//
// meshcode_to_latlong(meshcode, extension=false)
// : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
// meshcode_to_latlong_NW(meshcode, extension=false)
// : calculate northen western geographic position of the grid (latitude, longitude) from meshcode
// meshcode_to_latlong_SW(meshcode, extension=false)
// : calculate sourthern western geographic position of the grid (latitude, longitude) from meshcode
// meshcode_to_latlong_NE(meshcode, extension=false)
// : calculate northern eastern geographic position of the grid (latitude, longitude) from meshcode
// meshcode_to_latlong_SE(meshcode, extension=false)
// : calculate sourthern eastern geographic position of the grid (latitude, longitude) from meshcode
// meshcode_to_latlong_grid(meshcode, extension=false)
// : calculate northern western and sourthern eastern geographic positions of the grid (latitude0, longitude0, latitude1, longitude1) from meshcode
//
// 2.
//
// : calculate a basic (1km) grid square code (10 digits) from a geographical position (latitude, longitude)
// cal_meshcode1(latitude,longitude)
// : calculate an 80km grid square code (6 digits) from a geographical position (latitude, longitude)
// cal_meshcode2(latitude,longitude)
// : calculate a 10km grid square code (8 digits) from a geographical position (latitude, longitude)
// cal_meshcode3(latitude,longitude)
// : calculate a 1km grid square code (10 digits) from a geographical position (latitude, longitude)
// cal_meshcode4(latitude,longitude)
// : calculate a 500m grid square code (11 digits) from a geographical position (latitude, longitude)
// cal_meshcode5(latitude,longitude)
// : calculate a 250m grid square code (12 digits) from a geographical position (latitude, longitude)
// cal_meshcode6(latitude,longitude)
// : calculate a 125m grid square code (13 digits) from a geographical position (latitude, longitude)
//
// This grid square code set is not included in JIS X0410 directly but useful.
//
// cal_meshcode_ex100m_12(latitude,longitude)
// : calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude 
//
// cal_meshcode_ex100m_13(latitude,longitude)
// : calculate an extended 100m (500m / 5) grid square code (13 digits) from a geographical position (latitude, longitude) - 3 arc-second for latitude and 4.5 arc-second for longitude 
//
// cal_meshcode_ex10m_14(latitude,longitude)
// : calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) from a geographical position (latitude, longitude) - 0.3 arc-second for latitude and 0.45 arc-second for longitude
//
// cal_meshcode_ex1m_16(latitude,longitude)
// : calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) from a geographical position (latitude, longitude) - 0.03 arc-second for latitude and 0.045 arc-second for longitude
//
// Structure of the world grid square code with compatibility to JIS X0410
// A : area code (1 digit) A takes 1 to 8
// ABBBBB : 80km grid square code (40 arc-minutes for latitude, 1 arc-degree for longitude) (6 digits)
// ABBBBBCC : 10km grid square code (5 arc-minutes for latitude, 7.5 arc-minutes for longitude) (8 digits)
// ABBBBBCCDD : 1km grid square code (30 arc-seconds for latitude, 45 arc-secondes for longitude) (10 digits)
// ABBBBBCCDDE : 500m grid square code (15 arc-seconds for latitude, 22.5 arc-seconds for longitude) (11 digits)
// ABBBBBCCDDEF : 250m grid square code (7.5 arc-seconds for latitude, 11.25 arc-seconds for longitude) (12 digits)
// ABBBBBCCDDEFG : 125m grid square code (3.75 arc-seconds for latitude, 5.625 arc-seconds for longitude) (13 digits)
// ABBBBBCCDDHH : Extended 100m grid square code with 12 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (12 digits)
// ABBBBBCCDDEHH : Extended 100m grid square code with 13 digits (3 arc-seconds for latitude, 4.5 arc-seconds for longitude) (13 digits)
// ABBBBBCCDDHHGG : Extended 10m grid square code with 14 digits (0.3 arc-seconds for latitude, 0.45 arc-seconds for longitude) (14 digits)
// ABBBBBCCDDHHGGII : Extended 1m grid square code with 16 digits (0.03 arc-seconds for latitude, 0.045 arc-seconds for longitude) (16 digits)
//
// 3.
//
// Calculate geodesic distance and size of world grid square 
//
// T. Vincenty, ``Direct and Inverse Solutions of Geodesics
// on the Ellipsoid with application of nested equations'',
// Survey Review XXIII, Vol 176 (1975) Vol. 88-93.
//
// Vincenty(latitude1, longitude1, latitude2, longitude)
// : calculate geodesitc distance between two points (latitude1, longitude1) and (latitude2, longitude2) placed on the WGS84 Earth ellipsoid based on the Vincenty's formulae (1975)
// cal_area_from_meshcode(meshcode, extension=false)
// : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoide A) of world grid square indicated by meshcode
// cal_area_from_latlong(latlong)
// : calculate size (northern west-to-east span H1, sothern west-to-east span H2, north-to-south span W, and area approximated by trapezoid A) of a trapezoid on the WGS84 Earth ellipoid indicated by (latlong$lat0, latlong$long0, latlong$lat1, latlong$long1)
//

function meshcode_to_latlong_grid(meshcode, extension=false){
    // The extension option determines high resolution gird square code for 13 digits (100m x 100m grid squares 3 arc-seconds for latitude, 4.5 arc-seconds for longitude)
    //
    code = meshcode + '';
    ncode = code.length;
//    console.log("ncode ============"+ncode);
    if (ncode < 6){
	    return(null);
    }else{// more than 1st grid
	    code0 = parseInt(code.substr(0, 1)); // code0 : 1 to 8
	    code0 = code0 - 1; // transforming code0 from 0 to 7
	    code12 = code.substr(1, 3);
	    if(code12.substr(0,2)=='00'){
                code12 = parseInt(code.substr(3, 1));
	    }else{
            if(code12.substr(0,1)=='0'){
		code12 = parseInt(code.substr(2, 2));
            }
            else{
		code12 = parseInt(code.substr(1, 3));
            }
	    }
	    if(code.substr(4,1)=='0'){
            code34 = parseInt(code.substr(5, 1));
	    }
	    else{
            code34 = parseInt(code.substr(4, 2));
	    }
    }
    
    if (ncode >= 8) { // more than 2nd grid
	    code5 = parseInt(code.substr(6, 1));
	    code6 = parseInt(code.substr(7, 1));
	}

    if (ncode >= 10) { // more than 3rd grid
	    code7 = parseInt(code.substr(8, 1));
	    code8 = parseInt(code.substr(9, 1));
    }

    if(!extension){ // standard grid codes
        if (ncode >= 11) { // more than 4th grid
            code9 = parseInt(code.substr(10, 1));
        }
        if (ncode >= 12) { // more than 5th grid
            code10 = parseInt(code.substr(11, 1));
        }
        if (ncode >= 13) { // more than 6th grid
            code11 = parseInt(code.substr(12, 1)); 
        }
    } else {
        if(ncode >= 12){ // extended 100m grid square code (12digits) 
            codeex9 = parseInt(code.substr(10, 1));
            codeex10 = parseInt(code.substr(11, 1));
        }
        if(ncode >= 13){ // extended 100m grid square code (13digits) 
            code9 = parseInt(code.substr(10, 1));
            codeex10 = parseInt(code.substr(11, 1));
            codeex11 = parseInt(code.substr(12, 1));
        }
        if(ncode >= 14){ // extended 10m grid square code (14digits) 
            codeex11 = parseInt(code.substr(12, 1));
            codeex12 = parseInt(code.substr(13, 1));
        }
        if(ncode >= 16){ // extended 1m grid square code (16digits) 
            codeex13 = parseInt(code.substr(14, 1));
            codeex14 = parseInt(code.substr(15, 1));
        }
	}
    
    // 0'th grid
    z = code0 % 2;
    y = ((code0 - z)/2) % 2;
    x = (code0 - 2*y - z)/4;

    switch(ncode){
    case 6: // 1st grid (6 digits)
	lat0 = (code12-x+1) * 2 / 3;        
	long0 = (code34+y) + 100*z;
	lat0 = (1-2*x)*lat0;        
	long0 = (1-2*y)*long0;
	dlat = 2/3;
	dlong = 1;
    lat1 = myformat8(lat0-dlat);
    long1 = myformat8(long0+dlong);
    lat0 = myformat8(lat0);  
    long0 = myformat8(long0);
	break;
    case 8: // 2nd grid (8 digits)
	lat0 = code12 * 2 / 3;
	long0 = code34 + 100*z;
	lat0 = lat0  + ((code5-x+1) * 2 / 3) / 8; 
	long0 = long0 +  (code6+y) / 8;
	lat0 = (1-2*x) * lat0;
	long0 = (1-2*y) * long0;
	dlat = 2/3/8;
	dlong = 1/8;
    lat1 = myformat8(lat0-dlat);
    long1 = myformat8(long0+dlong);
    lat0 = myformat8(lat0);  
    long0 = myformat8(long0);
	break;
    case 10: // 3rd grid (10 digits)
	lat0 = code12 * 2 / 3;  
	long0 = code34 + 100*z;
	lat0 = lat0 + (code5 * 2 / 3) / 8; 
	long0 = long0 +  code6 / 8;
	lat0 = lat0 + ((code7-x+1) * 2 / 3) / 8 / 10;
	long0 = long0 + (code8+y) / 8 / 10;
	lat0 = (1-2*x)*lat0;
	long0 = (1-2*y)*long0;
	dlat = 2/3/8/10;
	dlong = 1/8/10;
    lat1 = myformat8(lat0-dlat);
    long1 = myformat8(long0+dlong);
    lat0 = myformat8(lat0);  
    long0 = myformat8(long0);
    break;
    case 11: // 4th grid (11 digits)
	// code 9
	//     N
	//   3 | 4
	// W - + - E
	//   1 | 2
	//     S
	lat0 = code12 * 2 / 3;  
	long0 = code34 + 100*z;
	lat0 = lat0 + (code5 * 2 / 3) / 8; 
	long0 = long0 + code6 / 8;
	lat0 = lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
	long0 = long0 + (code8+y) / 8 / 10;
	lat0 = lat0  + (Math.floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2;
	long0 = long0 + ((code9-1)%2-y) / 8 / 10 / 2;
	lat0 = (1-2*x)*lat0;
	long0 = (1-2*y)*long0;
	dlat = 2/3/8/10/2;
	dlong = 1/8/10/2;
    lat1 = myformat8(lat0-dlat);
    long1 = myformat8(long0+dlong);
    lat0 = myformat8(lat0);  
    long0 = myformat8(long0);
	break;
    case 12 : // 5th grid (12 digits)
        if(!extension){
    	    // code 10
    	    //     N
    	    //   3 | 4
    	    // W - + - E
    	    //   1 | 2
    	    //     S
    	    lat0 = code12 * 2 / 3;  
    	    long0 = code34 + 100*z;
    	    lat0 = lat0  + (code5 * 2 / 3) / 8; 
    	    long0 = long0 + code6 / 8;
    	    lat0 = lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
    	    long0 = long0 +  (code8+y) / 8 / 10;
    	    lat0 = lat0  + (Math.floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2;
    	    long0 = long0 + ((code9-1)%2-y) / 8 / 10 / 2;
    	    lat0 = lat0  + (Math.floor((code10-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2;
    	    long0 = long0 + ((code10-1)%2-y) / 8 / 10 / 2 / 2;
    	    lat0 = (1-2*x)*lat0;
    	    long0 = (1-2*y)*long0;
    	    dlat = 2/3/8/10/2/2;
	        dlong = 1/8/10/2/2;
        }else{ // extended 100m grid square code (12 digits)
            //   Code9
            //      9**********
            //      8**********
            //      7**********
            //      6**********
            //      5**********
            //      4**********
            //      3**********
            //      2**********
            //      1**********
            //      0**********
            //      0123456789 Code10
            lat0  = code12 * 2 / 3;  
            long0 = code34 + 100*z;
            lat0  = lat0  + (code5 * 2 / 3) / 8;
            long0 = long0 +  code6 / 8;
            lat0  = lat0  + (code7 * 2 / 3) / 8 / 10;
            long0 = long0 + code8 / 8 / 10;
            lat0  = lat0  + ((codeex9-x+1) * 2 / 3) / 8 / 10 / 10;
            long0 = long0 + (codeex10+y) / 8 / 10 / 10;
            lat0 = (1-2*x)*lat0;
            long0 = (1-2*y)*long0;
            dlat = 2/3/8/10/10;
            dlong = 1/8/10/10;
        }
        lat1 = myformat10(lat0-dlat);
        long1 = myformat10(long0+dlong);
        lat0 = myformat10(lat0);  
        long0 = myformat10(long0);    
	break;
    case 13: 
    	if(!extension){ // 6rd grid (13 digits)
	        // code 11
	        //     N
	        //   3 | 4
	        // W - + - E
	        //   1 | 2
	        //     S
	        lat0 = code12 * 2 / 3;  
	        long0 = code34 + 100*z;
	        lat0 = lat0  + (code5 * 2 / 3) / 8; 
	        long0 = long0 + code6 / 8;
	        lat0 = lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
	        long0 = long0 + (code8+y) / 8 / 10;
	        lat0 = lat0  + (Math.floor((code9-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2;
	        long0 = long0 + ((code9-1)%2-y) / 8 / 10 / 2;
	        lat0 = lat0  + (Math.floor((code10-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2;
	        long0 = long0 + ((code10-1)%2-y) / 8 / 10 / 2 / 2;
	        lat0 = lat0  + (Math.floor((code11-1)/2)+x-1) * 2 / 3 / 8 / 10 / 2 / 2 / 2;
	        long0 = long0 + ((code11-1)%2-y) / 8 / 10 / 2 / 2 / 2;
	        lat0 = (1-2*x)*lat0;
	        long0 = (1-2*y)*long0;
	        dlat = 2/3/8/10/2/2/2;
	        dlong = 1/8/10/2/2/2;
	    }else{ // Extended 100m grid square code (13 digits)
	        //   Code10
	        //      4*****
	        //      3*****
	        //      2*****
	        //      1*****
	        //      0*****
	        //       01234 Code11
            lat0  = code12 * 2 / 3;  
            long0 = code34 + 100*z;
            lat0  = lat0  + (code5 * 2 / 3) / 8; 
            long0 = long0 +  code6 / 8;
            lat0 = lat0  + ((code7-x+1) * 2 / 3) / 8 / 10;
            long0 = long0 +  (code8+y) / 8 / 10;
            lat0  = lat0  + (Math.floor((code9-1)/2)+2*x-2) * 2 / 3 / 8 / 10 / 2;
            long0 = long0 + ((code9-1)%2-2*y) / 8 / 10 / 2;
            lat0  = lat0  + (codeex10-x+1) * 2 / 3 / 8 / 10 / 2 / 5;
            long0 = long0 + (codeex11+y) / 8 / 10 / 2 / 5;
            lat0 = (1-2*x)*lat0;
            long0 = (1-2*y)*long0;
            dlat = 2/3/8/10/2/5;
            dlong = 1/8/10/2/5;
        }
        lat1 = myformat10(lat0-dlat);
        long1 = myformat10(long0+dlong);
        lat0 = myformat10(lat0);  
        long0 = myformat10(long0);    
	break;
    case 14:
        if(!extension){ 
        }else{ // Extended 10m grid square code (14digit)
            //   Code11
            //     9**********
            //     8**********
            //     7**********
            //     6**********
            //     5**********
            //     4**********
            //     3**********
            //     2**********
            //     1**********
            //     0**********
            //     0123456789 Code12
            lat0  = code12 * 2 / 3;  
            long0 = code34 + 100*z;
            lat0  = lat0  + (code5 * 2 / 3) / 8; 
            long0 = long0 +  code6 / 8;
            lat0  = lat0  + code7 * 2 / 3 / 8 / 10;
            long0 = long0 + code8 / 8 / 10;
            lat0  = lat0  + (codeex9 * 2 / 3) / 8 / 10 / 10;
            long0 = long0 + codeex10 / 8 / 10 / 10;
            lat0  = lat0  + ((codeex11-x+1) * 2 / 3) / 8 / 10 / 10 / 10;
            long0 = long0 + (codeex12+y) / 8 / 10 / 10 / 10;
            lat0 = (1-2*x)*lat0;
            long0 = (1-2*y)*long0;
            dlat = 2/3/8/10/10/10;
            dlong = 1/8/10/10/10;
        }
        lat1 = myformat12(lat0-dlat);
        long1 = myformat12(long0+dlong);
        lat0 = myformat12(lat0);  
        long0 = myformat12(long0);  
    break;
    case 16:
        if(!extension){ 
        }else{ // Extended 1m grid square code (16digit)
            //      Code13
            //      9**********
            //      8**********
            //      7**********
            //      6**********
            //      5**********
            //      4**********
            //      3**********
            //      2**********
            //      1**********
            //      0**********
            //      0123456789 Code14
            lat0  = code12 * 2 / 3;
            long0 = code34 + 100*z;
            lat0  = lat0  + (code5 * 2 / 3) / 8;
            long0 = long0 +  code6 / 8;
            lat0  = lat0  + (code7 * 2 / 3) / 8 / 10;
            long0 = long0 + code8 / 8 / 10;
            lat0  = lat0  + (codeex9 * 2 / 3) / 8 / 10 / 10;
            long0 = long0 + codeex10 / 8 / 10 / 10;
            lat0  = lat0  + (codeex11 * 2 / 3) / 8 / 10 / 10 / 10;
            long0 = long0 + codeex12 / 8 / 10 / 10 / 10;
            lat0  = lat0  + ((codeex13-x+1) * 2 / 3) / 8 / 10 / 10 / 10 / 10;
            long0 = long0 + (codeex14+y) / 8 / 10 / 10 / 10 / 10;
            lat0 = (1-2*x)*lat0;
            long0 = (1-2*y)*long0;
            dlat = 2/3/8/10/10/10/10;
            dlong = 1/8/10/10/10/10;           
        }
        lat1 = myformat14(lat0-dlat);
        long1 = myformat14(long0+dlong);
        lat0 = myformat14(lat0);  
        long0 = myformat14(long0);  
    break;
    default:
        lat1 = '99999';
        long1 = '99999';
        lat0 = '99999';
        long0 = '99999';
    }
    //
    var xx = new Object();
    xx['lat0'] = parseFloat(lat0);
    xx['long0'] = parseFloat(long0);
    xx['lat1'] = parseFloat(lat1);
    xx['long1'] = parseFloat(long1);
    return(xx);
}

function myformat(v,cnt){
   var s = v+'';
   var ss;
   if(v > 100) ss = s.substr(0,3+cnt);
   else if(v > 10) ss = s.substr(0,2+cnt);
   else ss = s.substr(0,1+cnt);
   return(ss);
}

function myformat8(v){
    return myformat(v,8);
}

function myformat10(v){
    return myformat(v,10);
}

function myformat12(v){
    return myformat(v,12);
}

function myformat14(v){
    return myformat(v,14);
}

function meshcode_to_latlong(meshcode, extension=false){
    var res = meshcode_to_latlong_grid(meshcode, extension);
    var xx = new Object;
    xx['lat'] = res['lat0'];
    xx['long'] = res['long0'];
    return(xx);
}

function meshcode_to_latlong_NW(meshcode, extension=false){
    var res = meshcode_to_latlong_grid(meshcode, extension);
    var xx = new Object;
    xx['lat'] = res['lat0'];
    xx['long'] = res['long0'];
    return(xx);
}

function meshcode_to_latlong_SW(meshcode, extension=false){
    var res = meshcode_to_latlong_grid(meshcode, extension);
    var xx = new Object;
    xx['lat'] = res['lat1'];
    xx['long'] = res['long0'];
    return(xx);
}

function meshcode_to_latlong_NE(meshcode, extension=false){
    var res = meshcode_to_latlong_grid(meshcode, extension);
    var xx = new Object;
    xx['lat'] = res['lat0'];
    xx['long'] = res['long1'];
    return(xx);
}

function meshcode_to_latlong_SE(meshcode, extension=false){
    var res = meshcode_to_latlong_grid(meshcode, extension);
    var xx = new Object;
    xx['lat'] = res['lat1'];
    xx['long'] = res['long1'];
    return(xx);
}

function cal_meshcode6(latitude, longitude){
    var mesh;
    if(latitude < 0){
          var o = 4;
    }
    else{
          var o = 0;
    }
    if(longitude < 0){
          var o = o + 2;
    }
    if(Math.abs(longitude) >= 100) o = o + 1;
    var z = o % 2;
    var y = ((o-z)/2) % 2;
    var x = (o - 2*y - z)/4;
    var o = o + 1;
    latitude = (1-2*x)*latitude;
    longitude = (1-2*y)*longitude;
    var p = Math.floor(latitude*60/40);
    var a = latitude*60-p*40;
    var q = Math.floor(a/5);
    var b = a-q*5;
    var r = Math.floor(b*60/30);
    var c = b*60-r*30;
    var s2u = Math.floor(c/15);
    var d = c-s2u*15;
    var s4u = Math.floor(d/7.5);
    var e = d-s4u*7.5;
    var s8u = Math.floor(e/3.75);
    var u = Math.floor(longitude-100*z);
    var f = longitude-100*z-u;
    var v = Math.floor(f*60/7.5);
    var g = f*60-v*7.5;
    var w = Math.floor(g*60/45);
    var h = g*60-w*45;
    var s2l = Math.floor(h/22.5);
    var i = h-s2l*22.5;
    var s4l = Math.floor(i/11.25);
    var j = i-s4l*11.25;
    var s8l = Math.floor(j/5.625);
    var s2 = s2u*2+s2l+1;
    var s4 = s4u*2+s4l+1;
    var s8 = s8u*2+s8l+1;
    if(u < 10){
       if(p < 10){
           var mesh = String(o)+'00'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
       }else{
           if(p < 100){
               var mesh = String(o)+'0'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
           }
           else{
               var mesh = String(o)+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
           }
       }
    }else{
       if(p < 10){
            var mesh = String(o)+'00'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
       }else{
           if(p < 100){
                var mesh = String(o)+'0'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
           }else{
                var mesh = String(o)+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(s4)+String(s8);
           }
       }
    }
    return(mesh);
}
// 

function cal_meshcode(latitude,longitude){
    return(cal_meshcode3(latitude,longitude));
}

function cal_meshcode1(latitude,longitude){
    var mesh = cal_meshcode6(latitude,longitude);
    return(mesh.substr(0,6));
}

function cal_meshcode2(latitude,longitude){
    var mesh = cal_meshcode6(latitude,longitude);
    return(mesh.substr(0,8));
}

function cal_meshcode3(latitude,longitude){
    var mesh = cal_meshcode6(latitude,longitude);
    return(mesh.substr(0,10));
}

function cal_meshcode4(latitude,longitude){
    var mesh = cal_meshcode6(latitude,longitude);
    return(mesh.substr(0,11));
}

function cal_meshcode5(latitude,longitude){
    var mesh = cal_meshcode6(latitude,longitude);
    return(mesh.substr(0,12));
}

// calculate an extended 1m (10m (14digits) / 10) grid square code (16 digits) 
// from a geographical position (latitude, longitude) 
// - 0.03 arc-second for latitude and 0.045 arc-second for longitude
function cal_meshcode_ex1m_16(latitude, longitude){
  if(latitude < 0){
    var o = 4;
  }else{
    var o = 0;
  }
  if(longitude < 0){
    var o = o + 2;
  }
  if(Math.abs(longitude) >= 100) o = o + 1;
  var z = o % 2;
  var y = ((o - z)/2) % 2;
  var x = (o - 2*y - z)/4;
  var o = o + 1;
//
  latitude = (1-2*x)*latitude;
  longitude = (1-2*y)*longitude;
//
  var p = Math.floor(latitude*60/40);
  var a = latitude*60-p*40;
  var q = Math.floor(a/5);
  var b = a-q*5;
  var r = Math.floor(b*60/30);
  var c = b*60-r*30;
  var s = Math.floor(c/3);
  var d = c-s*3;
  var t = Math.floor(d/0.3);
  var e = d-t*0.3;
  var tt = Math.floor(e/0.03);
  var ee = e-tt*0.03;
  var u = Math.floor(longitude-100*z);
  var f = longitude-100*z-u;
  var v = Math.floor(f*60/7.5);
  var g = f*60-v*7.5;
  var w = Math.floor(g*60/45);
  var h = g*60-w*45;
  var xx = Math.floor(h/4.5);
  var i = h-xx*4.5;
  var yy = Math.floor(i/0.45);
  var j = i-yy*0.45;
  var zz = Math.floor(j/0.045);
  var k = j-zz*0.045;
  //
  if(u < 10){
    if(p < 10){
        var mesh = String(o)+'00'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
    }else{
        if(p < 100){
            var mesh = String(o)+'0'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
        }
        else{
            var mesh = String(o)+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
        }
    }
  }else{
    if(p < 10){
        var mesh = String(o)+'00'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
    }else{
        if(p < 100){
             var mesh = String(o)+'0'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
        }else{
            var mesh = String(o)+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s)+String(xx)+String(t)+String(yy)+String(tt)+String(zz);
        }
    }
  }
  return(mesh)
}

// calculate an extended 100m (1km / 10) grid square code (12 digits) 
// from a geographical position (latitude, longitude) 
// - 3 arc-second for latitude and 4.5 arc-second for longitude
function cal_meshcode_ex100m_12(latitude,longitude){
    var mesh = cal_meshcode_ex1m_16(latitude,longitude);
    return(mesh.substr(0,12));
}

// calculate an extended 10m (100m (12digits) / 10) grid square code (14 digits) 
// from a geographical position (latitude, longitude) 
// - 0.3 arc-second for latitude and 0.45 arc-second for longitude
function cal_meshcode_ex10m_14(latitude,longitude){
    var mesh = cal_meshcode_ex1m_16(latitude,longitude);
    return(mesh.substr(0,14));
}

// calculate an extended 100m grid square code (13 digits) 
// from a geographical position (latitude, longitude) 
// - 3 arc-second for latitude and 4.5 arc-second for longitude
function cal_meshcode_ex100m_13(latitude, longitude){
    if(latitude < 0){
        var o = 4;
    }
    else{
	var o = 0;
    }
    if(longitude < 0){
	var o = o + 2;
    }
    if(Math.abs(longitude) >= 100) o = o + 1;
    var z = o % 2;
    var y = ((o - z)/2) % 2;
    var x = (o - 2*y - z)/4;
//
    var o = o + 1;
//
    latitude = (1-2*x)*latitude;
    longitude = (1-2*y)*longitude;
//
    var p = Math.floor(latitude*60/40);
    var a = latitude*60-p*40;
    var q = Math.floor(a/5);
    var b = a-q*5;
    var r = Math.floor(b*60/30);
    var c = b*60-r*30;
    var s2u = Math.floor(c/15);
    var d = c-s2u*15;
    var et = Math.floor(d/3);
    var e = d-et*3;
    var u = Math.floor(longitude-100*z);
    var f = longitude-100*z-u;
    var v = Math.floor(f*60/7.5);
    var g = f*60-v*7.5;
    var w = Math.floor(g*60/45);
    var h = g*60-w*45;
    var s2l = Math.floor(h/22.5);
    var i = h-s2l*22.5;
    var jt = Math.floor(i/4.5);
    var j = i*0.5-jt*4.5;
    var s2 = s2u*2+s2l+1;
//
    if(u < 10){
	    if(p < 10){
	        var mesh = String(o)+'00'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	    }else{
	        if(p < 100){
		        var mesh = String(o)+'0'+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	        }else{
		        var mesh = String(o)+String(p)+'0'+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	        }
	    }
    }else{
	    if(p < 10){
	        var mesh = String(o)+'00'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	    }else{
	        if(p < 100){
		        var mesh = String(o)+'0'+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	        }else{
		        var mesh = String(o)+String(p)+String(u)+String(q)+String(v)+String(r)+String(w)+String(s2)+String(et)+String(jt);
	        }
	    }
    }
    return(mesh)
}

function Vincenty(latitude1,longitude1,latitude2,longitude2){
    // WGS84
    var f = 1/298.257223563;
    var a = 6378137.0;
    var b = 6356752.314245;
    //
    if(latitude1 == latitude2 && longitude1 == longitude2) return(0);
    var L = (longitude1 - longitude2)/180*Math.PI;
    var U1 = Math.atan((1-f)*Math.tan(latitude1/180*Math.PI));
    var U2 = Math.atan((1-f)*Math.tan(latitude2/180*Math.PI));
    var lambda = L;
    var dlambda = 10;
    while(Math.abs(dlambda) > 1e-12){
       var cs = Math.cos(U2)*Math.sin(lambda);
       var cscc = Math.cos(U1)*Math.sin(U2)-Math.sin(U1)*Math.cos(U2)*Math.cos(lambda);
       var sinsigma = Math.sqrt(cs*cs + cscc*cscc);
       var cossigma = Math.sin(U1)*Math.sin(U2)+Math.cos(U1)*Math.cos(U2)*Math.cos(lambda);
       var sigma = Math.atan(sinsigma/cossigma);
       var sinalpha = Math.cos(U1)*Math.cos(U2)*Math.sin(lambda)/sinsigma;
       var cos2alpha = 1 - sinalpha*sinalpha;
       if(cos2alpha == 0.0){
         var C = 0.0;
         var lambda0 = L + f*sinalpha*sigma;
       }else{
         var cos2sigmam = cossigma - 2*Math.sin(U1)*Math.sin(U2)/cos2alpha;
         var C = f/16*cos2alpha*(4+f*(4-3*cos2alpha));
         var lambda0 = L + (1-C)*f*sinalpha*(sigma + C*sinsigma*(cos2sigmam + C*cossigma*(-1+2*cos2sigmam*cos2sigmam)));
       }
       dlambda = lambda0 - lambda;
       lambda = lambda0;
    }
    if(C == 0.0){
      var A = 1.0;
      var dsigma = 0.0;
    }else{
      var u2 = cos2alpha * (a*a-b*b)/(b*b);
      var A = 1 + u2/16384*(4096 + u2 * (-768 + u2*(320-175*u2)));
      var B = u2/1024*(256+u2*(-128+u2*(74-47*u2)));
      var dsigma = B*sinsigma*(cos2sigmam + 1/4*B*(cossigma*(-1+2*cos2sigmam*cos2sigmam)-1/6*B*cos2sigmam*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmam*cos2sigmam)));
    }
    var s = b*A*(sigma-dsigma);
    return(s);
}

function cal_area_from_meshcode(meshcode,extension=false){
    var latlong = meshcode_to_latlong_grid(meshcode,extension);
    return(cal_area_from_latlong(latlong));
}

function cal_area_from_latlong(latlong){
    var W1 = Vincenty(latlong['lat0'],latlong['long0'],latlong['lat0'],latlong['long1']);
    var W2 = Vincenty(latlong['lat1'],latlong['long0'],latlong['lat1'],latlong['long1']);
    var H = Vincenty(latlong['lat0'],latlong['long0'],latlong['lat1'],latlong['long0']);
    var A = (W1+W2)*H*0.5;
    var xx = new Object;
    xx['W1'] = W1;
    xx['W2'] = W2;
    xx['H'] = H;
    xx['A'] = A;
    return(xx);
}

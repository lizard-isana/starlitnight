Ephemeris = function(){}

Ephemeris.prototype.JulianDay = function(date){
  d = new Date();
  d.setTime(date.getTime())
  var year=d.getUTCFullYear();
  var month=d.getUTCMonth()+1;
  var day=d.getUTCDate();
  var calender="";

  if(month <= 2){
    var year = year - 1;
    var month = month + 12;
  } 

  var julian_day = Math.floor(365.25*(year+4716))+Math.floor(30.6001*(month+1))+day-1524.5;
  
  if (calender == "julian"){ 
    var transition_offset=0;
  }else if(calender == "gregorian"){
    var tmp = Math.floor(year/100);
    var transition_offset=2-tmp+Math.floor(tmp/4);
  }else if(julian_day<2299160.5){
    var transition_offset=0;
  }else{
    var tmp = Math.floor(year/100);
    var transition_offset=2-tmp+Math.floor(tmp/4);
  }
    var jd=julian_day+transition_offset;
    this.jd = jd;
    return jd
}


Ephemeris.prototype.GMST = function(date){
//load default values
  d = new Date();
  d.setTime(date.getTime())
  var hours=d.getUTCHours();
  var minutes=d.getUTCMinutes();
  var seconds=d.getUTCSeconds();
  var time_in_sec=hours*3600+minutes*60+seconds;
  var jd = this.JulianDay(date);
  var rad=Math.PI/180;
  
//gmst at 0:00
  var t = (jd-2451545.0)/36525;
  var gmst_at_zero = (24110.5484 + 8640184.812866*t+0.093104*t*t+0.0000062*t*t*t)/3600;
  if(gmst_at_zero>24){gmst_at_zero=gmst_at_zero%24;}
  this.gmst_at_zero=gmst_at_zero;

//gmst at target time

  var gmst=gmst_at_zero+(time_in_sec * 1.00273790925)/3600;

  //mean obliquity of the ecliptic
  e = 23+26.0/60+21.448/3600 -46.8150/3600*t -0.00059/3600*t*t +0.001813/3600*t*t*t;
  //nutation in longitude
  omega = 125.04452-1934.136261*t+0.0020708*t*t+t*t*t/450000;
  long1 = 280.4665 + 36000.7698*t;
  long2 = 218.3165 + 481267.8813*t;
  phai = -17.20*Math.sin(omega*rad)-(-1.32*Math.sin(2*long1*rad))-0.23*Math.sin(2*long2*rad) + 0.21*Math.sin(2*omega*rad);
  gmst =gmst + ((phai/15)*(Math.cos(e*rad)))/3600

  if(gmst<0){gmst=gmst%24+24;}
  if(gmst>24){gmst=gmst%24;}
  this.gmst=gmst
  return gmst
}

Ephemeris.prototype.Sun = function(date){
  //load default values
  d = new Date();
  d.setTime(date.getTime())
  var hours=d.getUTCHours();
  var minutes=d.getUTCMinutes();
  var seconds=d.getUTCSeconds();
  var time_in_day=hours/24+minutes/1440+seconds/86400;
  
  var jd = this.JulianDay(date);
  var rad=Math.PI/180;

  //ephemeris days from the epch J2000.0
  var t = (jd + time_in_day -2451545.0)/36525;

  //geometric_mean_longitude
  var mean_longitude = 280.46646 + 36000.76983*t + 0.0003032*t*t;
  if(mean_longitude<0){mean_longitude=mean_longitude%360+360;}
  if(mean_longitude>360){mean_longitude=mean_longitude%360;}

  //mean anomaly of the Sun
  var mean_anomaly =  357.52911+ 35999.05029*t - 0.0001537*t*t;
  if(mean_anomaly<0){mean_anomaly=mean_anomaly%360+360;}
  if(mean_anomaly>360){mean_anomaly=mean_anomaly%360;}

  //eccentricity of the Earth's orbit
  var eccentricity = 0.016708634 - 0.000042037*t - 0.0000001267*t*t;

  //Sun's equation of  the center
  var equation = (1.914602 - 0.004817*t - 0.000014*t*t)*Math.sin(mean_anomaly*rad);
  equation += (0.019993 - 0.000101*t)*Math.sin(2*mean_anomaly*rad);
  equation += 0.000289 *Math.sin(3*mean_anomaly*rad);

  //true longitude of the Sun
  var true_longitude = mean_longitude + equation;
  //for more accuracy
  //var year=d.getUTCFullYear();
  //true_longitude = true_longitude - 0.01397*(year - 2000);  
  
  //true anomary of the Sun
  var true_anomary = mean_anomaly + equation;

  //radius vector, distance between center of the Sun and the Earth
  var radius = (1.000001018*(1-eccentricity*eccentricity))/(1 + eccentricity*Math.cos(true_anomary*rad));
  this.radius=radius;
  var distance=radius*149597870
  this.distance = distance;
  //apparent longitude
  var w = 125.04 - 1934.136*t;
  var apparent_longitude = true_longitude - 0.00569 - 0.00478 * Math.sin(w*rad);
  var longitude = apparent_longitude;
  this.longitude = longitude;

  //obliquity of the ecliptic
  var obliquity = 23+26.0/60+21.448/3600 -(46.8150/3600)*t -(0.00059/3600)*t*t +(0.001813/3600)*t*t*t; 

  //correction for apperent position of the sun 
  obliquity = obliquity+0.00256*Math.cos(w*rad);

  //right asantion of the Sun
  var sun_ra = Math.atan2(Math.cos(obliquity*rad)*Math.sin(longitude*rad), Math.cos(longitude*rad))
  sun_ra = sun_ra/rad;
  if(sun_ra<0){sun_ra=sun_ra+360;}
  if(sun_ra>360){sun_ra=sun_ra%360;}
  this.ra=sun_ra/15;

  //declination of the Sun
  var sun_dec = Math.asin(Math.sin(obliquity*rad)*Math.sin(longitude*rad));
  this.dec=sun_dec/rad;

	this.x = distance*Math.cos(longitude*rad);
	this.y = distance*(Math.sin(longitude*rad)*Math.cos(obliquity*rad));
	this.z = distance*(Math.sin(longitude*rad)*Math.sin(obliquity*rad));
  return this;
}

Ephemeris.prototype.Moon = function(date){
  d = new Date();
  d.setTime(date.getTime())
  var hours=d.getUTCHours();
  var minutes=d.getUTCMinutes();
  var seconds=d.getUTCSeconds();
  var time_in_day=hours/24+minutes/1440+seconds/86400;
  var rad=Math.PI/180;
  var deg=180/Math.PI;
  var jd = this.JulianDay(date);


  //ephemeris days from the epch J2000.0
  var t = (jd + time_in_day -2451545.0)/36525;
  var e = 1- 0.002516*t - 0.0000074*t*t
  var L0 = (280.4665 + 36000.7698*t)*rad
  var L1 = (218.3164477 + 481267.88123421*t - 0.0015786*t*t + t*t*t/538841 - t*t*t*t/65194000)*rad
  var D0 = (297.8501921 + 445267.1114034*t - 0.0018819*t*t + t*t*t/545868 - t*t*t*t/113065000)*rad
  var M0 = (357.5291092 + 35999.0502909*t - 0.0001536*t*t + t*t*t/24490000)*rad
  var M1 = (134.9633964 + 477198.8675055*t + 0.0087414*t*t + t*t*t/69699 - t*t*t*t/14712000)*rad
  var F0 = (93.2720950 + 483202.0175233*t - 0.0036539 *t*t - t*t*t/3526000 + t*t*t*t/863310000)*rad
  var A1 = (119.75 + 131.849*t)*rad
  var A2 = (53.09 + 479264.290*t)*rad
  var A3 = (313.45 + 481266.484*t)*rad

  var sigma_l=0
sigma_l += 6288774*Math.sin((0*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_l += 1274027*Math.sin((2*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_l += 658314*Math.sin((2*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_l += 213618*Math.sin((0*D0)+(0*M0)+(2*M1)+(0*F0));
sigma_l += -185116*e*Math.sin((0*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_l += -114332*Math.sin((0*D0)+(0*M0)+(0*M1)+(2*F0));
sigma_l += 58793*Math.sin((2*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_l += 57066*e*Math.sin((2*D0)+(-1*M0)+(-1*M1)+(0*F0));
sigma_l += 53322*Math.sin((2*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_l += 45758*e*Math.sin((2*D0)+(-1*M0)+(0*M1)+(0*F0));
sigma_l += -40923*e*Math.sin((0*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_l += -34720*Math.sin((1*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_l += -30383*e*Math.sin((0*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_l += 15327*Math.sin((2*D0)+(0*M0)+(0*M1)+(-2*F0));
sigma_l += -12528*Math.sin((0*D0)+(0*M0)+(1*M1)+(2*F0));
sigma_l += 10980*Math.sin((0*D0)+(0*M0)+(1*M1)+(-2*F0));
sigma_l += 10675*Math.sin((4*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_l += 10034*Math.sin((0*D0)+(0*M0)+(3*M1)+(0*F0));
sigma_l += 8548*Math.sin((4*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_l += -7888*e*Math.sin((2*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_l += -6766*e*Math.sin((2*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_l += -5163*Math.sin((1*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_l += 4987*e*Math.sin((1*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_l += 4036*e*Math.sin((2*D0)+(-1*M0)+(1*M1)+(0*F0));
sigma_l += 3994*Math.sin((2*D0)+(0*M0)+(2*M1)+(0*F0));
sigma_l += 3861*Math.sin((4*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_l += 3665*Math.sin((2*D0)+(0*M0)+(-3*M1)+(0*F0));
sigma_l += -2689*e*Math.sin((0*D0)+(1*M0)+(-2*M1)+(0*F0));
sigma_l += -2602*Math.sin((2*D0)+(0*M0)+(-1*M1)+(2*F0));
sigma_l += 2390*e*Math.sin((2*D0)+(-1*M0)+(-2*M1)+(0*F0));
sigma_l += -2348*Math.sin((1*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_l += 2236*e*e*Math.sin((2*D0)+(-2*M0)+(0*M1)+(0*F0));
sigma_l += -2120*e*Math.sin((0*D0)+(1*M0)+(2*M1)+(0*F0));
sigma_l += -2069*e*e*Math.sin((0*D0)+(2*M0)+(0*M1)+(0*F0));
sigma_l += 2048*e*e*Math.sin((2*D0)+(-2*M0)+(-1*M1)+(0*F0));
sigma_l += -1773*Math.sin((2*D0)+(0*M0)+(1*M1)+(-2*F0));
sigma_l += -1595*Math.sin((2*D0)+(0*M0)+(0*M1)+(2*F0));
sigma_l += 1215*e*Math.sin((4*D0)+(-1*M0)+(-1*M1)+(0*F0));
sigma_l += -1110*Math.sin((0*D0)+(0*M0)+(2*M1)+(2*F0));
sigma_l += -892*Math.sin((3*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_l += -810*e*Math.sin((2*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_l += 759*e*Math.sin((4*D0)+(-1*M0)+(-2*M1)+(0*F0));
sigma_l += -713*e*e*Math.sin((0*D0)+(2*M0)+(-1*M1)+(0*F0));
sigma_l += -700*e*e*Math.sin((2*D0)+(2*M0)+(-1*M1)+(0*F0));
sigma_l += 691*e*Math.sin((2*D0)+(1*M0)+(-2*M1)+(0*F0));
sigma_l += 596*e*Math.sin((2*D0)+(-1*M0)+(0*M1)+(-2*F0));
sigma_l += 549*Math.sin((4*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_l += 537*Math.sin((0*D0)+(0*M0)+(4*M1)+(0*F0));
sigma_l += 520*e*Math.sin((4*D0)+(-1*M0)+(0*M1)+(0*F0));
sigma_l += -487*Math.sin((1*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_l += -399*e*Math.sin((2*D0)+(1*M0)+(0*M1)+(-2*F0));
sigma_l += -381*Math.sin((0*D0)+(0*M0)+(2*M1)+(-2*F0));
sigma_l += 351*e*Math.sin((1*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_l += -340*Math.sin((3*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_l += 330*Math.sin((4*D0)+(0*M0)+(-3*M1)+(0*F0));
sigma_l += 327*e*Math.sin((2*D0)+(-1*M0)+(2*M1)+(0*F0));
sigma_l += -323*e*e*Math.sin((0*D0)+(2*M0)+(1*M1)+(0*F0));
sigma_l += 299*e*Math.sin((1*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_l += 294*Math.sin((2*D0)+(0*M0)+(3*M1)+(0*F0));
sigma_l += 0*Math.sin((2*D0)+(0*M0)+(-1*M1)+(-2*F0));

  var sigma_r=0
sigma_r += -20905355*Math.cos((0*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_r += -3699111*Math.cos((2*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_r += -2955968*Math.cos((2*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_r += -569925*Math.cos((0*D0)+(0*M0)+(2*M1)+(0*F0));
sigma_r += 48888*e*Math.cos((0*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_r += -3149*Math.cos((0*D0)+(0*M0)+(0*M1)+(2*F0));
sigma_r += 246158*Math.cos((2*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_r += -152138*e*Math.cos((2*D0)+(-1*M0)+(-1*M1)+(0*F0));
sigma_r += -170733*Math.cos((2*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_r += -204586*e*Math.cos((2*D0)+(-1*M0)+(0*M1)+(0*F0));
sigma_r += -129620*e*Math.cos((0*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_r += 108743*Math.cos((1*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_r += 104755*e*Math.cos((0*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_r += 10321*Math.cos((2*D0)+(0*M0)+(0*M1)+(-2*F0));
sigma_r += 0*Math.cos((0*D0)+(0*M0)+(1*M1)+(2*F0));
sigma_r += 79661*Math.cos((0*D0)+(0*M0)+(1*M1)+(-2*F0));
sigma_r += -34782*Math.cos((4*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_r += -23210*Math.cos((0*D0)+(0*M0)+(3*M1)+(0*F0));
sigma_r += -21636*Math.cos((4*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_r += 24208*e*Math.cos((2*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_r += 30824*e*Math.cos((2*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_r += -8379*Math.cos((1*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_r += -16675*e*Math.cos((1*D0)+(1*M0)+(0*M1)+(0*F0));
sigma_r += -12831*e*Math.cos((2*D0)+(-1*M0)+(1*M1)+(0*F0));
sigma_r += -10445*Math.cos((2*D0)+(0*M0)+(2*M1)+(0*F0));
sigma_r += -11650*Math.cos((4*D0)+(0*M0)+(0*M1)+(0*F0));
sigma_r += 14403*Math.cos((2*D0)+(0*M0)+(-3*M1)+(0*F0));
sigma_r += -7003*e*Math.cos((0*D0)+(1*M0)+(-2*M1)+(0*F0));
sigma_r += 0*Math.cos((2*D0)+(0*M0)+(-1*M1)+(2*F0));
sigma_r += 10056*e*Math.cos((2*D0)+(-1*M0)+(-2*M1)+(0*F0));
sigma_r += 6322*Math.cos((1*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_r += -9884*e*e*Math.cos((2*D0)+(-2*M0)+(0*M1)+(0*F0));
sigma_r += 5751*e*Math.cos((0*D0)+(1*M0)+(2*M1)+(0*F0));
sigma_r += 0*e*e*Math.cos((0*D0)+(2*M0)+(0*M1)+(0*F0));
sigma_r += -4950*e*e*Math.cos((2*D0)+(-2*M0)+(-1*M1)+(0*F0));
sigma_r += 4130*Math.cos((2*D0)+(0*M0)+(1*M1)+(-2*F0));
sigma_r += 0*Math.cos((2*D0)+(0*M0)+(0*M1)+(2*F0));
sigma_r += -3958*e*Math.cos((4*D0)+(-1*M0)+(-1*M1)+(0*F0));
sigma_r += 0*Math.cos((0*D0)+(0*M0)+(2*M1)+(2*F0));
sigma_r += 3258*Math.cos((3*D0)+(0*M0)+(-1*M1)+(0*F0));
sigma_r += 2616*e*Math.cos((2*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_r += -1897*e*Math.cos((4*D0)+(-1*M0)+(-2*M1)+(0*F0));
sigma_r += -2117*e*e*Math.cos((0*D0)+(2*M0)+(-1*M1)+(0*F0));
sigma_r += 2354*e*e*Math.cos((2*D0)+(2*M0)+(-1*M1)+(0*F0));
sigma_r += 0*e*Math.cos((2*D0)+(1*M0)+(-2*M1)+(0*F0));
sigma_r += 0*e*Math.cos((2*D0)+(-1*M0)+(0*M1)+(-2*F0));
sigma_r += -1423*Math.cos((4*D0)+(0*M0)+(1*M1)+(0*F0));
sigma_r += -1117*Math.cos((0*D0)+(0*M0)+(4*M1)+(0*F0));
sigma_r += -1571*e*Math.cos((4*D0)+(-1*M0)+(0*M1)+(0*F0));
sigma_r += -1739*Math.cos((1*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_r += 0*e*Math.cos((2*D0)+(1*M0)+(0*M1)+(-2*F0));
sigma_r += -4421*Math.cos((0*D0)+(0*M0)+(2*M1)+(-2*F0));
sigma_r += 0*e*Math.cos((1*D0)+(1*M0)+(1*M1)+(0*F0));
sigma_r += 0*Math.cos((3*D0)+(0*M0)+(-2*M1)+(0*F0));
sigma_r += 0*Math.cos((4*D0)+(0*M0)+(-3*M1)+(0*F0));
sigma_r += 0*e*Math.cos((2*D0)+(-1*M0)+(2*M1)+(0*F0));
sigma_r += 1165*e*e*Math.cos((0*D0)+(2*M0)+(1*M1)+(0*F0));
sigma_r += 0*e*Math.cos((1*D0)+(1*M0)+(-1*M1)+(0*F0));
sigma_r += 0*Math.cos((2*D0)+(0*M0)+(3*M1)+(0*F0));
sigma_r += 8752*Math.cos((2*D0)+(0*M0)+(-1*M1)+(-2*F0));

  var sigma_b=0
sigma_b += 5128122*Math.sin((0*D0)+(0*M0)+(0*M1)+(1*F0));
sigma_b += 280602*Math.sin((0*D0)+(0*M0)+(1*M1)+(1*F0));
sigma_b += 277693*Math.sin((0*D0)+(0*M0)+(1*M1)+(-1*F0));
sigma_b += 173237*Math.sin((2*D0)+(0*M0)+(0*M1)+(-1*F0));
sigma_b += 55413*Math.sin((2*D0)+(0*M0)+(-1*M1)+(1*F0));
sigma_b += 46271*Math.sin((2*D0)+(0*M0)+(-1*M1)+(-1*F0));
sigma_b += 32573*Math.sin((2*D0)+(0*M0)+(0*M1)+(1*F0));
sigma_b += 17198*Math.sin((0*D0)+(0*M0)+(2*M1)+(1*F0));
sigma_b += 9266*Math.sin((2*D0)+(0*M0)+(1*M1)+(-1*F0));
sigma_b += 8822*Math.sin((0*D0)+(0*M0)+(2*M1)+(-1*F0));
sigma_b += 8216*e*Math.sin((2*D0)+(-1*M0)+(0*M1)+(-1*F0));
sigma_b += 4324*Math.sin((2*D0)+(0*M0)+(-2*M1)+(-1*F0));
sigma_b += 4200*Math.sin((2*D0)+(0*M0)+(1*M1)+(1*F0));
sigma_b += -3359*e*Math.sin((2*D0)+(1*M0)+(0*M1)+(-1*F0));
sigma_b += 2463*e*Math.sin((2*D0)+(-1*M0)+(-1*M1)+(1*F0));
sigma_b += 2211*e*Math.sin((2*D0)+(-1*M0)+(0*M1)+(1*F0));
sigma_b += 2065*e*Math.sin((2*D0)+(-1*M0)+(-1*M1)+(-1*F0));
sigma_b += -1870*e*Math.sin((0*D0)+(1*M0)+(-1*M1)+(-1*F0));
sigma_b += 1828*Math.sin((4*D0)+(0*M0)+(-1*M1)+(-1*F0));
sigma_b += -1794*e*Math.sin((0*D0)+(1*M0)+(0*M1)+(1*F0));
sigma_b += -1749*Math.sin((0*D0)+(0*M0)+(0*M1)+(3*F0));
sigma_b += -1565*e*Math.sin((0*D0)+(1*M0)+(-1*M1)+(1*F0));
sigma_b += -1491*Math.sin((1*D0)+(0*M0)+(0*M1)+(1*F0));
sigma_b += -1475*e*Math.sin((0*D0)+(1*M0)+(1*M1)+(1*F0));
sigma_b += -1410*e*Math.sin((0*D0)+(1*M0)+(1*M1)+(-1*F0));
sigma_b += -1344*e*Math.sin((0*D0)+(1*M0)+(0*M1)+(-1*F0));
sigma_b += -1335*Math.sin((1*D0)+(0*M0)+(0*M1)+(-1*F0));
sigma_b += 1107*Math.sin((0*D0)+(0*M0)+(3*M1)+(1*F0));
sigma_b += 1021*Math.sin((4*D0)+(0*M0)+(0*M1)+(-1*F0));
sigma_b += 833*Math.sin((4*D0)+(0*M0)+(-1*M1)+(1*F0));
sigma_b += 777*Math.sin((0*D0)+(0*M0)+(1*M1)+(-3*F0));
sigma_b += 671*Math.sin((4*D0)+(0*M0)+(-2*M1)+(1*F0));
sigma_b += 607*Math.sin((2*D0)+(0*M0)+(0*M1)+(-3*F0));
sigma_b += 596*Math.sin((2*D0)+(0*M0)+(2*M1)+(-1*F0));
sigma_b += 491*e*Math.sin((2*D0)+(-1*M0)+(1*M1)+(-1*F0));
sigma_b += -451*Math.sin((2*D0)+(0*M0)+(-2*M1)+(1*F0));
sigma_b += 439*Math.sin((0*D0)+(0*M0)+(3*M1)+(-1*F0));
sigma_b += 422*Math.sin((2*D0)+(0*M0)+(2*M1)+(1*F0));
sigma_b += 421*Math.sin((2*D0)+(0*M0)+(-3*M1)+(-1*F0));
sigma_b += -366*e*Math.sin((2*D0)+(1*M0)+(-1*M1)+(1*F0));
sigma_b += -351*e*Math.sin((2*D0)+(1*M0)+(0*M1)+(1*F0));
sigma_b += 331*Math.sin((4*D0)+(0*M0)+(0*M1)+(1*F0));
sigma_b += 315*e*Math.sin((2*D0)+(-1*M0)+(1*M1)+(1*F0));
sigma_b += 302*e*e*Math.sin((2*D0)+(-2*M0)+(0*M1)+(-1*F0));
sigma_b += -283*Math.sin((0*D0)+(0*M0)+(1*M1)+(3*F0));
sigma_b += -229*e*Math.sin((2*D0)+(1*M0)+(1*M1)+(-1*F0));
sigma_b += 223*e*Math.sin((1*D0)+(1*M0)+(0*M1)+(-1*F0));
sigma_b += 223*e*Math.sin((1*D0)+(1*M0)+(0*M1)+(1*F0));
sigma_b += -220*e*Math.sin((0*D0)+(1*M0)+(-2*M1)+(-1*F0));
sigma_b += -220*e*Math.sin((2*D0)+(1*M0)+(-1*M1)+(-1*F0));
sigma_b += -185*Math.sin((1*D0)+(0*M0)+(1*M1)+(1*F0));
sigma_b += 181*e*Math.sin((2*D0)+(-1*M0)+(-2*M1)+(-1*F0));
sigma_b += -177*e*Math.sin((0*D0)+(1*M0)+(2*M1)+(1*F0));
sigma_b += 176*Math.sin((4*D0)+(0*M0)+(-2*M1)+(-1*F0));
sigma_b += 166*e*Math.sin((4*D0)+(-1*M0)+(-1*M1)+(-1*F0));
sigma_b += -164*Math.sin((1*D0)+(0*M0)+(1*M1)+(-1*F0));
sigma_b += 132*Math.sin((4*D0)+(0*M0)+(1*M1)+(-1*F0));
sigma_b += -119*Math.sin((1*D0)+(0*M0)+(-1*M1)+(-1*F0));
sigma_b += 115*e*Math.sin((4*D0)+(-1*M0)+(0*M1)+(-1*F0));
sigma_b += 107*e*e*Math.sin((2*D0)+(-2*M0)+(0*M1)+(1*F0));

    sigma_l += 3958*Math.sin(A1)
    sigma_l += 1962*Math.sin(L1-F0)
    sigma_l += 318*Math.sin(A2)

    sigma_b += -2235*Math.sin(L1)
    sigma_b += 382*Math.sin(A3)
    sigma_b += 175*Math.sin(A1-F0)
    sigma_b += 175*Math.sin(A1+F0)
    sigma_b += 127*Math.sin(L1-M1)
    sigma_b += -115*Math.sin(L1+M1)


    var longitude = (L1/rad)%360  + sigma_l/1000000
    var latitude = sigma_b/1000000
    var distance = 385000.56 + sigma_r/1000

    var apparent_longitude = longitude + 0.004610

    var omega = (125.04452 - 1934.13261*t + 0.00220708*t*t + t*t*t/450000)*rad

    var obliquity_zero = 23+26.0/60+21.448/3600 -46.8150/3600*t -0.00059/3600*t*t +0.001813/3600*t*t*t; 
    var obliquity_delta = 9.20/3600*Math.cos(omega) + 0.57/3600*Math.cos(2*L0) +0.10/3600*Math.cos(2*L1) -0.09/3600*Math.cos(2*omega)
    var obliquity= obliquity_zero + obliquity_delta

    //var ra = (Math.atan2(Math.sin(apparent_longitude*rad)*Math.cos(obliquity*rad)-Math.tan(latitude*rad)*Math.sin(obliquity*rad),Math.cos(apparent_longitude*rad))/rad)/15
    //var dec = Math.asin(Math.sin(latitude*rad)*Math.cos(obliquity*rad) + Math.cos(latitude*rad)*Math.sin(obliquity*rad)*Math.sin(apparent_longitude*rad))/rad
    var ra = (Math.atan2(Math.sin(longitude*rad)*Math.cos(obliquity*rad)-Math.tan(latitude*rad)*Math.sin(obliquity*rad),Math.cos(longitude*rad))/rad)/15
    var dec = Math.asin(Math.sin(latitude*rad)*Math.cos(obliquity*rad) + Math.cos(latitude*rad)*Math.sin(obliquity*rad)*Math.sin(longitude*rad))/rad

    // equatiorial horizontal parallax
  var parallax = Math.asin(6378.14/distance)/rad
  this.t =t
  this.longitude = longitude;
  this.latitude = latitude;
  this.ra = ra;
  this.dec = dec;
  this.distance = distance;
  this.parallax = parallax;
  return this;
	
}

Ephemeris.prototype.MoonPhase = function(date){
var ephem = new Ephemeris();
var rad=Math.PI/180;
  var date_now = new Date();
  date_now.setTime(date.getTime())
  var hours=date_now.getUTCHours();
  var minutes=date_now.getUTCMinutes();
  var seconds=date_now.getUTCSeconds();
  var time_in_day=hours/24+minutes/1440+seconds/86400;
var jd = ephem.JulianDay(date) + time_in_day;
var yy = date_now.getYear();
if (yy < 2000) { yy += 1900;}
var date_first = new Date(yy, 0, 1, 0, 0, 0);
var date_last = new Date(yy, 11, 31, 11, 59, 59, 999);
var since_new_year = (date_now - date_first)/(date_last-date_first);
var year = yy+since_new_year;
var k = (year-2000) * 12.3685;
k = Math.floor(k);
var t = k/1236.85;
var t2 = t*t;
var t3 = t*t*t;
var t4 = t*t*t*t;
var jde0 = 2451550.09766 + 29.530588861*k + 0.00015437*t2 - 0.000000150*t3 + 0.00000000073*t4;

var e = 1-0.002516*t - 0.0000074*t2;
e = round_angle(e);
//Sun's mean anomary at the time;
var m0 = 2.5534 + 29.10535670*k - 0.0000014*t2 - 0.00000011*t3;
m0 = round_angle(m0);
//Moon's mean anomary at the time;
var m1 = 201.5643 + 385.81693528*k + 0.0107582*t2 + 0.00001238*t3 - 0.000000011*t4; 
m1 = round_angle(m1);
//Moon's argument of latitude
var f = 160.7108 + 390.67050284*k - 0.0016118*t2-0.00000227*t3 + 0.000000011*t4;
f = round_angle(f);
//Longitude of the ascending node of lunar orbit
var omega = 124.7746 -  1.56375588*k + 0.0020672*t2 + 0.00000215*t3;
omega = round_angle(omega);

var c1 = 0;
c1 = c1 - 0.40720 * Math.sin(m1*rad);
c1 = c1 + 0.17241 * e * Math.sin(m0*rad);
c1 = c1 + 0.01608 * Math.sin(2*m1*rad);
c1 = c1 + 0.01039 * Math.sin(2*f*rad);
c1 = c1 + 0.00739 * e * Math.sin((m1-m0)*rad);
c1 = c1 - 0.00514 * e * Math.sin((m1+m0)*rad);
c1 = c1 + 0.00208 * e * e * Math.sin(2*m0*rad); 
c1 = c1 - 0.00111 * Math.sin((m1-2*f)*rad)
c1 = c1 - 0.00057 * Math.sin((m1+2*f)*rad)
c1 = c1 + 0.00056 * e * Math.sin((2*m1+m0)*rad);
c1 = c1 - 0.00042 * Math.sin(3*m1*rad);
c1 = c1 + 0.00042 * e * Math.sin((m0+2*f)*rad)
c1 = c1 + 0.00038 * e * Math.sin((m0-2*f)*rad)
c1 = c1 - 0.00024 * e * Math.sin((2*m1-m0)*rad);
c1 = c1 - 0.00017 * Math.sin(omega*rad);
c1 = c1 - 0.00007 * Math.sin((m1+2*m0)*rad);
c1 = c1 + 0.00004 * Math.sin((2*m1-2*f)*rad);
c1 = c1 + 0.00004 * Math.sin(3*m0 *rad);
c1 = c1 + 0.00003 * Math.sin((m1+m0-2*f)*rad);
c1 = c1 + 0.00003 * Math.sin((2*m1+2*f)*rad);
c1 = c1 - 0.00003 * Math.sin((m1+m0+2*f)*rad);
c1 = c1 + 0.00003 * Math.sin((m1-m0+2*f)*rad);
c1 = c1 - 0.00002 * Math.sin((m1-m0-2*f)*rad);
c1 = c1 - 0.00002 * Math.sin((3*m1+m0)*rad);
c1 = c1 + 0.00002 * Math.sin(4*m1*rad);

var a1 = 299.77 + 0.107408*k-0.009173*t2;
var a2 = 251.88 + 0.016321*k;
var a3 = 251.83 + 26.651886*k;
var a4 = 349.42 + 36.412478 *k; 
var a5 =  84.66 + 18.206239*k;
var a6 =  141.74+53.303771*k;
var a7 =  207.14+2.453732*k;
var a8 =  154.84+7.306860*k;
var a9 =  34.52+27.261239*k;
var a10 =  207.19+0.121824*k;
var a11 =  291.34+1.844379*k;
var a12 =  161.72+24.198154*k;
var a13 =  239.56+25.513099*k;
var a14 =  331.55+3.592518*k;

var c2 = 0;
c2 = c2 + 0.000325 *Math.sin(a1*rad);
c2 = c2 + 0.000165 *Math.sin(a2*rad);
c2 = c2 + 0.000164 *Math.sin(a3*rad);
c2 = c2 + 0.000126 *Math.sin(a4*rad);
c2 = c2 + 0.000110 *Math.sin(a5*rad);
c2 = c2 + 0.000062 *Math.sin(a6*rad);
c2 = c2 + 0.000060 *Math.sin(a7*rad);
c2 = c2 + 0.000056 *Math.sin(a8*rad);
c2 = c2 + 0.000047 *Math.sin(a9*rad);
c2 = c2 + 0.000042 *Math.sin(a10*rad);
c2 = c2 + 0.000040 *Math.sin(a11*rad);
c2 = c2 + 0.000037 *Math.sin(a12*rad);
c2 = c2 + 0.000035 *Math.sin(a13*rad);
c2 = c2 + 0.000023 *Math.sin(a14*rad);

var jde = jde0 + c1 + c2;

var phase_of_the_moon = jd - jde;
return phase_of_the_moon;

}

function round_angle(angle){
  if(angle>360){
    angle= angle%360
  }else if(angle<0){
    angle= angle%360+360
  }else{
    angle = angle;
  }
return angle;
}



//to_eci
Ephemeris.prototype.to_eci = function(ra,dec,distance){
  var rad=Math.PI/180;
	var x = distance*Math.cos(dec*rad)*Math.cos(ra*rad);
	var y = distance*Math.cos(dec*rad)*Math.sin(ra*rad);
	var z = distance*Math.sin(dec*rad);
	this.x = x;
	this.y = y;
	this.z = z;	
  return this;
}

Ephemeris.prototype.to_horizontal = function(date,longitude,latitude,ra,dec){
  var rad=Math.PI/180;
  var dec = dec*rad;  
  var gmst = this.GMST(date);
  var hour_angle = (gmst*15 + longitude - (ra*15));
  var h = hour_angle*rad;
  var lat = latitude*rad;
  var azimuth = (Math.atan2(-Math.cos(dec)*Math.sin(h),Math.sin(dec)*Math.cos(lat)-Math.cos(dec)*Math.sin(lat)*Math.cos(h)))/rad;
  var altitude = (Math.asin(Math.sin(dec)*Math.sin(lat)+Math.cos(lat)*Math.cos(dec)*Math.cos(h)))/rad;
  if (azimuth<0){azimuth = azimuth%360 +360}
  this.azimuth = azimuth;
  this.altitude = altitude;
  this.hour_angle = hour_angle;
  return this;
}

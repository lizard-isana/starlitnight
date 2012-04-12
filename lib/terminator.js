// Sub script (terminator overlay) for GoogleSatTrack var 3.5 RC
// Copyright:: Copyright (C) 2009 Isana Kashiwai<isana at lizard-tail.com>
// License::  GPL
// Update:: 2009.03.24
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
  var ephemeris = new Ephemeris()
  var jd = ephemeris.JulianDay(date);
  var rad=Math.PI/180;
  var deg=180/Math.PI;
  
//gmst at 0:00
  var t = (jd-2451545.0)/36525;
  var gmst_at_zero = (24110.5484 + 8640184.812866*t+0.093104*t*t+0.0000062*t*t*t)/3600;
  if(gmst_at_zero>24){gmst_at_zero=gmst_at_zero%24;}
  this.gmst_at_zero=gmst_at_zero;

//gmst at target time
  var gmst=gmst_at_zero+(time_in_sec * 1.00273790925)/3600;
  if(gmst<0){gmst=gmst%24+24;}
  if(gmst>24){gmst=gmst%24;}
  this.gmst=gmst
  return gmst
}

Ephemeris.prototype.Sun = function(date){
  //load default values
  d = new Date();
  d.setTime(date.getTime())
  var year=d.getUTCFullYear();
  var hours=d.getUTCHours();
  var minutes=d.getUTCMinutes();
  var seconds=d.getUTCSeconds();
  var time_in_day=hours/24+minutes/1440+seconds/86400;
  
  var ephemeris = new Ephemeris()
  var jd = ephemeris.JulianDay(date);
  var rad=Math.PI/180;
  var deg=180/Math.PI;

  //ephemeris days from the epch J2000.0
  var t = (jd + time_in_day -2451545.0)/36525;

  //geometric_mean_longitude
  var mean_longitude = 280.46646 + 36000.76938*t + 0.0003032*t*t;
   if(mean_longitude<0){mean_longitude=mean_longitude%360+360;}
   if(mean_longitude>360){mean_longitude=mean_longitude%360;}

  //mean anomaly of the Sun
  var mean_anomaly =  357.52911+ 35999.05029*t + 0.0001537*t*t;
   if(mean_anomaly<0){mean_anomaly=mean_anomaly%360+360;}
   if(mean_anomaly>360){mean_anomaly=mean_anomaly%360;}

  //eccentricity of the Earth's orbit
  var eccentricity = 0.01678634 + 0.000042037*t + 0.0000001267*t*t;
  //Sun's equation of  the center
  c = (1.914602 - 0.004817*t + 0.000014*t*t)*Math.sin(mean_anomaly*rad);
  c =c+ (0.019993 - 0.000101*t)*Math.sin(2*mean_anomaly*rad);
  c =c+ 0.000289 *Math.sin(3*mean_anomaly*rad);
  //true longitude of the Sun
  var true_longitude = mean_longitude + c;
  var true_anomary = mean_anomaly + c;
  //radius vector, distance between center of the Sun and the Earth
  var radius = (1.000001018*(1-eccentricity*eccentricity))/(1 + eccentricity*Math.cos(true_anomary*rad));
  this.radius=radius;
  //apparent longitude
  w = 125.04 - 1934.136*t;
  apparent_longitude = true_longitude - 0.00569 - 0.00478 * Math.sin(w*rad);
  //obliquity of the ecliptio
  obliquity = 23+26.0/60+21.448/3600 -46.8150/3600*t -0.00059/3600*t*t +0.001813/3600*t*t*t; 
  //correction for apperent position of the sun 
  obliquity = obliquity+0.00256*Math.cos(w*rad);
  //right asantion of the Sun
  var ra = Math.atan2(Math.cos(obliquity*rad)*Math.sin(apparent_longitude*rad), Math.cos(apparent_longitude*rad))
  ra = ra *deg;
 if(ra<0){ra=ra+360;}
 if(ra>360){ra=ra%360;}
  this.ra=ra/15;
 //declination of the Sun
 var dec = Math.asin(Math.sin(obliquity*rad)*Math.sin(apparent_longitude*rad));
 this.dec=dec * deg;
 return this;
}

Terminator=function(date){
 this.date = date;
}

Terminator.prototype.generate = function(){
  date = this.date;
  var rad=Math.PI/180;
  var ephem = new Ephemeris();
  var sun = ephem.Sun(date);
  var sun_ra = sun.ra;
  var sun_dec = sun.dec;
  var gmst = ephem.GMST(date);
  var sun_long = -(gmst*15 - sun_ra*15);
  if (sun_long>360){sun_long=sun_long%360};
  if (sun_long<0){sun_long=sun_long%360+360};
  var sun_lat = sun_dec;

  if(sun_lat>5){
    var polar_night_flag = "south";
  }else if(sun_lat<-5){
    var polar_night_flag = "north";    
  }else{
    var polar_night_flag = "none";
  }

  var NightCircle = function(start,end,rev){
    var lat_array=[];
    var lng_array=[];
    for (i=start; i<= end; i+=1){
      var delta_lat = Math.asin(Math.cos(sun_lat*rad)*Math.sin(i*rad))/rad;
      if(Math.abs(delta_lat)<85){
        var x = -Math.cos(sun_long*rad)*Math.sin(sun_lat*rad)*Math.sin(i*rad)-Math.sin(sun_long*rad)*Math.cos(i*rad);
        var y = -Math.sin(sun_long*rad)*Math.sin(sun_lat*rad)*Math.sin(i*rad)+Math.cos(sun_long*rad)*Math.cos(i*rad);
        var delta_long = Math.atan2(y,x)/rad;
        if(delta_long>360){delta_long=delta_long%360}
        if(delta_long<0){delta_long=delta_long%360+360}
        lat_array.push(delta_lat);
        lng_array.push(delta_long);
      }
    }
    if(rev==true){
      this.lat_array= lat_array.reverse();
      this.lng_array= lng_array.reverse();
    }else{
      this.lat_array= lat_array;
      this.lng_array= lng_array;
    }
    return this;
  }
  var qtr=NightCircle(0,90,true);
  var latitude_array_east = qtr.lat_array;
  var longitude_array_east = qtr.lng_array;
  
  var qtr=NightCircle(270,360,true);
  latitude_array_east = latitude_array_east.concat(qtr.lat_array);
  longitude_array_east = longitude_array_east.concat(qtr.lng_array);
  
  var qtr=NightCircle(90,180,false);
  var latitude_array_west = qtr.lat_array;
  var longitude_array_west = qtr.lng_array;
  
  var qtr=NightCircle(180,270,false);
  latitude_array_west = latitude_array_west.concat(qtr.lat_array);
  longitude_array_west = longitude_array_west.concat(qtr.lng_array);

  if(sun_lat < 0){
    latitude_array_west.reverse();
    longitude_array_west.reverse();
  }else if(sun_lat > 0){
    latitude_array_east.reverse();
    longitude_array_east.reverse();
  }
  var latitude_array=latitude_array_east.concat(latitude_array_west);
  var longitude_array=longitude_array_east.concat(longitude_array_west);

  this.latitude_array=latitude_array;
  this.longitude_array=longitude_array;
  this.polar_night_flag = polar_night_flag;  
  return this;
}


Terminator.prototype.split = function(){
  var trm =this.generate();
  var lng_array = trm.longitude_array;
  var lat_array = trm.latitude_array;
  var array_length=lng_array.length;
  var start_lng = lng_array[0];
  var last_lng = lng_array[array_length-1];

  var split_size = 180; //2 parts
  //var split_size = 90; //4 parts
  //var split_size = 45; //8 parts
  //var split_size = 30; //12 parts
  //var split_size = 22.5; //16 parts
  //var split_size = 15; //24 parts
  //var split_size = 11.25; //32 parts
  //var split_size = 10; //36 parts
  //var split_size = 5; //72 parts
  //var split_size = 3; //120 parts
  if(start_lng>split_size){
    var first_split_point = Math.floor(start_lng/split_size)*split_size+split_size;
    if(first_split_point<=0){first_split_point=360}
  }else{
    var first_split_point = split_size;
  }
   var last_split_point = Math.floor(last_lng/split_size)*split_size+split_size;

  var split_latitude = function(lat1,lng1,lat2,lng2,split_lng){
    if(split_lng==360){
      var offset = 360
    }else{
      var offset=0
    }
    var  split_lat = lat1 + (split_lng- lng1)*((lat2-lat1)/((lng2+offset - lng1)));
    return split_lat;
  }

  var splitted_lng_array = [];
  var splitted_lat_array = [];
  var split_latlng_array = function(start,end,split_count,next){

    for (sp = start;sp<=end;sp=sp+split_size){
      splitted_lng_array[split_count]=[];
      splitted_lat_array[split_count]=[];
      for (i = next; i < array_length; i++) {
        if(sp>lng_array[i]&&sp-split_size<lng_array[i]){
          splitted_lng_array[split_count].push(lng_array[i])
          splitted_lat_array[split_count].push(lat_array[i])
        }else{
          var spl_lat = split_latitude(lat_array[i-1],lng_array[i-1],lat_array[i],lng_array[i],sp);
          splitted_lng_array[split_count].push(sp);
          splitted_lat_array[split_count].push(spl_lat);
          next = i;
          split_count++;
          break;
        }
      }
    }
    this.splitted_lng_array = splitted_lng_array;
    this.splitted_lat_array = splitted_lat_array;
    this.split_point = sp;
    this.split_count = split_count;
    this.next = next;
  }

  if(start_lng<last_lng){
    var sp = new split_latlng_array(first_split_point,last_split_point,0,0);
  }else{
    var sp = new split_latlng_array(first_split_point,360,0,0);
    sp = new split_latlng_array(split_size,last_split_point,sp.split_count,sp.next);
  }


  var splitted_lng_array = sp.splitted_lng_array;
  var splitted_lat_array = sp.splitted_lat_array;
  
 var sp_length = splitted_lng_array.length
 for(i=1;i<sp_length;i++){
   var last_lng = splitted_lng_array[i-1][splitted_lng_array[i-1].length-1];
   var last_lat = splitted_lat_array[i-1][splitted_lat_array[i-1].length-1];
   if(last_lng == 360){last_lng=0}
   splitted_lng_array[i].unshift(last_lng)
   splitted_lat_array[i].unshift(last_lat)
 }
  this.splitted_lng_array = splitted_lng_array;
  this.splitted_lat_array = splitted_lat_array;
  return this;

}

Terminator.prototype.show = function(map,boundary,night_shade){
  var split = this.split();
  var splitted_lng_array = split.splitted_lng_array;
  var splitted_lat_array = split.splitted_lat_array;
  var polar_night_flag = this.polar_night_flag;  

//Generate Porygon Array
  if(polar_night_flag=="north"){
    var lat_limit = 85;
  }else if(polar_night_flag=="south"){
    var lat_limit = -85;
  }

  var PorigonMaker = function(lng_array,lat_array,start_lat){
    var polygon = []
    var polyline = []
    var lng_array_length =lng_array.length;
    var lat_array_length =lat_array.length;
    var pref=new Pref();
    if(lng_array[0]>180){
      var true_lng = lng_array[0]-360;
    }else{
      var true_lng = lng_array[0];    
    }
    if(pref.move_mode=="fixed" && true_lng<0){
      lat_array.reverse();
      lng_array.reverse();
    }
    delete pref;
    if (lng_array_length>0){
      if(polar_night_flag!="none"){polygon.push(new GLatLng(lat_limit,lng_array[0]));}
      polygon.push(new GLatLng(start_lat,lng_array[0]));
      for (var i = 0; i < lat_array_length; i++) {
        polygon.push(new GLatLng(lat_array[i],lng_array[i]));    
        polyline.push(new GLatLng(lat_array[i],lng_array[i]));    
      }
      polygon.push(new GLatLng(start_lat,lng_array[lng_array_length-1]));
     
      if(polar_night_flag!="none"){polygon.push(new GLatLng(lat_limit,lng_array[lng_array_length-1]));}
      polygon.push(polygon[0]);
  
      if((lat_array[lat_array_length-1]==lat_array[lat_array_length-2])&&(lng_array[lng_array_length-1]!=lng_array[lng_array_length-2])){
          polyline.pop();  
      }  
      if((lat_array[0]==lat_array[1])&&(lng_array[0]!=lng_array[1])){
          polyline.shift();  
      }
    }
    this.polygon=polygon;
    this.polyline=polyline;
    return this;
  }
  
  var night_polygon_array = [];
  var night_polyline_array = [];
  var start_lat = splitted_lat_array[0][0];
  var al = splitted_lng_array.length;
  for(var j=0;j<al;j++){
    var p = new PorigonMaker(splitted_lng_array[j],splitted_lat_array[j],start_lat);
    night_polygon_array[j] = p.polygon;
    night_polyline_array[j] = p.polyline;
  }

  if (night_shade){
  var npl1 = night_polygon_array.length;
   for(var k1=0;k1<npl1;k1++){
      map.addOverlay(new GPolygon(night_polygon_array[k1], "#000000",1, 0, "#000000", 0.35));
    }
  }

  if (boundary){
    var npl2 = night_polyline_array.length;
    for(var k2=0;k2<npl2;k2++){
      map.addOverlay(new GPolyline(night_polyline_array[k2], "#000000", 2, 0.5));
    }
  }

}

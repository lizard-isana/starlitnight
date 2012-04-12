
# Orbital Propagation Model 'SGP4' for NORAD mean element sets
# Algorism from "Space Track Report No.3 : Models for Propagation of NORAD Elements Sets"
# update 2009.02.24

require 'date'
include Math

class Numeric
  # from degrees to radian 
   def to_rad
     self.to_f * Math::PI / 180.0
   end

   # from radian to degrees
   def to_deg
     self.to_f * 180.0 / Math::PI
   end
end

class TLE
  #decode two line elements to hash
  def decode(line1,line2)
    orbital_elements={
    "line_number_1" => line1[0..0],
    "catalog_no_1" =>  line1[2..6],
    "security_classification" =>  line1[7..7],
    "international_identification" => line1[10..16],
    "epoch_year" => line1[18..19],
    "epoch" => line1[20..31],
    "first_derivative_mean_motion" => line1[33..42],
    "second_derivative_mean_motion" => line1[44..51],
    "bstar" => line1[53..60],
    "ephemeris_type" =>  line1[62..62],
    "element_number" =>  line1[64..67],
    "check_sum_1" =>   line1[68..68],
    "line_number_2" =>   line2[0..0],
    "catalog_no_2" =>  line2[2..6],
    "inclination" => line2[8..15],
    "right_ascension" => line2[17..24],
    "eccentricity" => line2[26..32],
    "argument_of_perigee" => line2[34..41],
    "mean_anomaly" => line2[43..50],
    "mean_motion" => line2[52..62],
    "rev_number_at_epoch" => line2[63..67],
    "check_sum_2" =>   line1[68..68]
    }
    return self.trim(orbital_elements)
  end

  #trim orbital elements to proper format
  def trim(orbital_elements)
    orbital_elements["line_number_1"]=orbital_elements["line_number_1"].to_s.strip
    orbital_elements["catalog_no_1"]=orbital_elements["catalog_no_1"].to_s.strip
    orbital_elements["security_classification"]=orbital_elements["security_classification"].to_s.strip
    orbital_elements["international_identification"]=orbital_elements["international_identification"].to_s.strip
    epy = orbital_elements["epoch_year"]
    #epoch_year should be smaller than 2057.
    if epy.to_i<57 then
      orbital_elements["epoch_year"]=epy.to_i+2000
    else
      orbital_elements["epoch_year"]=epy.to_i+1900
    end
    orbital_elements["epoch"] = orbital_elements["epoch"].to_f
    orbital_elements["first_derivative_mean_motion"] = orbital_elements["first_derivative_mean_motion"].to_f
    second_derivative_mean_motion_mantissa = orbital_elements["second_derivative_mean_motion"][0..5].to_s.strip
    second_derivative_mean_motion_exponent = "1e" + orbital_elements["second_derivative_mean_motion"][6..7].to_s.strip
    orbital_elements["second_derivative_mean_motion"] = second_derivative_mean_motion_mantissa.to_f*second_derivative_mean_motion_exponent.to_f
    bstar_mantissa = (orbital_elements["bstar"][0..5]).to_i*1e-5
    bstar_exponent = "1e" + orbital_elements["bstar"][6..7].to_s.strip
    orbital_elements["bstar"] = bstar_mantissa.to_f*bstar_exponent.to_f
    orbital_elements["ephemeris_type"]=orbital_elements["ephemeris_type"].to_s.strip
    orbital_elements["element_number"]=orbital_elements["element_number"].to_i
    orbital_elements["check_sum_1"]=orbital_elements["check_sum_1"].to_s.strip
    
    orbital_elements["line_number_2"]=orbital_elements["line_number_2"].to_s.strip
    orbital_elements["catalog_no_2"]=orbital_elements["catalog_no_2"].to_s.strip
    orbital_elements["inclination"]=orbital_elements["inclination"].to_f
    orbital_elements["right_ascension"]=orbital_elements["right_ascension"].to_f
    orbital_elements["eccentricity"] = orbital_elements["eccentricity"].to_i*1e-7
    orbital_elements["argument_of_perigee"]=orbital_elements["argument_of_perigee"].to_f
    orbital_elements["mean_anomaly"]=orbital_elements["mean_anomaly"].to_f
    orbital_elements["mean_motion"]=orbital_elements["mean_motion"].to_f
    orbital_elements["rev_number_at_epoch"]=orbital_elements["rev_number_at_epoch"].to_i
    orbital_elements["check_sum_2"]=orbital_elements["check_sum_2"].to_s.strip

    return orbital_elements
  end
end

class SGP4

  def initialize(orbital_elements)
    @orbital_elements = orbital_elements
    @xincl = orbital_elements["inclination"].to_rad
    @xnodeo = orbital_elements["right_ascension"].to_rad
    @eo = orbital_elements["eccentricity"]
    @omegao  = orbital_elements["argument_of_perigee"].to_rad
    @xmo = orbital_elements["mean_anomaly"].to_rad
    @xno  = orbital_elements["mean_motion"]*2.0*PI/1440.0;
    @bstar = orbital_elements["bstar"]
    
    #set constant
    @xj2 = 1.082616e-3
    @xj4 = -1.65597e-6
    @xj3 = -0.253881e-5
    @xke = 0.743669161e-1
    @xkmper = 6378.135
    @ae = 1.0
    @qo = @ae +120.0/@xkmper
    @ck2 =5.413080e-4
    @ck4 = 0.62098875e-6
    @qoms2t = 1.88027916e-9
    @s = 1.0+78.0/@xkmper
    @de2ra = 0.174532925e-1
    @pio2 = 1.57079633
    @twopi = 6.2831853
    @x3pio2 = 4.71238898
    @x=@y=@z=@xdot=@ydot=@zdot=@perigee=@latitude=@longitude=@altitude=@velocity = 0
    
    a1 = (@xke/@xno)**(2.0/3.0)
    cosio = cos(@xincl)
    theta2=cosio**2
    x3thm1=3*theta2-1.0
    betao2=1-@eo**2
    betao=betao2**0.5
    del1=1.5*@ck2*x3thm1/(a1*a1*betao*betao2)
    ao=a1*(1-del1*((1.0/3.0)+del1*(1.0+(134.0/81.0)*del1)))
    delo=1.5*@ck2*x3thm1/(ao*ao*betao*betao2)
    xnodp=@xno/(1.0+delo) #original_mean_motion
    aodp=ao/(1.0-delo) #semi_major_axis

    isimp=0;
    if (aodp*(1.0-@eo)/@ae) < (220.0/@xkmper+@ae) then
      isimp=1;
    end

    s4=@s;
    qoms24=@qoms2t;
    perigee=(aodp*(1.0-@eo)-@ae)*@xkmper

    if perigee < 156.0 then
       s4 = perigee-78.0 
       if perigee <= 98.0 then
         s4 = 20.0
       else
         qoms24=((120.0-s4)*@ae/@xkmper)**4
         s4 = s4/@xkmper+@ae;
       end
    end
    pinvsq=1.0/(aodp*aodp*betao2*betao2);
    tsi=1.0/(aodp-s4);
    eta=aodp*@eo*tsi;
    etasq=eta*eta;
    eeta=@eo*eta;
    psisq=(1.0-etasq).abs
    coef=qoms24*(tsi**4)
    coef1=coef/(psisq**3.5);

    c2=coef1*xnodp*(aodp*(1.0+1.5*etasq+eeta*(4.0+etasq))+0.75*@ck2*tsi/psisq*x3thm1*(8.0+3.0*etasq*(8.0+etasq)))
    c1=@bstar*c2;
    sinio=sin(@xincl)
    a3ovk2=-@xj3/@ck2*(@ae**3)
    c3=coef*tsi*a3ovk2*xnodp*@ae*sinio/@eo;
    x1mth2=1.0-theta2;
    c4=2.0*xnodp*coef1*aodp*betao2*(eta*(2.0+0.5*etasq)+@eo*(0.5+2.0*etasq)-2.0*@ck2*tsi/(aodp*psisq)*(-3.0*x3thm1*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+0.75*x1mth2*(2.0*etasq-eeta*(1.0+etasq))*cos(2.0*@omegao)))
    c5=2.0*coef1*aodp*betao2*(1.0+2.75*(etasq+eeta)+eeta*etasq)
    theta4=theta2*theta2
    temp1=3.0*@ck2*pinvsq*xnodp
    temp2=temp1*@ck2*pinvsq
    temp3=1.25*@ck4*pinvsq*pinvsq*xnodp
    xmdot=xnodp+0.5*temp1*betao*x3thm1+0.0625*temp2*betao*(13.0-78.0*theta2+137.0*theta4);

    x1m5th=1.0-5.0*theta2;
    omgdot=-0.5*temp1*x1m5th+0.0625*temp2*(7.0-114.0*theta2+395.0*theta4)+temp3*(3.0-36.0*theta2+49.0*theta4);
    xhdot1=-temp1*cosio;
    xnodot=xhdot1+(0.5*temp2*(4.0-19.0*theta2)+2.0*temp3*(3.0-7.0*theta2))*cosio;
    omgcof=@bstar*c3*cos(@omegao);
    xmcof=-(2.0*3.0)*coef*@bstar*@ae/eeta;
    xnodcf=3.5*betao2*xhdot1*c1;
    t2cof=1.5*c1;
    xlcof=0.125*a3ovk2*sinio*(3.0+5.0*cosio)/(1.0+cosio);
    aycof=0.25*a3ovk2*sinio;
    delmo=(1.0+eta*Math.cos(@xmo))**3
    sinmo=sin(@xmo)
    x7thm1=7.0*theta2-1.0;

    if isimp != 1 then
      c1sq=c1*c1
      d2=4.0*aodp*tsi*c1sq
      temp=d2*tsi*c1/3.0
      d3=(17.0*aodp+s4)*temp
      d4=0.5*temp*aodp*tsi*(221.0*aodp+31.0*s4)*c1
      t3cof=d2+2.0*c1sq
      t4cof=0.25*(3.0*d3+c1*(12.0*d2+10.0*c1sq))
      t5cof=0.2*(3.0*d4+12.0*c1*d3+6.0*d2*d2+15.0*c1sq*(2.0*d2+c1sq))
     end

    @xmdot = xmdot;
    @omgdot = omgdot;
    @xnodot = xnodot;
    @xnodcf = xnodcf
    @t2cof = t2cof;
    @omgcof = omgcof;
    @isimp = isimp;
    @xmcof = xmcof;
    @eta = eta;
    @delmo = delmo;
    @c1 = c1;
    @c4 = c4;
    @c5 = c5;
    @d2 = d2;
    @d3 = d3;
    @d4 = d4;
    @sinmo = sinmo;
    @t3cof = t3cof;
    @t4cof = t4cof;
    @t5cof = t5cof;
    @aodp = aodp;
    @xnodp = xnodp;
    @xlcof = xlcof;
    @aycof = aycof;
    @x3thm1 = x3thm1;
    @x1mth2 = x1mth2;
    @cosio = cosio;
    @sinio = sinio;
    @x7thm1 = x7thm1;

  end

  def calc(date)

    xmdot = @xmdot;
    omgdot = @omgdot;
    xnodeo = @xnodeo;
    xnodot = @xnodot;
    xnodcf = @xnodcf
    t2cof = @t2cof;
    omgcof = @omgcof;
    isimp = @isimp;
    xmcof = @xmcof;
    eta = @eta;
    delmo = @delmo;
    c1 = @c1;
    c4 = @c4;
    c5 = @c5;
    d2 = @d2;
    d3 = @d3;
    d4 = @d4;
    sinmo = @sinmo;
    t3cof = @t3cof;
    t4cof = @t4cof;
    t5cof = @t5cof;
    aodp = @aodp;
    eo = @eo;
    xnodp = @xnodp;
    xke = @xke;
    xlcof = @xlcof;
    aycof = @aycof;
    x3thm1 = @x3thm1;
    x1mth2 = @x1mth2;
    xincl = @xincl;
    cosio = @cosio;
    sinio = @sinio;
    x7thm1 = @x7thm1;

    t = Clock.new(date)
    @tsince = t.elapsed_time(@orbital_elements["epoch_year"],@orbital_elements["epoch"])/60
    @gmst = t.gmst*15

    xmdf=@xmo+xmdot*@tsince
    omgadf=@omegao+omgdot*@tsince
    xnoddf=@xnodeo+xnodot*@tsince
    omega=omgadf
    xmp=xmdf
    tsq=@tsince*@tsince
    xnode=xnoddf+xnodcf*tsq
    tempa=1.0-c1*@tsince
    tempe=@bstar*c4*@tsince
    templ=t2cof*tsq



    if isimp != 1 then
      delomg=omgcof*@tsince
      delm=xmcof*(((1.0+eta*Math.cos(xmdf))**3)-delmo)
      temp=delomg+delm
      xmp=xmdf+temp
      omega=omgadf-temp
      tcube=tsq*@tsince
      tfour=@tsince*tcube
      tempa=tempa-d2*tsq-d3*tcube-d4*tfour
      tempe=tempe+@bstar*c5*(sin(xmp)-sinmo)
      templ=templ+t3cof*tcube+tfour*(t4cof+@tsince*t5cof)
    end
    a=aodp*tempa*tempa
    e=@eo-tempe
    xl=xmp+omega+xnode+xnodp*templ
    beta=(1.0-e*e)**0.5
    xn=@xke/(a**1.5)


  #long period periodics
    axn=e*cos(omega);
    temp=1.0/(a*beta*beta);
    xll=temp*xlcof*axn;
    aynl=temp*aycof;
    xlt=xl+xll;
    ayn=e*sin(omega)+aynl;

    #solve keplers equation
    capu = (xlt-xnode)%(2.0*PI);
    temp2=capu;
    for i in 1..10
      sinepw=sin(temp2);
      cosepw=cos(temp2);
      temp3=axn*sinepw;
      temp4=ayn*cosepw;
      temp5=axn*cosepw;
      temp6=ayn*sinepw;
      epw=(capu-temp4+temp3-temp2)/(1.0-temp5-temp6)+temp2;
      if (epw-temp2).abs <= 1.0e-6 then
        break
      end
      temp2=epw;
    end

    #short period preliminary quantities

    ecose=temp5+temp6;
    esine=temp3-temp4;
    elsq=axn*axn+ayn*ayn;
    temp=1.0-elsq;
    pl=a*temp;
    r=a*(1.0-ecose);
    temp1=1.0/r;
    rdot=@xke*(a**0.5)*esine*temp1;
    rfdot=@xke*(pl**0.5)*temp1;
    temp2=a*temp1;
    betal=temp**0.5;
    temp3=1.0/(1.0+betal);
    cosu=temp2*(cosepw-axn+ayn*esine*temp3);
    sinu=temp2*(sinepw-ayn-axn*esine*temp3);
    u=atan2(sinu,cosu);
    if u<0 then u+= 2*PI end
    sin2u=2.0*sinu*cosu;
    cos2u=2.0*cosu*cosu-1.0;
    temp=1.0/pl;
    temp1=@ck2*temp;
    temp2=temp1*temp;

    #update for short periodics

    rk=r*(1.0-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
    uk=u-0.25*temp2*x7thm1*sin2u;
    xnodek=xnode+1.5*temp2*cosio*sin2u;
    xinck=@xincl+1.5*temp2*cosio*sinio*cos2u;
    rdotk=rdot-xn*temp1*x1mth2*sin2u;
    rfdotk=rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);

    #orientation vectors

    sinuk=sin(uk);
    cosuk=cos(uk);
    sinik=sin(xinck);
    cosik=cos(xinck);
    sinnok=sin(xnodek);
    cosnok=cos(xnodek);
    xmx=-sinnok*cosik;
    xmy=cosnok*cosik;
    ux=xmx*sinuk+cosnok*cosuk;
    uy=xmy*sinuk+sinnok*cosuk;
    uz=sinik*sinuk
    vx=xmx*cosuk-cosnok*sinuk
    vy=xmy*cosuk-sinnok*sinuk
    vz=sinik*cosuk
    x=rk*ux
    y=rk*uy
    z=rk*uz
    xdot=rdotk*ux+rfdotk*vx
    ydot=rdotk*uy+rfdotk*vy
    zdot=rdotk*uz+rfdotk*vz

    xkm = x*@xkmper;
    ykm = y*@xkmper;
    zkm = z*@xkmper;
    xdotkmps = xdot*@xkmper/60;
    ydotkmps = ydot*@xkmper/60;
    zdotkmps = zdot*@xkmper/60;


    #Convert Earth-Centered Inertial (ECI) Coordinate System to Equatorial coordinate system
    velocity = ((xdot**2 + ydot**2 + zdot**2)**0.5)*@xkmper/60;

    us_axis=xkm*cos(@gmst.to_rad)+ykm*sin(@gmst.to_rad);
    vs_axis=-xkm*sin(@gmst.to_rad)+ykm*cos(@gmst.to_rad);
    ws_axis=zkm;

    #longitude
    longitude=atan(vs_axis/us_axis).to_deg
    if us_axis<0 then
      longitude=longitude+180;
    end
    if longitude>=180 then
      longitude=longitude-360;
    end

    #latitude
    t_lat1=ws_axis/((us_axis*us_axis+vs_axis*vs_axis)**0.5)
    t_lat2=(6377.397155*0.006674372230614)/((us_axis*us_axis+vs_axis*vs_axis)**0.5)

    tan_lat=0
    for i in 1..10
      tmp=tan_lat;
      tan_lat= t_lat1+t_lat2*(tan_lat/((1+0.99332*tan_lat*tan_lat)**0.5))
      if (tan_lat-tmp).abs<0.00000001 then
        break
      end
    end
    latitude=atan(tan_lat).to_deg
    
    #altitude
    altitude=((1+tan_lat*tan_lat)**0.5)*(ws_axis/tan_lat-6377.397155*0.993325628/((1+0.993325628*tan_lat*tan_lat)**0.5));

    @x,@y,@z = xkm,ykm,zkm
    @xdot,@ydot,@zdot = xdotkmps,ydotkmps,zdotkmps
    @perigee = perigee
    @latitude,@longitude,@altitude,@velocity = latitude,longitude,altitude,velocity

    return latitude,longitude,altitude,velocity;

  end

  attr_accessor :x, :y, :z, :xdot, :ydot, :zdot, :perigee, :latitude, :longitude, :altitude, :velocity

end


class Clock
  def initialize(date)
      date = date.utc
      @date = date
      @year,@month,@day,@hour,@min,@sec = date.year,date.month,date.day,date.hour,date.min,date.sec
      @time_in_day = @hour.to_f/24 + @min.to_f/1140 + @sec.to_f/86400
      @time_in_hour = @hour.to_f + @min.to_f/60 + @sec.to_f/3600
      @time_in_sec = @hour*3600+@min*60+@sec
  end

  def jd
    
    if @month<=2 then
      y=(@year-1).to_i
      m=(@month+12).to_i
    else
     y = @year
     m = @month
    end
    
    julian_day = (365.25*(y+4716)).floor+(30.6001*(m+1)).floor+@day-1524.5

    if julian_day<2299160.5 then
      transition_offset=0
    else
      tmp = (@year/100).floor
      transition_offset=2-tmp+(tmp/4).floor  
    end
 
    return julian_day=julian_day+transition_offset
  end
 
  def gmst
    jd = jd()
    t = (jd-2451545.0)/36525
    gmst0 = (24110.54841+8640184.812866*t+0.093104*t*t+0.0000062*t*t*t)/3600
    gmst = gmst0 + (@time_in_sec * 1.00273790925)/3600

    if gmst<0 then gmst=gmst%24+24 end
    if gmst>24 then gmst=gmst%24 end
    return gmst

  end
 
  def elapsed_time(epoch_year,epoch)
    year2=epoch_year-1;
    #now_in_sec=Date.UTC(@year, @month, @day, @hour, @min, @sec);
    epoch_date=Time.utc(year2, 12, 31, 0, 0, 0)+(epoch*24*60*60);
    return elapsed_time = @date-epoch_date
  end

end

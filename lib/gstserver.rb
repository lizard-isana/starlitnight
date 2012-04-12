#!/usr/local/bin/ruby

#Server Script for GoogleSatTrack ver. 09
require 'time'
require 'date'
require 'parsedate'
require 'sgp4'
require 'cgi'
require 'net/http'
require 'rexml/document' 

include REXML  
include Math

class LoadTLE < TLE
  attr_reader :update, :satname, :first_line, :second_line
  attr_accessor :sat_name, :first_line, :second_line
  def initialize(target="iss")
    @target = target
    @epoch_expired = 3600
    @data_type ="TLE"
    @sat_name = ""
    @first_line = ""
    @second_line = ""
    case
    when @target == "sts" then
      @tle_file = "./data/sts_tle.xml"
    when @target == "atv" then
      @tle_file = "./data/atv_tle.xml"
    when @target == "htv" then
      @tle_file = "./data/htv_tle.xml"
    when @target == "soyuz" then
      @tle_file = "./data/soyuz_tle.xml"
    when @target == "hst" then
      @tle_file = "./data/hst_tle.xml"

    else
      @tle_file = "./data/iss_tle.xml"
    end
  end
  
  def load
    if @target == "custom" then
     tle = custom_tle
    else
      local_tle = get_local
      expire = expire?
      if expire == true then
        tle = get_remote
        save
      else
        tle = local_tle
      end
    end
    return tle
  end

  def custom_tle
      @tle_from = "custom" 
  end

  def get_local
    #puts "access to local file"
    file_exist_flag = FileTest.exist?(@tle_file)
    if(file_exist_flag) then
      file = open(@tle_file,"r")

      doc = Document.new file
      update=doc.elements['/tle/update'].text
      sat_name=doc.elements['/tle/sat_name'].text
      first_line=doc.elements['/tle/first_line'].text
      second_line=doc.elements['/tle/second_line'].text
      file.close
      @update = update
      @sat_name = sat_name
      @first_line = first_line
      @second_line = second_line
      @tle_from = "local"
    else
      tle = get_remote
      save
    end
  end
  
  def get_remote
    #puts "access to external server"
    if @target == "sts" || @target == "iss" then
      from_nasa
      #from_celestrak
    elsif @target == "atv" then
      from_celestrak
    elsif @target == "htv" then
      from_celestrak
    elsif @target == "soyuz" then
      from_celestrak
    elsif @target == "hst" then
      from_celestrak

    end
    @update = Time.now.utc
    @tle_from = "remote"
  end

  def from_nasa
     host='spaceflight.nasa.gov'
    if @target=="sts"
      @sat_name = 'Space Shuttle'
      dir_file='/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/SHUTTLE/SVPOST.html'
    else
      sat_name = 'International Space Station'
      dir_file='/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html'
    end
  begin
    http = Net::HTTP.new(host, 80)
    response = http.get(dir_file)
  rescue
    from_celestrak
  else
    data = response.body
    num_array=[]
    counter=0
    now = Time.now.utc
    data_array_full=[]
    data.each {|line| 
      if /TWO LINE MEAN ELEMENT SET/=~line then
        num_array.push(counter+3)
        num_array.push(counter+4)
      end
    counter+=1
    data_array_full << line
    }
    data_array=[]
    num_array.each{ |num|
      data_array << data_array_full[num].strip
    }    
    t=0
    while data_array.size > t
      epoch_year="20"+ data_array[t+2].slice(18,2)
      epoch_day_next=data_array[t+2].slice(20,12)
      epoch_day_carrent=data_array[t].slice(20,12)
      epoch=Time.gm(epoch_year.to_i-1, 12, 31, 0, 0, 0)
      epoch_next = epoch + epoch_day_next.to_f*86400
      
      if epoch_next-now>0 then
        @first_line=data_array[t]
        @second_line=data_array[t+1]
        break
      end
      t += 2
    end
  end # end rescue
  end

  def from_celestrak
    host='celestrak.com'
    case
    when @target=="iss"
       @sat_name = 'Internationl Space Station'
        dir_file='/NORAD/elements/stations.txt'
        key = "ISS"
    when @target=="sts"
       @sat_name = 'SpaceShuttle'
        dir_file='/NORAD/elements/stations.txt'
        key = "STS"
    when @target=="atv"
       @sat_name = 'ATV/Jules Verne'
        dir_file='/NORAD/elements/stations.txt'
        key = "ATV"
    when @target=="htv"
       @sat_name = 'HTV'
        dir_file='/NORAD/elements/stations.txt'
        key = "HTV"
    when @target=="hst"
       @sat_name = 'Hubble Space Telescope'
        dir_file='/NORAD/elements/science.txt'
        key = "HST"
    when @target=="soyuz"
       @sat_name = 'Soyuz-TMA'
        dir_file='/NORAD/elements/stations.txt'
        key = "SOYUZ-TMA"

    end

      http = Net::HTTP.new(host, 80)
      response = http.get(dir_file)
      data = response.body
      data_array_full=[]
      counter = 0
      line_num = 0
      data.each {|line|
        if /^#{key}/ =~ line then
          line_num = counter
        end
        counter += 1
        data_array_full << line
      }
      @first_line=data_array_full[line_num+1].strip
      @second_line=data_array_full[line_num+2].strip
  end

  def expire?
    now = Time.new.utc
    update = Time.parse(@update.to_s)
    if now - update > @epoch_expired  then
      return true
    else
      return false
    end
    
  end

  def save
    now = Time.now.utc
    xml = File.open(@tle_file,'w')
    xml.puts <<EOS
<?xml version="1.0" encoding="utf-8"?> 
<tle>
<update>#{now}</update>
<data_type>#{@data_type}</data_type>
<sat_name>#{@sat_name}</sat_name>
<first_line>#{@first_line}</first_line>
<second_line>#{@second_line}</second_line>
</tle>
EOS
    xml.close
  end

  def show
print "Content-type:text/xml\n"
print "Pragma: no-cache\n"
print "Cache-Control: no-cache\n\n"
puts <<EOS
<?xml version="1.0" encoding="utf-8"?> 
<tle>
<update>#{update}</update>
<data_type>#{@data_type}</data_type>
<sat_name>#{@sat_name}</sat_name>
<first_line>#{@first_line}</first_line>
<second_line>#{@second_line}</second_line>
</tle>
EOS

  end

  def show_test
    now = Time.now.utc
    update = Time.parse(@update.to_s)
    elapsed = now - update
    epy="20"+ @first_line.slice(18,2)
    ep=@first_line.slice(20,12)
    epoch=Time.gm(epy.to_i-1, 12, 31, 0, 0, 0) + ep.to_f*86400

print "Content-type:text/plain\n"
print "Pragma: no-cache\n"
print "Cache-Control: no-cache\n\n"
puts <<EOS
data from   : #{@tle_from}

now         : #{now}
update      : #{@update}
epoch       : #{epoch.to_s}
elapsed     : #{elapsed}
target      : #{@target}
sat_name    : #{@sat_name}
first_line  : #{@first_line}
second_line : #{@second_line}

EOS
if @tle_from != "custom" then
tle = get_remote
puts <<EOS

remote
first_line  : #{@first_line}
second_line : #{@second_line}

EOS
end
end

def update_tle
tle = get_remote
save

print "Content-type:text/plain\n"
print "Pragma: no-cache\n"
print "Cache-Control: no-cache\n\n"
puts <<EOS

remote
first_line  : #{@first_line}
second_line : #{@second_line}

EOS


end
  def show_text
    now = Time.now.utc
    update = Time.parse(@update.to_s)
    elapsed = now - update
    epy="20"+ @first_line.slice(18,2)
    ep=@first_line.slice(20,12)
    epoch=Time.gm(epy.to_i-1, 12, 31, 0, 0, 0) + ep.to_f*86400
  
print "Content-type:text/plain\n"
print "Pragma: no-cache\n"
print "Cache-Control: no-cache\n\n"
puts <<EOS
#{@sat_name}
#{@first_line}
#{@second_line}

===
now         : #{now}
epoch       : #{epoch.to_s}

update      : #{update}
elapsed     : #{elapsed}

EOS
  end

end

class GroundTrack
   attr_reader :gt_array1, :gt_array2, :gt_array3, :gt_array4, :gt_array5, :gt_array6 
  attr_accessor :first_line, :second_line
  def initialize(target="iss")
    @target = target
    @now = Time.now.utc
    @gt_array=[]
    @first_line = "";
    @second_line = "";
  end
  
  def generate
    load_tle
    for i in 0..275
      tmp_time = @now + (i*60)
      position = calc_position(tmp_time)
        @gt_array.push(position)
    end
      show
  end
   
def show
print "Content-type:text/plain\n"
print "Pragma: no-cache\n"
print "Cache-Control: no-cache\n\n"
puts "{'data':["
@gt_array.each{ |point| puts "\{\"lat\":#{point[0]},\"lng\":#{point[1]}\},"}
puts "]}"
end

  def load_tle
    if @target == "custom" then
      tle= LoadTLE.new(@target)
      @orbital_elements = tle.decode(@first_line.to_s,@second_line.to_s)
      
    else
      tle= LoadTLE.new(@target)
      tle.load
      @orbital_elements = tle.decode(tle.first_line.to_s,tle.second_line.to_s)
    end
  end
  
  def calc_position(time)
    orbit = SGP4.new(@orbital_elements);
    position = orbit.calc(time)
    return position
  end
end

#set defalult target & mode
begin
input = CGI.new
target = input["target"].to_s
mode = input["mode"].to_s
if target =="" then target = "iss" end
if mode =="" then mode = "tle" end
sat_name =  CGI.unescape(CGI.escapeHTML(input["sat_name"]))
first_line = CGI.unescape(CGI.escapeHTML(input["1st"]))
second_line = CGI.unescape(CGI.escapeHTML(input["2nd"]))


if mode =="gt" then
  gt = GroundTrack.new(target)
  if target == "custom" then
    gt.first_line = first_line
    gt.second_line = second_line
  end
  gt.generate
elsif mode =="test" then
  tle= LoadTLE.new(target)
  if target == "custom" then
    tle.sat_name = sat_name
    tle.first_line = first_line
    tle.second_line = second_line
  end
  tle.load
  tle.show_test
elsif mode =="text" then
  tle= LoadTLE.new(target)
  if target == "custom" then
    tle.sat_name = sat_name
    tle.first_line = first_line
    tle.second_line = second_line
  end
  tle.load
  tle.show_text

elsif mode =="update" then
  tle= LoadTLE.new(target)
  tle.load
  tle.update_tle

else
  tle= LoadTLE.new(target)
  if target == "custom" then
    tle.sat_name = sat_name
    tle.first_line = first_line
    tle.second_line = second_line
  end
  tle.load
  tle.show
end
rescue  => e
      puts 'Content-Type: text/plain'
      puts
      puts e
      puts e.backtrace
  end #end rescue
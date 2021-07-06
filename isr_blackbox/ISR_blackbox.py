"""ISR_blackbox.py is a module that calculates the performance of an ISR given input parameters.

Contains three classes:
    ISR - fundamental ISR calculator
    Virtual_ISR - allows you to simulate any possible ISR anywhere
    Existing_ISR - allows you to simulate an existing ISR
    
Written by Bill Rideout (brideout@haystack.mit.edu) and Phil Erickson (perickson@haystack.mit.edu)

$Id: ISR_blackbox.py 119 2021-05-12 15:38:08Z brideout $
"""
# standard python imports
import sys
import math
import urllib.request, urllib.error, urllib.parse
import datetime
import copy

# third party imports
import numpy
import numpy.ma
import matplotlib.pyplot
import matplotlib.cm

def run_iri(month,day,hour,lat,lon):
    """
    run_iri returns a numpy recarray with columns (height, ne, ti, te)
    
    Inputs:
        month,dy, hour (hour can be float) - time as UT
        lat, lon - of site
        
    Note uses year 2015 (which as of now is in the future) to turn storm model off.
        
    Runs IRI via web site.
    """
    # make sure lon between 0 and 360
    if lon < 0.0:
        lon += 360.0

    url =  'https://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi'
    req_template = '?model=iri&year=2015&month=%02i&day=%02i&time_flag=0&hour=%f&geo_flag=0.&latitude=%f&longitude=%f&height=100.&profile=1&start=100.&stop=2000.&step=50.&format=0&vars=06&vars=11&vars=14&vars=15'
    req = url + req_template % (month,day,hour,lat,lon)
    response =  urllib.request.urlopen(req)
    text =  response.read().decode('utf-8')

    # parse text
    # first pass is to count the lines
    lineCount = 0
    inLines = False # state variable
    for line in text.split('\n'):
        items = line.split()
        if len(items) != 4:
            continue
        if items[0] == '1' and items[1] == '2':
            inLines = True
            continue
        if inLines:
            try:
                float(items[0])
                lineCount += 1
            except:
                inLines = False
                
    retArr = numpy.zeros((lineCount,), dtype=[('height', float), ('ne', float),
                                              ('ti', float), ('te', float)])
    
    # next pass is to fill out retArr
    lineCount = 0
    inLines = False # state variable
    for line in text.split('\n'):
        items = line.split()
        if len(items) != 4:
            continue
        if items[0] == '1' and items[1] == '2':
            inLines = True
            continue
        if inLines:
            try:
                retArr[lineCount]['height'] = float(items[0])
                retArr[lineCount]['ne'] = float(items[1])
                retArr[lineCount]['ti'] = float(items[2])
                retArr[lineCount]['te'] = float(items[3])
                lineCount += 1
            except:
                inLines = False
                
    return(retArr)


def get_alt_from_range_el(this_range, el):
    """get_alt_from_range_el returns an altitude in km for a given this_range in km and
    el in degrees for a round earth of radius 6,371 km
    """
    re = 6371.0
    
    # h is distance from center of earth to point at end of range
    # c is angle at radar in degrees
    c = 90.0 + el
    # convert to radians
    c = math.radians(c)
    
    h2 = this_range*this_range + re*re - 2*this_range*re*math.cos(c)
    h = math.sqrt(h2)
    
    # alt is h - re
    return(h-re)


def get_range_from_alt_el(alt, el):
    """get_range_from_alt_el returns a range in km for a given alt in km and
    el in degrees for a round earth of radius 6,371 km
    """
    # zoom in by smaller and smaller steps
    this_range = 0.0
    for step in (1000.0,100.0,10.0,1.0):
        while (True):
            this_alt = get_alt_from_range_el(this_range, el)
            if this_alt > alt:
                this_range -= step
                break
            else:
                this_range += step
                
    return(this_range)
        
    



class ISR(object):
    """ISR performs calculations for a particular ISR with particular characteristics
    """
    
    def __init__(self, freq, aperture, t_eff=None):
        """__init__ creates a new ISR instance
        
        Inputs:
            freq - ISR frequency in Hz.
            aperture - Radar aperture in m^2. Antenna efficiency of 0.5 assumed, so give full area
            t_eff - effective temperature of the receiver.  Default is None, in which case, uses sky noise model.
        """
        self._ant_eff = 0.5
        self.freq = float(freq)
        self.aperture = float(aperture) * self._ant_eff
        sky_noise = self.get_sky_noise()
        if t_eff:
            if float(t_eff) < sky_noise and self.freq > 100.0E6:
                raise ValueError('Input system temp %f less than estimated sky noise temp %f' % (float(t_eff), 
                                                                                                  sky_noise))
            elif float(t_eff)*3.0 < sky_noise and self.freq < 100.0E6:
                raise ValueError('Input system temp %f too far below sky noise temp %f' % (float(t_eff), 
                                                                                                  sky_noise))
            self.t_eff = float(t_eff)
        else:
            self.t_eff = self.get_sky_noise()
        self.bw = 50E3 * (freq/440.0E6)
        self._c = 3.0E8 # speed of light in m/s
        self._ka = 0.7 # assume aperture efficiency = 0.7
        self._sigma_e = 1.0E-28 # electron cross section
        self._kb = 1.3806488E-23 # boltzmann constant in m2 kg s-2 K-1
        
        
    def get_signal(self, txp, radar_range, ne, pl):
        """get_signal returns the total signal power in Watts
        
        Inputs:
            txp - Transmit power in Watts
            radar_range - range in km
            ne - electron density in m^-3 at that range
            pl - pulse length in seconds
        """
        cross_section = self.get_cross_section(radar_range, ne, pl)
        return(txp * self.get_radar_loss(cross_section, radar_range))
        
        
        
    def get_radar_loss(self, cross_section, radar_range):
        """get_radar_loss returns the Power received/Power transmitted ratio
        
        Inputs:
            cross_section - target cross section in m^2
            radar_range - range in km
        """
        gain = self.get_gain()
        num = gain * cross_section * self.aperture * self._ka
        dem = math.pow((4*numpy.pi), 2) * math.pow(radar_range*1000.0, 4)
        return(num/dem)
    
    
    def get_cross_section(self, radar_range, ne, pl):
        """get_cross_section calculates the effective cross section as a function of
        range, electron density, pulse length, and gain.  Because this is a soft target,
        cross_section goes up with range squared assuming constant electron density.
        
        Inputs:
            radar_range - range in km
            ne - electron density in m^-3 at that range
            pl - pulse length in seconds
            gain - radar gain
            
        Returns: the effective cross section im m^2
        """
        pdist = pl * self._c
        return(4*numpy.pi*math.pow(radar_range*1000.0,2)*pdist*ne*self._sigma_e/self.get_gain())
        
        
    def get_gain(self):
        """get_gain returns gain as a number (not dB) for this antenna.  Raise ValueError is
        aperture or freq not set
        """
        wavelength = self.get_wavelength()
        return((4*numpy.pi*self.aperture*self._ka)/math.pow(wavelength, 2.0))
        
        
    def get_wavelength(self):
        return(self._c / self.freq)
    
    
    def get_noise(self):
        """get_noise returns the noise power for a given ISR
        """
        return(self._kb * self.t_eff * self.bw)
    
    
    def get_acf_err_ratio(self, count, txp, radar_range, ne, pl, pulse_type):
        """get_acf_err_ratio returns error ratio for ACF with count pulses
        
        Inputs:
            count - the number of pulses combined in the ACF
            txp - Transmit power in Watts
            radar_range - range in km
            ne - electron density in m^-3 at that range
            pl - pulse length in seconds
            pulse_type - a string describing pulse type
                see (_validate_pulse_types)
        """
        self._validate_pulse_types(pulse_type)
        signal = self.get_signal(txp, radar_range, ne, pl)
        if pulse_type == 'mracf':
            # mracf effectively reduces signal by 7 and increases count by 7 by freq multiplexing
            count *= 7
            signal /= 7.0
        noise = self.get_noise()
        if pulse_type == 'alternating_code':
            # very crude model of alternating code added noise - multiply existing noise times 2.0
            noise *= 2.0
        term1 = (signal+noise)/signal
        term2 = 1.0/math.sqrt(count)
        return(term1*term2)
    
    
    def get_range_resolution(self, pl, chip_len=None):
        """get_range_resolution returns a range resolution for a given pl and chip_len
        
        pl - pulse length in seconds
        chip_len - chip length in seconds (chip is modulation on pulse). Default is None,
            meaning no modulation in pulse (single pulse)
            
        Returns range resolution in km
        """
        c = 3.0E5 # speed of light in km/s
        this_len = chip_len
        if this_len == None:
            this_len = float(pl)
        return((this_len*c)/2.0)
    
    
    def get_sky_noise(self):
        """get_sky_noise returns noise temperature based on sky noise model I got from
        Modern Antenna Design by Thomas A. Milligan
        """
        x = math.log10(self.freq/1.0E6)
        y = -2.11*x + 7.7
        return(math.pow(10.0, y))
    
    
    def _validate_pulse_types(self, pulse_type):
        """_validate_pulse_types raises an error if passed in an unknown pulse_type
        """
        valid_types = ['single_pulse',
                       'alternating_code',
                       'mracf']
        if pulse_type not in valid_types:
            raise ValueError('pulse_type must be one of <%s>, not <%s>' % (str(valid_types),
                                                                            pulse_type))
    
    
    
class Virtual_ISR(ISR):
    """Virtual_ISR allows you to test the performance of an ISR with any characteristics you set
    """
    
    def __init__(self, freq, aperture, lat, lon, t_eff=None, txp=1.0E6):
        """__init__ creates a new Virtual_ISR instance
        
        Inputs:
            freq - ISR frequency in Hz.
            aperture - Radar aperture in m^2.
            lat - ISR latitude. Used in IRI modeling of ionosphere.
            lon - ISR longitude (-180 to 180). Used in IRI modeling of ionosphere.
            t_eff - effective temperature of the receiver.  Default is None, in which case, uses sky noise model.
            txp - Transmit power in Watts - Default is 1 MW
        """
        ISR.__init__(self, freq, aperture, t_eff)
        if lat >= -90.0 and lat <= 90.0:
            self.lat = lat
        else:
            raise ValueError('lat must be between -90 and 90, not %f' % (float(lat)))
        if lon > 180.0:
            lon -= 360.0
        if lon >= -180.0 and lon <= 180:
            self.lon = lon
        else:
            raise ValueError('lon must be between -180 and 180, not %f' % (float(lon)))
        self.txp = float(txp)
        self.min_el = None
        
        
    def set_min_el(self, min_el):
        """sets min_el
        """
        if min_el > 90 or min_el < 0:
            raise ValueError('min_el must be between 0 and 90, not %s' % (str(min_el)))
        self.min_el = min_el
        
        
    def set_t_eff(self, new_t_eff):
        """set_t_eff modifies self.t_eff
        """
        self.t_eff = float(new_t_eff)
        
        
        
    def get_acf_err_ratio(self, month, day, hour, el, count, pl, pulse_type):
        """get_acf_err_ratio gets the acf_err_ratio at this virtual radar for the
        given time and configuration.
        
        Inputs:
            month,dy, hour (hour can be float) - time as UT
            el - beam elevation (90 vertical)
            count - the number of pulses combined in the ACF
            pl - pulse length in seconds
            pulse_type - a string describing pulse type
                see (ISR._validate_pulse_types)
                
        Returns:
             numpy.recarray with columns (height (km), range (km), err_ratio, ne)
        """
        if self.min_el != None:
            if el < self.min_el:
                raise ValueError('el %f below minimum elevation %f' % (el, self.min_el))
            
        # first run iri to determine ne versus alt
        iri = run_iri(month, day, hour, self.lat, self.lon)
        
        retArr = numpy.zeros((len(iri['height']),), dtype=[('height', float), ('range', float),
                                              ('error_ratio', float), ('ne', float)])
        
        for i, height in enumerate(iri['height']):
            this_range = get_range_from_alt_el(height, el)
            this_err_ratio = ISR.get_acf_err_ratio(self, count, self.txp, this_range, 
                                                   iri['ne'][i], pl, pulse_type)
            # set all values
            retArr['height'][i] = height
            retArr['range'][i] = this_range
            retArr['error_ratio'][i] = this_err_ratio
            retArr['ne'][i] = iri['ne'][i]
            
        return(retArr)
    
    
    def plot_acf_err_ratio(self, plot_file, start_month, start_day, start_hour, total_hours, el, count, pl, pulse_type,
                           text=None, ipp=None):
        """get_acf_err_ratio gets the acf_err_ratio at this virtual radar for the
        given time and configuration.
        
        Inputs:
            plot_file - path to plot file to save.  If None, do not save
            start_month,start_day, start_hour - start time as UT
            total_hours - total hours to plot (int)
            el - beam elevation (90 vertical)
            count - the number of pulses combined in the ACF
            pl - pulse length in seconds
            pulse_type - a string describing pulse type
                see (ISR._validate_pulse_types)
            text - text to add to figure, if not None (the default)
            ipp - interpulse period in seconds.  If not None, will mask any altitude with interference
                from adjacent pulses.
                
        Returns:
             None
        """
        total_hours = int(total_hours)
        if total_hours < 3:
            raise ValueError('total_hours must be at least three')
        
        if ipp:
            max_range = (3.0E5 * ipp)/2.0
            max_alt = get_alt_from_range_el(max_range, el)
        else:
            max_alt = None
        
        # create a false time
        dt = datetime.datetime(2015, start_month, start_day, int(start_hour))
        for hour in range(total_hours):
            # create a false time
            this_dt = dt + datetime.timedelta(hours=hour)
            result = self.get_acf_err_ratio(this_dt.month, this_dt.day, this_dt.hour, 
                                            el, count, pl, pulse_type)
            if hour == 0:
                y = result['height']
                y = numpy.reshape(y, (len(y),1))
                x = numpy.ones((len(y),1)) * hour
                c = result['error_ratio'] * 100.0
                c = numpy.reshape(c, (len(c), 1))
                ne = numpy.log10(result['ne'])
                ne = numpy.reshape(ne, (len(ne), 1))
            else:
                new_y = result['height']
                new_y = numpy.reshape(new_y, (len(new_y),1))
                y = numpy.concatenate((y, new_y), 1)
                new_x = numpy.ones((len(new_y),1)) * hour
                x = numpy.concatenate((x, new_x), 1)
                new_c = result['error_ratio'] * 100.0
                new_c = numpy.reshape(new_c, (len(new_c), 1))
                c = numpy.concatenate((c, new_c), 1)
                new_ne = numpy.log10(result['ne'])
                new_ne = numpy.reshape(new_ne, (len(new_ne), 1))
                ne = numpy.concatenate((ne, new_ne), 1)
        
        thisCmap = copy.copy(matplotlib.cm.jet)
        if max_alt:
            if max_alt < 2000.0:
                nan_indexes = numpy.where(new_y > max_alt)[0]
                c[nan_indexes,:] = -1.0
                c = numpy.ma.masked_less(c, 0.0)
                thisCmap = copy.copy(thisCmap.set_bad('white', alpha=None))
                
        titleStr = 'Error percentage for ISR ACF'
        xlabel = 'Hours since UT month=%i, day=%i, hour=%i' % (start_month, start_day, int(start_hour))
        matplotlib.pyplot.clf()
        fig, (ax1, ax2) = matplotlib.pyplot.subplots(1, 2)
        fig.set_size_inches(12.8, 9.6)
        pmap = ax1.pcolormesh(x,y,c, vmin=0.0, vmax=25.0, cmap=thisCmap, shading='auto')
        fig.colorbar(pmap, ax=ax1)
        ax1.set_title(titleStr)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel('Altitude in km')
        if text:
            matplotlib.pyplot.figtext(0.2, 0.2, text, color='k')
        if max_alt:
            if max_alt < 2000.0:
                matplotlib.pyplot.figtext(0.2, 0.85, 'Alt limited by IPP')
                
        pmap2 = ax2.pcolormesh(x,y,ne, vmin=9.0, vmax=13.0, shading='auto')
        fig.colorbar(pmap2, ax=ax2)
        ax2.set_title('IRI electron density (log10(m^-3))')
        ax2.set_xlabel(xlabel)
        
        if not plot_file is None:
            matplotlib.pyplot.savefig(plot_file)
        else:
            matplotlib.pyplot.show()
                
        
class ISR_Mode:
    """ISR_mode is a class that describes the particular mode an ISR can run in
    """
    def __init__(self, pl, ipp, chip_len, pulse_type, t_eff=None):
        """
        All inputs are set as attributes.
        
        Inputs:
            pl - pulse length in seconds
            ipp - interpulse period in seconds
            chip_len - sub pulse length in seconds - may be None if no sub pulse
            pulse_type - a string describing pulse type
                see (_validate_pulse_types)
            t_eff - effective system temperature of this mode.  Used if different from
                default radar system temperature.  If None (the default), use default
                radar system temperature.
        """
        self.pl = float(pl)
        self.ipp = float(ipp)
        if ipp <= pl:
            raise ValueError('ipp %f must be greater than pl %f' % (self.ipp, self.pl))
        self.chip_len = chip_len
        self.pulse_type = pulse_type
        self.t_eff = t_eff
            
        

class Existing_ISR_Radars:
    """Existing_ISR_Radars allows you to predict the performance of all existing ISR radars
    """
    # freq, area, lat, lon, nf, min el, txp
    _isr_characteristics = { \
        'millstone_zenith': (440.0E6, 3.14*math.pow(68.0/2,2), 42.61950, -71.5, 150, 88.0, 1.5E6),
        'millstone_misa': (440.0E6, 3.14*math.pow(46.0/2,2), 42.61950, -71.5, 150, 4.0, 1.5E6),
        'poker_flat': (450.0E6, 900, 65.13, -147.47, None, 35, 1.7E6),
        'resolute_bay': (442.0E6, 900, 74.73, -94.9, None, 16, 1.7E6),
        'arecibo': (430.0E6, 3.14*math.pow(305.0/2,2), 18.34, -66.75, 150, 70.0, 1.5E6),
        'jicamarca':(49.92E6, 85000, -11.95, -143.62, 5000.0, 85, 2E6),
        'svalbard_fixed': (500.0E6, 3.14*math.pow(42.0/2,2), 78.09, 16.02, None, 81, 1.0E6),
        'svalbard_steerable': (500.0E6, 3.14*math.pow(32.0/2,2), 78.09, 16.02, None, 30, 1.0E6),
        'eiscat_tromso_vhf': (223.0E6, 120.0*40.0, 69.583, 19.21, None, 30, 1.6E6),
        'eiscat_tromso_uhf': (929.0E6, 3.14*math.pow(32.0/2,2), 69.583, 19.21, None, 22, 1.6E6),
        'sondrestrom': (1290.0E6, 3.14*math.pow(32.0/2,2), 67.0, -51.0, None, 22, 3.2E6)
        }
    
    _isr_modes = {('millstone_zenith', 'sp_480'): ISR_Mode(480.0E-6, 8910.0E-6, None, 'single_pulse'),
                  ('millstone_zenith', 'sp_640'): ISR_Mode(640.0E-6, 11600.0E-6, None, 'single_pulse'),
                  ('millstone_zenith', 'sp_960'): ISR_Mode(960.0E-6, 17000.0E-6, None, 'single_pulse'),
                  ('millstone_zenith', 'sp_1280'): ISR_Mode(1280.0E-6, 22400.0E-6, None, 'single_pulse'),
                  ('millstone_zenith', 'sp_2000'): ISR_Mode(2000.0E-6, 34600.0E-6, None, 'single_pulse'),
                  ('millstone_zenith', 'ac_480'): ISR_Mode(480.0E-6, 8910.0E-6, 30.0E-6, 'alternating_code'),
                  ('millstone_misa', 'sp_480'): ISR_Mode(480.0E-6, 8910.0E-6, None, 'single_pulse'),
                  ('millstone_misa', 'sp_640'): ISR_Mode(640.0E-6, 11600.0E-6, None, 'single_pulse'),
                  ('millstone_misa', 'sp_960'): ISR_Mode(960.0E-6, 17000.0E-6, None, 'single_pulse'),
                  ('millstone_misa', 'sp_1280'): ISR_Mode(1280.0E-6, 22400.0E-6, None, 'single_pulse'),
                  ('millstone_misa', 'sp_2000'): ISR_Mode(2000.0E-6, 34600.0E-6, None, 'single_pulse'),
                  ('millstone_misa', 'ac_480'): ISR_Mode(480.0E-6, 8910.0E-6, 30.0E-6, 'alternating_code'),
                  ('arecibo', 'mracf'): ISR_Mode(196.0E-6, 10000.0E-6, None, 'mracf'),
                  ('arecibo', 'topside'): ISR_Mode(500.0E-6, 20000.0E-6, None, 'single_pulse'),
                  ('sondrestrom', 'sp_320'): ISR_Mode(320.0E-6, 8910.0E-6, None, 'single_pulse'),
                  ('sondrestrom', 'ac_320'): ISR_Mode(320.0E-6, 8910.0E-6, 20.0E-6, 'alternating_code'),
                  ('poker_flat', 'sp_480'): ISR_Mode(480.0E-6, 9600.0E-6, None, 'single_pulse'),
                  ('poker_flat', 'ac_480'): ISR_Mode(480.0E-6, 9600.0E-6, 30.0E-6, 'alternating_code'),
                  ('eiscat_tromso_vhf', 'beata'): ISR_Mode(640.0E-6, 5580.0E-6, 20.0E-6, 'alternating_code'),
                  ('eiscat_tromso_vhf', 'bella'): ISR_Mode(1350.0E-6, 11250.0E-6, 45.0E-6, 'alternating_code'),
                  ('eiscat_tromso_vhf', 'manda'): ISR_Mode(146.4E-6, 1500.0E-6, 2.4E-6, 'alternating_code', 350.0),
                  ('eiscat_tromso_vhf', 'tau7'): ISR_Mode(1536.0E-6, 15624.0E-6, 96.0E-6, 'alternating_code'),
                  ('eiscat_tromso_uhf', 'beata'): ISR_Mode(640.0E-6, 5580.0E-6, 20.0E-6, 'alternating_code'),
                  ('eiscat_tromso_uhf', 'bella'): ISR_Mode(1350.0E-6, 11250.0E-6, 45.0E-6, 'alternating_code'),
                  ('eiscat_tromso_uhf', 'manda'): ISR_Mode(146.4E-6, 1500.0E-6, 2.4E-6, 'alternating_code', 150.0),
                  ('svalbard_fixed', 'ipy'): ISR_Mode(900.0E-6, 3750.0E-6, 30E-6, 'alternating_code'),
                  ('svalbard_fixed', 'beata'): ISR_Mode(1500.0E-6, 6250.0E-6, 50.0E-6, 'alternating_code'),
                  ('svalbard_fixed', 'manda'): ISR_Mode(256E-6, 1250.0E-6, 4E-6, 'alternating_code', 150.0),
                  ('svalbard_fixed', 'tau7'): ISR_Mode(1920.0E-6, 20000.0E-6, 120.0E-6, 'alternating_code'),
                  ('svalbard_steerable', 'ipy'): ISR_Mode(900.0E-6, 3750.0E-6, 30E-6, 'alternating_code'),
                  ('svalbard_steerable', 'beata'): ISR_Mode(1500.0E-6, 6250.0E-6, 50.0E-6, 'alternating_code'),
                  ('svalbard_steerable', 'manda'): ISR_Mode(256E-6, 1250.0E-6, 4E-6, 'alternating_code', 150.0),
                  ('svalbard_steerable', 'tau7'): ISR_Mode(1920.0E-6, 20000.0E-6, 120.0E-6, 'alternating_code'),
                  ('svalbard_steerable', 'folke'): ISR_Mode(960.0E-6, 5000.0E-6, 60.0E-6, 'alternating_code'),
                  ('jicamarca', 'long pulse'): ISR_Mode(1600.0E-6, 41667.0E-6, None, 'single_pulse')
                  }
    
    def __init__(self):
        """__init__ creates a new Existing_ISR_Radars instance
        
        """
        self._radar_dict = {} # a dict of key = radar name, value of Virtual_ISR
        self._radar_modes = {} # a dict of key = radar name, value = dict with key=mode name, value = ISR_Mode
        for name in list(self._isr_characteristics.keys()):
            freq, aperature, lat, lon, t_eff, min_el, txp = self._isr_characteristics[name]
            self._radar_dict[name] = Virtual_ISR(freq, aperature, lat, lon, t_eff, txp=txp)
            self._radar_dict[name].set_min_el(min_el)
            self._radar_modes[name] = {}
            for key in list(self._isr_modes.keys()):
                isr_name, mode_name = key
                if isr_name == name:
                    self._radar_modes[name][mode_name] = self._isr_modes[key]
                    
                    
    def get_radar_names(self):
        """get_radar_names returns a list of defined radar names
        """
        return(list(self._radar_dict.keys()))
    
    
    def get_mode_names(self, radar_name):
        """get_mode_names returns a list of mode names for a given radar
        """
        if radar_name not in self.get_radar_names():
            raise ValueError('Radar %s unknown' % (str(radar_name)))
        return(list(self._radar_modes[radar_name].keys()))
                    
                    
    def get_signal(self, radar_name, radar_range, ne, mode_name):
        """get_signal returns the total signal power in Watts
        
        Inputs:
            radar_name - radar name
            radar_range - range in km
            ne - electron density in m^-3 at that range
            mode_name
        """
        if radar_name not in self.get_radar_names():
            raise ValueError('Radar %s unknown' % (str(radar_name)))
        if mode_name not in self.get_mode_names(radar_name):
            raise ValueError('For radar %s, mode %s unknown' % (str(mode_name), str(radar_name)))
        vir_radar = self._radar_dict[radar_name]
        txp = vir_radar.txp
        pl = self._radar_modes[radar_name][mode_name].pl
        return(vir_radar.get_signal(txp, radar_range, ne, pl))
    
    
    def get_range_resolution(self, radar_name, mode_name):
        """get_range_resolution returns a range resolution
        
        Inputs:
            radar_name - radar name
            mode_name
            
        Returns range resolution in km
        """
        if radar_name not in self.get_radar_names():
            raise ValueError('Radar %s unknown' % (str(radar_name)))
        if mode_name not in self.get_mode_names(radar_name):
            raise ValueError('For radar %s, mode %s unknown' % (str(mode_name), str(radar_name)))
        vir_radar = self._radar_dict[radar_name]
        pl = self._radar_modes[radar_name][mode_name].pl
        chip_len = self._radar_modes[radar_name][mode_name].chip_len
        return(vir_radar.get_range_resolution(pl, chip_len))
    
    
    def get_acf_err_ratio(self, radar_name, mode_name, month, day, hour, el, integration_seconds):
        """get_acf_err_ratio gets the acf_err_ratio at this virtual radar for the
        given time and configuration.
        
        Inputs:
            radar_name - radar name
            mode_name
            month,dy, hour (hour can be float) - time as UT
            el - beam elevation (90 vertical)
            integration_seconds - the number of seconds to integrate
                
        Returns:
             numpy.recarray with columns (height (km), range (km), err_ratio)
        """
        if radar_name not in self.get_radar_names():
            raise ValueError('Radar %s unknown' % (str(radar_name)))
        if mode_name not in self.get_mode_names(radar_name):
            raise ValueError('For radar %s, mode %s unknown' % (str(mode_name), str(radar_name)))
        vir_radar = self._radar_dict[radar_name]
        pl = self._radar_modes[radar_name][mode_name].pl
        ipp = self._radar_modes[radar_name][mode_name].ipp
        pulse_type = self._radar_modes[radar_name][mode_name].pulse_type
        count = int(integration_seconds/ipp)
        return(vir_radar.get_acf_err_ratio(month, day, hour, el, count, pl, pulse_type))
    
    
    def plot_acf_err_ratio(self, radar_name, mode_name, plot_file, start_month, start_day, start_hour, total_hours, el, integration_seconds):
        """get_acf_err_ratio gets the acf_err_ratio at this virtual radar for the
        given time and configuration.
        
        Inputs:
            radar_name - radar name
            mode_name
            plot_file - path to plot file to save.  If None, do not save
            start_month,start_day, start_hour - start time as UT
            total_hours - total hours to plot (int)
            el - beam elevation (90 vertical)
            integration_seconds - the number of seconds to integrate
                
        Returns:
             None
        """
        if radar_name not in self.get_radar_names():
            raise ValueError('Radar %s unknown' % (str(radar_name)))
        if mode_name not in self.get_mode_names(radar_name):
            raise ValueError('For radar %s, mode %s unknown' % (str(mode_name), str(radar_name)))
        vir_radar = self._radar_dict[radar_name]
        pl = self._radar_modes[radar_name][mode_name].pl
        ipp = self._radar_modes[radar_name][mode_name].ipp
        duty_cycle = pl/ipp
        pulse_type = self._radar_modes[radar_name][mode_name].pulse_type
        count = int(integration_seconds/ipp)
        if self._radar_modes[radar_name][mode_name].t_eff != None:
            vir_radar.set_t_eff(self._radar_modes[radar_name][mode_name].t_eff)
        text = 'Radar: %s\nMode: %s\npl: %g, ipp: %g\nStart month=%i, day=%i, hour=%i\nInt secs=%i, el=%3.1f\nduty cycle: %5.3f, num pulses: %i' % (radar_name, 
                                                               mode_name, pl, ipp, start_month,start_day, start_hour, 
                                                               integration_seconds, el, duty_cycle, count)
        return(vir_radar.plot_acf_err_ratio(plot_file, start_month, start_day, start_hour, total_hours, 
                                            el, count, pl, pulse_type, text, ipp))
        
        
class ISR_Blackbox_Service:
    """this class support the ISR_blackbox demo applications
    """

    def getRadarNames(self):
        """getRadarNames returns a list of radar names with modes defined
        """
        isrs = Existing_ISR_Radars()
        radar_names = isrs.get_radar_names()
        ret_list = [] # only accept names with modes defined
        for radar_name in radar_names:
            if len(isrs.get_mode_names(radar_name)) > 0:
                ret_list.append(radar_name)
                
        return(ret_list)
    
    
    def getModeNames(self, radar_name):
        """getModeNames returns a list of modes names for a given radar
        """
        isrs = Existing_ISR_Radars()
        return(isrs.get_mode_names(radar_name))
    
    
    def plotAcfErrRatio(self, radar_name, mode_name, start_month, start_day, start_hour, total_hours, el, integration_seconds):
        """plotAcfErrRatio gets the acf_err_ratio at this existing radar for the
        given time and configuration.
        
        Inputs:
            radar_name - radar name
            mode_name
            start_month,start_day, start_hour - start time as UT
            total_hours - total hours to plot (int)
            el - beam elevation (90 vertical)
            integration_seconds - the number of seconds to integrate
                
        Returns:
             None
        """
        # force date to be rational
        try:
            datetime.datetime(2015, start_month, start_day)
        except:
            if start_month == 2:
                start_day = 29
            else:
                start_day = 30
                
        
            
        imagePlot = None
        
        isrs = Existing_ISR_Radars()
        isrs.plot_acf_err_ratio(radar_name, mode_name, imagePlot, start_month, start_day, start_hour, 
                                total_hours, el, integration_seconds)
            
        
        
        
    def plotAcfErrRatioVirtualISR(self, freq, aperture, lat, lon, txp, pl,
                                  start_month, start_day, start_hour, total_hours, el, 
                                  integration_seconds, ipp):
        """plotAcfErrRatio gets the acf_err_ratio at this virtual radar for the
        given time and configuration.
        
        Inputs:
            freq - ISR frequency in Hz.
            aperture - Radar aperture in m^2.
            lat - ISR latitude. Used in IRI modeling of ionosphere.
            lon - ISR longitude (-180 to 180). Used in IRI modeling of ionosphere.
            txp - Transmit power in Watts
            pl - single pulse length in seconds
            start_month,start_day, start_hour - start time as UT
            total_hours - total hours to plot (int)
            el - beam elevation (90 vertical)
            integration_seconds - the number of seconds to integrate
            ipp - ipp length in seconds
                
        Returns:
             None
        """
        # force date to be rational
        try:
            datetime.datetime(2015, start_month, start_day)
        except:
            if start_month == 2:
                start_day = 29
            else:
                start_day = 30
                
            
        imagePlot = None
            
        count = int(integration_seconds/ipp)
        
        text = """
        freq=%g aperture=%g m^2,
        lat=%3.1f, lon=%4.1f, 
        txp=%g, pl=%g,
        start_month=%i, start_day=%i,
        start_hour=%i, pulse count=%i,
        total_hours=%i, el=%i
        int_sec=%i, duty cycle=%5.3f,""" % (freq, aperture, lat, lon, txp, pl, 
                                                start_month, start_day, start_hour, count,
                                                int(total_hours), int(el), int(integration_seconds),
                                                pl/ipp)
        
        v_isr = Virtual_ISR(freq, aperture, lat, lon, txp=txp)
        v_isr.plot_acf_err_ratio(imagePlot, start_month, start_day, start_hour, total_hours, el, count, pl, 
                                 'single_pulse', text, ipp=ipp)
        
        
        return(None)
        

        


    
        
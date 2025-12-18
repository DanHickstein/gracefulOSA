# -*- coding: utf-8 -*-
"""
Software for collecting spectra using a Yokogowa/Ando/Tohorlabs spectrum 
analyzer (OSA) while rotating a waveplate using a Thorlabs K10CR1. So far, 
this has been tested with the following spectrum analyzers:

 - ANDO AQ-6315A
 - Yokogawa AQ6375 (set to AQ6317 compatibility mode)

This should work with any of the modern Yokogawa OSAs that have the AQ6317
legacy compatiblity mode. To activate this mode press 

   System >> More >> GPIB Settings >> Command Format >> AQ6317.

This software assumes GPIB communication. It is necessary to manually set the
GPIB address in the script to match the optical spectrum analyzer.

This program requires that the Thorlabs Kinesis drivers are in the same folder
as this program.
pyVISA is used to communicate with the ANDO OSA using GPIB.

To calibrate the power, use the Kinesis software to move the stage between the
maximum and minimum power levels. Record the results in the first parameter
block below.

For the Ando/Yokogawa OSAs, the sensitivity can be selected using the val_sens
variable. The following options are available, with the time to acquire each
spectrum for the Ando AQ6317 in parentheses:
    
    0 = Norm Range Hold (7   sec) - not recommended, spectra are terrible
    1 = Norm Range Auto (10  sec)
    2 = High 1          (1.2 min)
    3 = High 2          (2.4 min)
    4 = High 3          (~17 min)
    

TODO: 
    Save acquisition parameters in the Thorlabs datafiles.
    Automatically check that power array doesn't exceed limits
    It would be cool if the homing and OSA warmup could occur at the same time.
      So, non-blocking stage homing routine.

"""

from __future__ import print_function
import ctypes as c
import numpy as np
import os, time
import sys
import pyvisa as visa
import datetime
import platform
import serial



# ------------- ENTER PARAMETERS HERE: ----------------------------------------
# set which OSA type will be used:
# OSA_type = 'ando'  # Also works with Yokogawa OSAs in AQ6317 compatibilty mode
OSA_type = 'thorlabs'  # Tested with OSA205

osa_address = "GPIB0::23::INSTR"
# osa_address = "GPIB1::23::INSTR"

# Set how the power will be scanned. 
# K10CR1 rotation stage or EDFA100P amplifier

power_type = 'K10CR1'

# power_type = 'EDFA100P'
edfa_comport = 16


# if the EDFA100P is used for power scan, these parameters are important.
# Measure the power at some low current (but where you still get gain, say 500 mA)
# then measure the power at some higher current, say the full current at 1000 mA.
# the current-to-power function is linearly interpolated/extrapolated from this. 
currentLo = 400
currentHi = 1000

powerLo = 32
powerHi = 93

# is the K10CR1 is used for power scan, these parameters are important:
extinguishingAngle = 44 # Angle at which minimum power occurs. Only for rotation
maxPower = 37.5       # Maximum Power through half wave plate
minPower = 0.1       # Minimum power through half wave plate

# Scan parameters:
numPoints = 50   # number of data points (different power) to be collected
scanPowerLo = 0.2  #0.1 #high numbers for the menhir
scanPowerHi = 37   #73 #high numbers for the menhir
powerarray = np.linspace(scanPowerLo, scanPowerHi, numPoints)  # adjust range of scan
powerarray = powerarray[::-1]

# Ando/Yokogawa parameters:
val_sens = int(1)

# Thorlabs OSA parameters:
avenum = 1

# -------------- End of parameters---------------------------------------------


def main():
    
    if power_type == 'EDFA100P':
        power_scanner = EDFA(edfa_comport, currentLo=currentLo, 
                             currentHi=currentHi, powerLo=powerLo, powerHi=powerHi)
        units = 'mA'
        
    elif power_type == 'K10CR1':
        power_scanner = K10CR1(maxPower=maxPower, minPower=minPower, 
                               extinguishingAngle=extinguishingAngle)
        units = 'deg'
    else:
        raise ValueError('power_type not recognized: %s'%power_type)
    

    # initialize the OSA
    osa = OSA(OSA_type, val_sens, avenum)

    #---Create Base Directory for saving data
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    cwd = os.getcwd()
    data_dir = os.path.join(cwd, 'data')
    if not(os.path.isdir(data_dir)):
        os.mkdir(data_dir)

    base_dir = os.path.join(data_dir, today)
    if not(os.path.isdir(base_dir)):
        os.mkdir(base_dir)

    run_counter = 1
    run_filename  = 'run_%04i.txt'%(run_counter)

    # find the first available file name:
    while os.path.exists(os.path.join(base_dir, run_filename)):
        run_counter = run_counter+1
        run_filename = 'run_%04i.txt'%(run_counter)
    filename = os.path.join(base_dir, run_filename)

    print('Saving to:   %s\n' %(filename))

    with open(filename, 'w') as file:
        time_now  = datetime.datetime.now().strftime("%Y-%m-%d %X")
        
        file.write('OSA type: %s\n'%OSA_type)
        file.write('Power type: %s\n'%power_type)
        file.write('Time\t'+time_now+'\n')
        file.write('Num points: %i\n'%numPoints)
        file.write('Powers: ')
        for p in powerarray:
            file.write('%.2f\t'%p)
        file.write('\nThe first line of data has the wavelengths in nm\n')
        file.write('All subsequent lines consist of\n')
        file.write('FileNum\tPower (mW)\t Current (mA)\n')
        file.write('followed by the dBm.\n')
        file.write('Data:\n')

        wls, levels = osa.get_data()
        
        for wl in wls:
            file.write('%.2f\t'%wl)
        file.write('\n')
        
        deg_mA = power_scanner.set_power(powerarray[0])
        t_start = time.time()

        for count, power in enumerate(powerarray):
            deg_mA = power_scanner.set_power(power)
            print('Run %04i % 3i of % 3i - %5.3f %s - %6.2f mW - '%(run_counter, count+1, powerarray.size, deg_mA, units, powerarray[count]), end='')
            sys.stdout.flush()
    
            wls, levels = osa.get_data()
            
            file.write('%03i\t% .4f\t% .3f\t'%(count, power, deg_mA))
            
            for level in levels:
                file.write('%.2f\t'%level)
            file.write('\n')
            
            elapsed = time.time()-t_start
            total_time = elapsed/(count+1) * powerarray.size
            print('%.1f/%.1f sec.'%(time.time()-t_start, total_time))
            
    power_scanner.close()

   
class OSA:
    def __init__(self, OSA_type, val_sens, ave_num):
        self.OSA_type = OSA_type
        self.val_sens = val_sens
        self.ave_num = ave_num
        if self.OSA_type == 'ando':
            rm = visa.ResourceManager()
            print('ando OSAs detected: ')
            print(rm.list_resources())
            self.osa = rm.open_resource(osa_address)
            sens = ["SNHD", "SNAT", "SHI1", "SHI2", "SHI3"]
            self.osa.write(sens[val_sens])
                
        elif OSA_type == 'thorlabs':
            # # dllnameOSA = os.path.join('C:\\Program Files\\Thorlabs\\Thorlabs OSA\\FTSLib.dll')
            dllnameOSA = os.path.join('C:\\Program Files\\Thorlabs\\ThorSpectra\\FTSLib.dll')
    
            print(dllnameOSA)
            if not os.path.exists(dllnameOSA):
                print("ERROR: DLL not found")
            #o = c.windll.LoadLibrary(dllnameOSA)  #Alternate between dll loading method
            self.o = c.CDLL(dllnameOSA) # Altnate between dll loading method
            
            self.sn = c.c_int(0) # index of the OSA to be used. 
            try:
                self.o.FTS_CloseSpectrometer(self.sn)
                print('Previous OSA connection closed.')
            except:
                pass
            print('Initializing spectrometer(s)...')
            numSpectrometers = self.o.FTS_InitializeSpectrometers()
            print('(%i) Thorlabs OSAs connected'%numSpectrometers)
            
            
            self.o.FTS_OpenSpectrometer(self.sn)
            self.o.FTS_SetAcquisitionOption_AverageSpectrum(self.sn, c.c_int(1))    
            
            # get a few spectra to get the Thorlabs OSA "warmed up"
            for e in range(10):
                print('Warming up the OSA by getting spectrum %i/%i'%(e, 10))
                self.get_data()
            
        else:
            raise ValueError('OSA_type not recognized.')
            
    
    def get_data(self):
        if self.OSA_type == 'ando':
            #Tells OSA to begin sweep
            self.osa.write("SGL")
            query = int(self.osa.query('SWEEP?')) # greater that zero means OSA is currently performing a sweep
        
            #Checking if OSA Done Sweeping
            while query > 0:
                time.sleep(.2) # in seconds
                query = int(self.osa.query('SWEEP?'))
        
            ### Capturing Data Trace from OSA
        
            # Active Trace
            t_active = int(self.osa.query("ACTV?"))
            trace = "ABC"[t_active]
        
            # # Instrument ID
            # osa_ID = ''.join([i if ord(i) < 128 else ' ' for i in osa.read_raw().rstrip()]) # strips non-ASCII characters
            osa_ID = str(self.osa.read_raw().rstrip()) # strips non-ASCII characters
        
            # Time Stamp
            time_now = datetime.datetime.now().strftime("%Y-%m-%d %X")
        
            # Measurement Characteristics
            t_list_hds = "Center Wvl:,Span Range:,REF Level:,Level Scale:,Wvl Resolution:,Avg Count:,Sampl Count:,Sensitivity:,Monochro:,Waveform Type:".split(',')
            t_list_cmds = ['CTRWL', 'SPAN', 'REFL', 'LSCL', 'RESLN', 'AVG', 'SEGP', 'SENS', 'MONO', 'TR' + trace]
            t_list = []
            
            d_mono = {0: 'SGL', 1: 'DBL'}
            d_sens = {1: 'HIGH 1', 2: 'HIGH 2', 3: 'HIGH 3', 4: 'NORMAL RANGE HOLD', 5: 'NORMAL RANGE AUTO'}
            d_trace = {0: {0: 'MEAS', 1: 'FIX', 2: 'MAX HOLD', 3: 'ROL AVG'},
                       1: {0: 'MEAS', 1: 'FIX', 2: 'MIN HOLD', 3: 'ROL AVG'},
                       2: {0: 'MEAS', 1: 'FIX', 2: 'A-B', 3: 'B-A', 4: 'A-B (LIN)', 5: 'B-A (LIN)',6: 'A+B (LIN)', 
                           7: 'NORMALIZE', 8: 'DOMINANT', 10: 'CURVE FIT', 110: 'CURVE FIT PK'}}
            d_commands = {'MONO': d_mono, 'SENS': d_sens, 'TRA': d_trace[0], 'TRB': d_trace[1], 'TRC': d_trace[2]}
            
            for cmd in t_list_cmds:
                cmd_response = self.osa.query(cmd+'?').rstrip()
                try:
                    t_dic = d_commands[cmd]
                    try:
                        t_list.append(t_dic[int(cmd_response)])
                    except KeyError:
                        t_list.append('N/A')
                except KeyError:
                    t_list.append(cmd_response)
                                 
            # Spectral Data
            self.osa.write("LDTDIG3") #sets retrieval the maximum of 3 decimal places
            level_unit = ["W","dBm"][bool(float(self.osa.query("LSCL?")))]
            abs_or_dens = ["","/nm"][int(self.osa.query("LSUNT?"))]
            t_wave = self.osa.query("WDAT"+trace).rstrip().split(',')[1:] #discards the sample count
            t_level = self.osa.query("LDAT"+trace).rstrip().split(',')[1:]
             
            # # Format Data String:
            # col_1 = ["Instrument:"] + ["Time Stamp:"] + t_list_hds + ["", "Wavelength(nm)"] + t_wave
            # col_2 = [osa_ID] + [time_now] + t_list + ["", "Level("+level_unit+abs_or_dens+")"] + t_level
            # col_comb = zip(col_1, col_2)
            # data_list = []
            # for data_row in col_comb:
            #     data_list.append('\t'.join(data_row))
            # data_string = "\n".join(data_list)
            wls = np.array([float(x) for x in t_wave])
            lvl = np.array([float(x) for x in t_level])

            return wls, lvl
            
            
        elif self.OSA_type == 'thorlabs':
            def GrabSpectrum():
                global grabbing, x, y, xcfloat, ycfloat, err1, err2
                grabbing = True
                #t3 = time.time()
                
                def retreiveSpectrum():
                     global grabbing, x, y, xcfloat, ycfloat, err1, err2
                     xcfloat = (c.c_float*16774)()
                     ycfloat = (c.c_float*16774)()
                     p = (c.c_float*16774)()
                     f = c.c_int()
                     l = c.c_uint()
                     
                     err1 = 10
                     err2 = 10
             
                     t2 = time.time()
                     self.o.FTS_GetLastSpectrumArrays(self.sn, c.pointer(ycfloat), c.pointer(xcfloat), c.pointer(p), c.pointer(f), c.pointer(l))  #The spectrum has the x-axis unit of wavenumbers and the y-axis unit of mW. 
                     err1 = time.time()-t2
              
                     y = np.array(list(ycfloat))
                     x = np.array(list(xcfloat))        
                     
                     grabbing = False
                     
                         
                def callback(p1,p2,p3):
                    if p2 == 1:
                        #print('Good news from the OSA: Interferogram now available!')
                        a=0
                    elif p2 == 2:
                        #print('Good news from the OSA: Spectrum now available!')
                        a=0        
                        retreiveSpectrum()
                    elif p2 == 3:
                        #print('Good news from the OSA: One of the average spectra has been collected, but the averaged spectrum is not yet available.')
                        a=0
                    else:
                        #print('Very strange news from the spectrometer. p1,p2,p3:',p1,p2,p3)
                        a=0
                    return a
                
                callback_format = c.CFUNCTYPE(c.c_voidp ,c.c_ushort, c.c_uint, c.c_uint)
                CALL_BACK = callback_format(callback)
    
                self.o.FTS_AcquireSingleSpectrum(self.sn, CALL_BACK)
               
                while grabbing==True:
                    time.sleep(0.005)
                return x,y,xcfloat,ycfloat, err1, err2
            
            
            count1 = 0
            xall = []
            yall = []
            while count1<avenum:
                
                x,y,xcfloat,ycfloat,err1,err2 = GrabSpectrum()
                
               # x,y,xcfloat,ycfloat = GrabSpectrum()  #Gets spectrum from previously defined function
                count1=count1+1
                y = y[:-1] # get rid of the zero mW (absolute power) component
                x = x[:-1] # get rid of the zero cm-1 (wavenumber) component        
                grid_size =  x[1]-x[0]
                
                #Converts data from cm-1 & mW to nm & dBm/nm
                y_mW_per_wavenumber = y/grid_size
                nmx=(1E7)/x  #cm-1 to nm
            
                y_mW_per_nm = y_mW_per_wavenumber /( 1e7/(x**2))
                #ynew=y/((5E-8)*nmx*nmx)  #absolute power to relative power
                dby=10*np.log10(y_mW_per_nm) #log scale of relative power
                xall.append(nmx)
                yall.append(dby)
                xa=x
                ya=y
                
            #averages all data collected from while loop 
            yave = np.mean(yall,axis=0,dtype=np.float64)
            xave = np.mean(xall,axis=0,dtype=np.float64)
            
            return xave, yave
        
            # data_list = []
            # data_list.append('Wavelength(nm)\t  (dBm/nm)')
            # for x,y in zip(xave, yave):
            #     data_list.append('%.4f\t%.4f'%(x,y))    
            
            # data_string = "\n".join(data_list)
    
        else:
            raise ValueError('OSA_type not recognized!')
    


################### Commands for the Thorlabs EDFA100P
class EDFA:
    """"Class for interfacing with the Thorlabs EDFA100P and similar. 
    Parameters:
        com_port is the COM port for communications. Use device manager and
         unplug/plug the EDFA to determine which COM port the EDFA is assigned. 
        int baudrate is for serial communications, default is 115200
        int bytesize is for serial communications, defauly is 8
        float delay is the delay between issuing commands. It takes the EDFA a little
         while to respond, so use this to let it settle to the new power before 
         doing something else. 0.5 is the default."""
         
    def __init__(self, com_port, currentLo, currentHi, powerLo, powerHi,
                 baudrate=115200, bytesize=8, delay=0.5):
        
        try:
            # connect to EDFA
            self.ser = serial.Serial(
                port=f"COM{edfa_comport}", baudrate=baudrate, bytesize=bytesize, timeout=2)
            # check if EDFA is connected
            if self.ser.is_open:
                print('Connection to EDFA already open.')
            else:
                self.ser.open()
                print('Connected to EDFA.')
        except:
            raise NameError(
                '''Could not connect to EDFA. Just try again. 
                If that doesn't fix is, check if the EDFA is on 
                and connected via USB. Also try power cycling the EDFA and 
                starting a new python consol.''')
        
        delta_current = currentHi - currentLo
        
        if delta_current < 0:
            raise ValueError('currentLo must be smaller than currentHi')
            
        delta_power = powerHi - powerLo
        
        if delta_power < 0:
            raise ValueError('powerLo must be smaller than powerHi')
            
        self.slope = delta_power / delta_current
        self.intercept = powerLo - self.slope*currentLo
        
    def set_power(self, power):
        # calculate current!
        current = (power - self.intercept)/self.slope
        command = f'current=%.3f\r'%current
        self.ser.write(command.encode('utf-8'))
        return current

    def close(self):
        self.ser.close()
        
        
################### Commands for the Rotation Stage:
class K10CR1:
    def __init__(self, maxPower, minPower, extinguishingAngle):
        self.maxPower = maxPower
        self.minPower = minPower
        self.extinguishingAngle = extinguishingAngle
        print('Setting up the Thorlabs K10CR1:')
        bits, version = platform.architecture()
        print('    Detected %s Python on %s. Loading %s DLLs' % (bits, version, bits))
        
        dllname = os.path.join(os.path.dirname(__file__), 'dll%s' %
                               bits[:2], 'Thorlabs.MotionControl.IntegratedStepperMotors.dll')
        
        os.environ['PATH'] = os.environ['PATH'] + ';' + \
            os.path.join(os.path.dirname(__file__), 'dll%s' % bits[:2])
        
        if not os.path.exists(dllname):
            raise ValueError('DLL Not found! dllname=%s' % dllname)
        
        if bits == '32bit':
            self.p = c.CDLL(dllname)  # Alternate between dll loading method
        else:
            self.p = c.windll.LoadLibrary(dllname)
    
        try:
            serialNumber = self.getDeviceList()[0]
        except:
            raise ValueError(
                'Couldn\'t get the list of serial numbers! Is your stage plugged in?',
                ' Or is Thorlabs Kinesis/APT open?')
    
        self.SN = c.c_buffer(serialNumber.encode('UTF-8'))
        print('    Stage found! Serial number %s'%serialNumber)
    
        try:
            self.p.ISC_Close(self.SN)
            print('    Previous stage connection closed.')
        except:
            pass
    
        self.p.ISC_Open(self.SN)
        print('    New stage connection opened.')
            
        hardwareinfoval = self.getHardwareInfo()
        self.p.ISC_StartPolling(self.SN, c.c_int(20))

        # Calculate the conversion between "Device units" and degrees
        stepsPerRev, gearBoxRatio, pitch = self.getMotorParamsExt()
        # from https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=8750
        microstepsPerFullstep = 2048
        self.conversion = stepsPerRev * microstepsPerFullstep * \
            gearBoxRatio / pitch  # convert to degrees
        # conversion is in "Device units per degree"

        # Begin moving stage to Home- defines where zero is
        print('\nHoming...', end='')
        sys.stdout.flush()
        self.p.ISC_Home(self.SN)
        time.sleep(0.2)

        # While Motor is moving.  Stopped is -2147482624
        while (self.p.ISC_GetStatusBits(self.SN)) != (-2147482624):
            # print(p.ISC_GetStatusBits(SN))
            # print 'in the process of homing'
            print('.', end='')
            sys.stdout.flush()
            time.sleep(0.5)

        print('   Done! K10CR1 is ready.\n')
        

    def getHardwareInfo(self):
        modelNo = c.c_buffer(255)
        sizeOfModelNo = c.c_ulong(255)
        hardwareType = c.c_ushort()
        numChannels = c.c_short()
        notes = c.c_buffer(255)
        sizeOfNotes = c.c_ulong(255)
        firmwareVersion = c.c_ulong()
        hardwareVersion = c.c_ushort()
        modState = c.c_ushort()
        # p.PCC_GetHardwareInfo(SN)
    
        self.p.ISC_GetHardwareInfo(self.SN,
                              c.pointer(modelNo),
                              c.pointer(sizeOfModelNo),
                              c.pointer(hardwareType),
                              c.pointer(numChannels),
                              c.pointer(notes),
                              c.pointer(sizeOfNotes),
                              c.pointer(firmwareVersion),
                              c.pointer(hardwareVersion),
                              c.pointer(modState))
    
        return [x.value for x in (modelNo, sizeOfModelNo, hardwareType,
                                  numChannels, notes, sizeOfNotes, firmwareVersion,
                                  hardwareVersion, modState)]
    
    
    def getMotorParamsExt(self):
        # p.ISC_ClearMessageQueue(SN);
    
        stepsPerRev = c.c_double()
        gearBoxRatio = c.c_double()
        pitch = c.c_double()
    
        self.p.ISC_GetMotorParamsExt(self.SN, c.pointer(stepsPerRev),
                                      c.pointer(gearBoxRatio),
                                      c.pointer(pitch))
    
        if stepsPerRev.value < 1 or gearBoxRatio.value < 1 or pitch.value < 1:
            print('    Failed to get motor params, using default values!')
            print('        stepsPerRev=200, gearBoxRatio=120, pitch=360')
    
            return 200, 120, 360
        else:
            return stepsPerRev.value, gearBoxRatio.value, pitch.value
    
    
    def getDeviceList(self):
        self.p.TLI_BuildDeviceList()
        receiveBuffer = c.c_buffer(200)
        sizeOfBuffer = c.c_ulong(255)
        self.p.TLI_GetDeviceListExt(c.pointer(receiveBuffer), c.pointer(sizeOfBuffer))
        ser_num = [x.replace('b\'', '') for x in
                (str(receiveBuffer.value)).split(',')[:-1]]
        return ser_num
    
    
    def MoveToPosition(self, deviceUnits, timeout=20, queryDelay=0.01, tolerance=1):
        """
        Moves the rotation stage to a certain position (given by device units).
        This call blocks future action until the move is complete.
        The timeout is in seconds
    
        SN is a c_buffer of the serial number string
        deviceUnits shold be a int.
        tolerance is when the blocking should end (device units)
        """
    
        self.p.ISC_RequestStatus(self.SN)  # order the stage to find out its location
        self.p.ISC_MoveToPosition(self.SN, c.c_int(int(deviceUnits)))
    
        t = time.time()
    
        while time.time() < (t+timeout):
            self.p.ISC_RequestStatus(self.SN)  # order the stage to find out its location
            currentPosition = self.p.ISC_GetPosition(self.SN)
            error = currentPosition - deviceUnits
            if np.abs(error) < tolerance:
                return
            else:
                time.sleep(queryDelay)
        raise ValueError('Did not reach position!',
                         'Increase timeout from %.3f seconds?' % timeout)
    
    def write_header(self, logfile):
        logfile.write('Min_pow angle: %.4f\n'%self.extinguishAngle)
        logfile.write('Max Power: %.4f\n'%self.maxPower)
        logfile.write('Min Power: %.4f\n'%self.minPower)
        logfile.write('FileNum\tPower\t Angle (deg)\n')
            
    def set_power(self, power):
        if power>=maxPower or power<=minPower:
            raise ValueError('Power requested outside of range.')
    
        angle = self.angleFromPower(power)
        
        # convert the desired position to integer "Device units" to be passed to the stage
        # NOTE: this involves rounding, and could introduce errors, especially if you are making
        # steps of just a few device units.
        try:
            deviceUnits = abs(int(angle*self.conversion))  # -deviceUnitsZero)
        except ValueError:
            raise ValueError(('Could not get the position from the stage.',
                              ' This typically means that you need to unplug the',
                              ' stage and plug it back in. And restart your',
                              ' python terminal.'))

        self.MoveToPosition(deviceUnits)
        
        return angle

    def angleFromPower(self, power):
        return -(np.arcsin( ((power-self.minPower)/(self.maxPower-self.minPower))**0.5))*(180/np.pi)/2.  + self.extinguishingAngle
    
    def close(self):
        print('Moving Back to Max Power')
        self.MoveToPosition(abs(int((self.extinguishingAngle+45)*self.conversion)))
        print('Power scan complete!')
        self.p.ISC_Close(self.SN)
        


if __name__ == '__main__':
    main()

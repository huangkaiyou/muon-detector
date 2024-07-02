import pyvisa as pv
import time


##### find and get the information of the oscilloscope #####
rm = pv.ResourceManager()
de = rm.list_resources()     # find the address of the devices that connected to the rpi

inst = rm.open_resource(f"{de[0]}") # open the oscilloscope


print("Initializing...")

##### initialize the oscilloscope #####
trigger_source = 2

def initialize():
    length = 1e+4               # record length
    scale = 50e-9              # timescale
    pos = scale * 4          # tiembase position
    hold = scale                # holdoff
    level = 1               # trigger level
    
    
    
    inst.write(":ACQuire:MODe SAMPle")               # set the acquisition mode to "sample mode"
    inst.write(f":ACQuire:RECOrdlength {length}")    # set the record length
    inst.write(":TRIGger:MODe NORMal")               # triger mode set as normal
    
    inst.write(f":TIMebase:SCALe {scale}")           # set the time scale
    inst.write(f":TIMebase:POSition {pos}")          # set the postion of the timebase
    
    inst.write(f":TRIGger:HOLDoff {hold}")           # set the holdoff 
    inst.write(f":TRIGger:LEVel {level}")            # set the triger level
    inst.write(f":TRIGger:SOURce CH{trigger_source}")                # set the trigger source at channel 2
    
    inst.write(":SAVe:WAVEform:FILEFormat FCSV")     # set the save fileformat as fast CSV
    
   
    
##### set up the channel 1 and channel 2 #####
def setup():
    scale = 2          # vertical scale
    pos = -scale * 3.5               # vertical position

    
    
    inst.write(f":CHANnel1:POSition {pos}")         # set vertical position of the channel 
    inst.write(f":CHANnel2:POSition {pos}")         # set vertical position of the channel 

    inst.write(f":CHANnel1:SCALe {scale}")         # set vertical scale of the channel 1 
    inst.write(f":CHANnel2:SCALe {scale}")         # set vertical scale of the channel 2 


##### read the waveform data #####
initialize()
setup()


inst.write(":SINGle")
t0 = time.time()   # the start time as reference
duration = 60 * 60 * 2     # how long we wanna record data (second/minunte/hour)

print("Reading...")
time.sleep(0.1)

i=0
while time.time() - t0 <= duration:
    try:
        if inst.query(f":ACQuire{trigger_source}:STATe?") == '1\n':
            
            i+=1
            inst.query(f":SAVe:WAVEform ALL, 'USB:/energy_dis/raw{i}.CSV'")
            time.sleep(0.01)
            inst.write(":SINGle")
            time.sleep(0.01)
            
    except:
        inst.write(":SINGle")
        time.sleep(0.01)
        pass
  
print("Finish")
          




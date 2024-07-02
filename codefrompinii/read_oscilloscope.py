import pyvisa as pv
import time
import pandas as pd


##### find and get the information of the oscilloscope #####
rm = pv.ResourceManager()
de = rm.list_resources()     # find the address of the devices that connected to the rpi

inst = rm.open_resource(f"{de[0]}") # open the oscilloscope
inst.timeout = 30000


##### initialize the oscilloscope #####
def initialize():
    length = 1e+4               # record length
    scale = 1.0e-6              # timescale
    pos = scale * 5          # tiembase position
    hold = scale                # holdoff
    level = 0.5               # trigger level
    
    
    
    inst.write(":ACQuire:MODe SAMPle")               # set the acquisition mode to "sample mode"
    inst.write(f":ACQuire:RECOrdlength {length}")    # set the record length
    inst.write(":TRIGger:MODe NORMal")               # triger mode set as normal
    
    inst.write(f":TIMebase:SCALe {scale}")           # set the time scale
    inst.write(f":TIMebase:POSition {pos}")          # set the postion of the timebase
    
    inst.write(f":TRIGger:HOLDoff {hold}")           # set the holdoff 
    inst.write(f":TRIGger:LEVel {level}")            # set the triger level
    inst.write(":TRIGger:SOURce CH2")                # set the trigger source at channel 2
    
   
    
##### set up the channel 1 and channel 2 #####
def setup():
    pos = 0               # vertical position
    scale = 5e-1          # vertical scale
    
    
    
    inst.write(f":CHANnel1:POSition {pos}")         # set vertical position of the channel 
    inst.write(f":CHANnel2:POSition {pos}")         # set vertical position of the channel 

    inst.write(f":CHANnel1:SCALe {scale}")         # set vertical scale of the channel 1 
    inst.write(f":CHANnel2:SCALe {scale}")         # set vertical scale of the channel 2 


##### read the waveform data #####
initialize()
setup()


inst.write(":SINGle")
t0 = time.time()   # the start time as reference
duration = 60      # how long we wanna record data

i = 1
print("start")
while time.time() - t0 <= duration:
    try:
        if inst.query(":ACQuire2:STATe?") == '1\n':
        
            inst.query(":ACQuire1:MEMory")   # get the waveform information of the channel 1
            inst.read_bytes(7)               # get rid off the information about the data
            col1 = inst.read_bytes(2e4)     # save the acutal data from channel 1 into a list
        
        
            inst.query(":ACQuire2:MEMory")   # get the waveform information of the channel 2
            inst.read_bytes(7)               # get rid off the information about the data
            col2 = inst.read_bytes(2e4)     # save the acutal data from channel 2 into a list
        
            data = {"ch1" : col1,
                    "ch2" : col2
                    }
        
            df = pd.DataFrame(data)
            df.to_csv(f"/home/chip/Desktop/data/rawdata{i}.csv")
        
            i+=1
            inst.write(":SINGle")
            print("1")
    except:
        inst.write(":SINGle")
        print("0")
        pass
        


  
print("Finish")
          



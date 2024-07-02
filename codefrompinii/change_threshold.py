import pyvisa as pv
import time
import pandas as pd


##### find and get the information of the oscilloscope #####
rm = pv.ResourceManager()
de = rm.list_resources()     # find the address of the devices that connected to the rpi

inst = rm.open_resource(f"{de[0]}") # open the oscilloscope


print("Initializing...")

##### initialize the oscilloscope #####
trigger_source = 1    # the channel as the trigger source

def initialize():
    length = 1e+4               # record length
    scale = 1.0e-6              # timescale
    pos = scale * 4          # tiembase position
    hold = scale                # holdoff
    level = 0.3               # trigger level
    
    
    
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
    scale = 0.1          # vertical scale
    pos = -scale * 3               # vertical position

    
    
    
    inst.write(f":CHANnel1:POSition {pos}")         # set vertical position of the channel 
    inst.write(f":CHANnel2:POSition {pos}")         # set vertical position of the channel 

    inst.write(f":CHANnel1:SCALe {scale}")         # set vertical scale of the channel 1 
    inst.write(f":CHANnel2:SCALe {scale}")         # set vertical scale of the channel 2 


##### measure the counts of each threshold #####
initialize()
setup()


duration = 60 * 20     # how long we wanna record data (second/minunte/hour)


threshold = [0.001*k for k in range(50, 520, 50)]    # the threshold we wanna test
counts = []    # the counts of the events of each set of threshold 

print("Measuring...")


for i in range(len(threshold)): 
    inst.write(f":TRIGger:LEVel {threshold[i]}")    # set the triger level
    inst.write(":SINGle")
    i=0     # count the number of events 
    t0 = time.time()      # the reference time
    duration = 60 * 20     # the duration of the measurement of each threshold (second/minuite/hour)

    print("Measuring the threshold = " + str(threshold[i]))

    while time.time() - t0 <= duration:
        try:
            if inst.query(f":ACQuire{trigger_source}:STATe?") == '1\n':
                a = time.time() - t0
                b = 100 * a//duration
                print("Process: " + str(b) + "%")
                inst.write(":SINGle")
                i+=1
                
        except:
            inst.write(":SINGle")
            pass
    print("Counts: " + str(i))
    counts.append(i)


data = {
        "Threshold": threshold,
        "Counts": counts
        }


df = pd.DataFrame(data)
df.to_csv("/home/chip/Desktop/data/counts_threshold_new.csv")
 
print("Finish")
          




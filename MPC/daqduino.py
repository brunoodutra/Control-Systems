# Bruno Dutra

import matplotlib.pyplot as plt
import serial
import time
import numpy as np
var=0.0;
def start(port,baud):
    port=str(port)
    global ser
    ser= serial.Serial(port,baud, timeout=None, xonxoff=False, rtscts=False, dsrdtr=False)
    
    #2000000

    time.sleep(2.5)
    ser.write(bytearray(str(0),'ascii'))
    print("Connected!!")
    
def write(u,Ts):
  #for u in val:
    #print(u)
  u=np.round(u,2)
  ser.write(bytearray(str(u),'ascii'))
  time.sleep(Ts);
   
def read():
        Str=ser.readline().rstrip().decode();
        STR=Str.split()
        var=np.zeros(len(STR))
        for i,Read in enumerate(STR):
            var[i]=float(Read)
        ser.flushInput()
        # #print(Str)
        # global var
        # if len(Str)>4:
        #     try:
        #          v=Str
        #          var=float(v)
        #     except ValueError:
        #          print ("Oops! erro")
     
           
        return(var)
def readStr():

      try:
         Str=str((ser.readline().rstrip()),'ascii');

      except ValueError:
        print ("Oops! error")
        Str=" "
        
      return(Str)
def fastRead():

      try:
          bytesToRead = ser.inWaiting()
         # print (bytesToRead)
          Str=str((ser.read(191).rstrip()),'ascii');
        # Str=str((ser.readline().rstrip()),'ascii');
        
      except ValueError:
        print ("Oops! error")
        Str=" "

      return(Str)    
    
def end():
    ser.write(bytearray(str(0),'ascii'))
    ser.close();  
    print("Serial Connection Finished!!")

     

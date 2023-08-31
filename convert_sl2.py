import numpy as np
import struct
import os 
import csv

"""
Convert Lowrance sl2 files to lat, lon, waterdepth csv file
Reference from R, Ruby, js code:
    R: https://gitlab.com/hrbrmstr/arabia/-/blob/master/R/read-sl2.R?ref_type=heads (index starts at 1)
    Ruby: https://github.com/Chris78/sl2decode/blob/master/sl2decode.rb (index starts at 0)
    js: https://github.com/KW-M/sl2-csv-converter/blob/master/ruby_conversion.js
"""

#   name of block : {offset in block, type, length}
#       type: byte=UInt8="C"
#             short=Uint16LE="v" (integer, precision of 15 bits)
#             int=UInt32LE="V" (integer, precision of 31 bits)
#             float=FloatLE (32 bits IEEE754 floating point num)="e" (single precision floating-point)

#header = first 10 bytes
blockOffset = 10 #offset header

info = {'current blocksize' : {26 , "v", 2} , 
        'previous blocksize' : {28 , "v", 2} ,
        'waterDepth' : {62, "e", 4}, 
        'longitude' : {106, 'V', 4}, 
        'latitude' : {110, "V", 4}}

#for lon, lat, depth conversion
POLAR_EARTH_RADIUS = 6356752.3142
MAX_UINT4 = 4294967295
FT2M = 1/3.2808399

def bytes_to_int(n, chunk, seek, add):
    """
    converts bytes to integers
    input:
        n, how many bytes look at (don't really need this, all are set to 1)
        chunk, how many bytes to read in
        seek, where the info begins in the byte string block
        add, number of bytes to add for the next block
    output:
        integer
    """
    r = []
    N = seek + add
    file.seek(N) #this advances N bytes into byte stream
    byte = file.read(chunk)
    count = 0
    while count < n and byte:
        integ = int.from_bytes(byte, byteorder='little')
        r.append(integ)
        byte = file.read(chunk)
        count += 1
    return r #returns list, need to specify element when operating on value

def bytes_to_float(n, chunk, seek, add):
    """
    converts bytes to float
    input:
        n, how many bytes look at (don't really need this, all are set to 1)
        chunk, how many bytes to read in
        seek, where the info begins in the byte string block
        add, number of bytes to add for the next block
    output:
        float
    """
    r = []
    N = seek + add
    file.seek(N) #this advances N bytes into byte stream
    byte = file.read(chunk)
    count = 0
    while count < n and byte:
        integ = struct.unpack('<f', byte)[0]
        r.append(integ)
        byte = file.read(chunk)
        count += 1
    return r #returns list, need to specify element when operating on value

def convertLon(lon):
    """
    input:
        lon, raw longitude from bytes function
    output:
        conventional longitude value 
    """
    r = lon[0] - MAX_UINT4
    return r/POLAR_EARTH_RADIUS * (180/np.pi)

def convertLat(lat):
    """
    input:
        lat, raw latitude from bytes function
    output:
        conventional latitude value
    """
    return ((2*np.arctan(np.exp(lat[0]/POLAR_EARTH_RADIUS)))-(np.pi/2)) * (180/np.pi)

filename = '/Users/mtan/Desktop/MT_cores/Sonar_2023-07-30_12.49.06.sl2' #sl2 filepath
project = 'Sonar_2023-07-30_12.49.06.sl2_conv3' #name the converted csv file
filesize = os.path.getsize(filename) #get size of file in bytes
file = open(filename, "rb") #read binary file using 'rb'

blockSize = bytes_to_int(1, 2, 4, 0)[0] #determine if sl2 file; if not, code terminates
if blockSize == 3200:
    print("Blocksize: Sidescan")
elif blockSize == 1600:
    print("Blocksize: Downscan")
else:
    print("ERROR!  Blocksize: N/A; this file may be corrupted and/or not an sl2 file")
    exit()

with open("/Users/mtan/Desktop/MT_cores/{}.csv".format(project),'a+') as f: #create or append to csv file
    fieldnames = ['lat', 'lon', 'waterDepth']
    writeHeader = csv.DictWriter(f, fieldnames=fieldnames)
    writeHeader.writeheader()
    writer = csv.writer(f)
    counter = 0
    updateBlock = 0

    prevLon = 0
    prevLat = 0

    while updateBlock < filesize:
        try:
            currentBlock = bytes_to_int(1, 2, 26+blockOffset, updateBlock)[0]
        except:
            #throws back IndexError when it gets to the last block. Probably better way to implement...
            print("End of byte string, sl2 conversion complete.")
            exit()

        previousBlock = bytes_to_int(1, 2, 28+blockOffset, updateBlock)[0]

        depth = bytes_to_float(1, 4, 62+blockOffset, updateBlock)[0] #water depth [ft]
        lon = bytes_to_int(1, 4, 106+blockOffset, updateBlock) #longitude
        lat = bytes_to_int(1, 4, 110+blockOffset, updateBlock) #latitude

        #avoid writing duplicate data or data with depth=0
        if (prevLat != lat[0] and prevLon != lon[0]) or depth != 0:  
            #lat and lon conversions
            longitude = convertLon(lon)
            latitude = convertLat(lat)

            #write to csv file
            attr = [latitude, longitude, depth*FT2M]
            writer.writerow(attr)

            #update checker
            prevLat = lat
            prevLon = lon
        else:
            pass

        #tell code how to advance in byte string
        if counter == 0:
            updateBlock = 0
        elif counter == 1:
            updateBlock = 3216
        else:
            updateBlock += currentBlock

        counter += 1

#uncomment the below code (and comment above code) to test code on snippet of sl2 file
"""
counter = 0
updateBlock = 0
while counter < 30: #prints out first 30 blocks
    try:
        currentBlock = bytes_to_int(1, 2, 26+blockOffset, updateBlock)[0]
    except:
        print("End of byte string, sl2 conversion complete.")
        exit()

    previousBlock = bytes_to_int(1, 2, 28+blockOffset, updateBlock)[0]

    depth = bytes_to_float(1, 4, 62+blockOffset, updateBlock)[0] #water depth [ft]
    lon = bytes_to_int(1, 4, 106+blockOffset, updateBlock) #longitude
    lat = bytes_to_int(1, 4, 110+blockOffset, updateBlock) #latitude

    #lat and lon conversions
    longitude = convertLon(lon)
    latitude = convertLat(lat)

    print(latitude, longitude, depth)

    #tell code how to advance in byte string
    if counter == 0:
        print("0- block, disregard")
        updateBlock = 0
    elif counter == 1:
        print("1- block")
        updateBlock = 3216
    else:
        print("{}- block".format(counter))
        updateBlock += currentBlock

    print("BLOCK PREV=", previousBlock)
    print("BLOCK NEXT=", currentBlock)
    print("BLOCK SIZE=", updateBlock)
    print("--------------------------")
    counter += 1
"""
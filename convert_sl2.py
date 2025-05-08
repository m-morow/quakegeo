import numpy as np
import struct
import os 
import csv
import sys

"""
Original Author: M. Morow Tan, M.Tan@colorado.edu

Convert Lowrance sl2 files to lat, lon, waterdepth csv file
Reference from R, Ruby, js code:
    R: https://gitlab.com/hrbrmstr/arabia/-/blob/master/R/read-sl2.R?ref_type=heads (index starts at 1)
    Ruby: https://github.com/Chris78/sl2decode/blob/master/sl2decode.rb (index starts at 0)
    js: https://github.com/KW-M/sl2-csv-converter/blob/master/ruby_conversion.js

Updating Author: Jocelyn Reahl, jocelyn.reahl@colorado.edu
Last updated: 4/1/2025
"""

#   name of block : {offset in block, type, length}
#       type: byte=UInt8="C"
#             short=Uint16LE="v" (integer, precision of 15 bits)
#             int=UInt32LE="V" (integer, precision of 31 bits)
#             float=FloatLE (32 bits IEEE754 floating point num)="e" (single precision floating-point)

def bytes_to_int(file, n=1, chunk=None, seek=None, add=None):
    """
    Converts bytes to integers.
    
    Parameters
    ----------
    file : _io.BufferedReader
        Opened sl2 file to convert bytes to int.
    
    n : int, optional
        How many bytes look at. The default is 1.
        
    chunk : int, optional
        How many bytes to read in. The default is None.
        
    seek : int, optional
        Where the info begins in the byte string block. The default is None.
        
    add : int, optional
        Number of bytes to add for the next block. The default is None.
    
    Returns
    -------
    r : list of int
        List of converted values.
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


def bytes_to_float(file, n=1, chunk=None, seek=None, add=None):
    """
    Converts bytes to floats.
    
    Parameters
    ----------
    file : _io.BufferedReader
        Opened sl2 file to convert bytes to float.
    
    n : int, optional
        How many bytes look at. The default is 1.
        
    chunk : int, optional
        How many bytes to read in. The default is None.
        
    seek : int, optional
        Where the info begins in the byte string block. The default is None.
        
    add : int, optional
        Number of bytes to add for the next block. The default is None.
    
    Returns
    -------
    r : list of float
        List of converted values.
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

def convertLon(lon, MAX_UINT4=4294967295, POLAR_EARTH_RADIUS=6356752.3142):
    """
    Parameters
    ----------
    lon : list of int
        Raw longitude output from bytes_to_int() function.
    
    MAX_UINT4 : int, optional
        Maximum integer to subtract from lon for conversion.
        The default is 4294967295.
    
    POLAR_EARTH_RADIUS : float, optional
        The polar radius of the Earth in meters.
        The default is 6356752.3142.
    
    Returns
    -------
    lon_deg : float
        Converted longitude value in decimal degrees.
    """
    r = lon[0] - MAX_UINT4
    lon_deg = r/POLAR_EARTH_RADIUS * (180/np.pi)
    return lon_deg


def convertLat(lat, POLAR_EARTH_RADIUS=6356752.3142):
    """
    Parameters
    ----------
    lat : list of int
        Raw latitude output from bytes_to_int() function.
    
    POLAR_EARTH_RADIUS : float, optional
        The polar radius of the Earth in meters.
        The default is 6356752.3142.
    
    Returns
    -------
    lat_deg : float
        Converted latitude value in decimal degrees.
    """
    lat_deg = (((2*np.arctan(np.exp(lat[0]/POLAR_EARTH_RADIUS)))-(np.pi/2))\
               * (180/np.pi))
    return lat_deg


def check_blockSize(file, n=1, chunk=2, seek=4, add=0):
    '''
    Checks if the file is an sl2 file. Prints the blocksize style if an sl2
    file (e.g. Sidescan or Downscan) and exits the run if the file is not an
    sl2 file.

    Parameters
    ----------
    file : _io.BufferedReader
        Opened file to to check.
        
    n : int, optional
        How many bytes look at. The default is 1.
        
    chunk : int, optional
        How many bytes to read in. The default is 2.
        
    seek : int, optional
        Where the info begins in the byte string block. The default is 4.
        
    add : int, optional
        Number of bytes to add for the next block. The default is 0.

    '''
    blockSize = bytes_to_int(file, n=n, chunk=chunk, seek=seek, add=add)[0]
    if blockSize == 3200:
        print("Blocksize: Sidescan")
    elif blockSize == 1600:
        print("Blocksize: Downscan")
    else:
        print("ERROR!  Blocksize: N/A; this file may be corrupted and/or not an sl2 file")
        sys.exit()


def get_file(sl2_filename, cwd='', sl2_subfolder=''):
    '''
    Opens file for sl2 to csv conversion.

    Parameters
    ----------
    sl2_filename : str
        Name of the individual file.
        
    cwd : str, optional
        General file path to the main directory. The default is ''.
        
    sl2_subfolder : str, optional
        File path to a subfolder where the sl2 file is stored.
        The default is ''.

    Returns
    -------
    file : _io.BufferedReader
        Opened file.
    '''
    if cwd == '':
        cwd = os.getcwd()
    path = os.path.join(cwd, sl2_subfolder, sl2_filename)
    file = open(path, "rb") #read binary file using 'rb'
    return file


def get_filesize(sl2_filename, cwd='', sl2_subfolder=''):
    '''
    Gets size of the file in bytes.

    Parameters
    ----------
    sl2_filename : str
        Name of the individual file.
        
    cwd : str, optional
        General file path to the main directory. The default is ''.
        
    sl2_subfolder : str, optional
        File path to a subfolder where the sl2 file is stored.
        The default is ''.

    Returns
    -------
    filesize : int
        Size of the file in bytes.
    '''
    if cwd == '':
        cwd = os.getcwd()
    path = os.path.join(cwd, sl2_subfolder, sl2_filename)
    filesize = os.path.getsize(path) #get size of file in bytes
    return filesize


def sl2_to_csv(sl2_filename, csv_filename, cwd='',
               sl2_subfolder='', csv_subfolder='',
               blockOffset=10, POLAR_EARTH_RADIUS=6356752.3142,
               MAX_UINT4=4294967295, DEPTH_UNIT_CONVERSION=1/3.2808399):
    '''
    Converts an sl2 file to a csv file that contains latitude and longitude
    (in decimal degrees) and the water depth (in original depth units from
    sl2 file).

    Parameters
    ----------
    sl2_filename : str
        Name of the original sl2 file.
        
    csv_filename : str
        Name of the exported csv file.
        
    cwd : str, optional
        General file path to the main directory. The default is ''.
        
    sl2_subfolder : str, optional
        File path to a subfolder where the sl2 file is stored.
        The default is ''.
        
    csv_subfolder : str, optional
        File path to the subfolder where the csv file will be saved.
        The default is ''.
    
    blockOffset : int, optional
        Byte offset to use for conversions.
        The default is 10.
        
    POLAR_EARTH_RADIUS : float, optional
        The polar radius of the Earth in meters.
        The default is 6356752.3142.
        
    MAX_UINT4 : TYPE, optional
        Maximum integer to subtract from lon for conversion.
        The default is 4294967295.
        
    DEPTH_UNIT_CONVERSION : float, optional
        Constant to convert water depth from default feet (ft) to new unit.
        The default is 1/3.2808399, which converts feet to meters.

    '''
    file = get_file(sl2_filename, cwd=cwd, sl2_subfolder=sl2_subfolder)
    filesize = get_filesize(sl2_filename, cwd=cwd, sl2_subfolder=sl2_subfolder)
    check_blockSize(file)
    
    with open(os.path.join(cwd, csv_subfolder,
                           '{}.csv'.format(csv_filename)),
              'a+') as f:
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
                currentBlock = bytes_to_int(file, n=1, chunk=2,
                                            seek=26+blockOffset,
                                            add=updateBlock)[0]
            except:
                print("End of byte string, sl2 conversion complete.")
                file.close()
                break
            
            previousBlock = bytes_to_int(file, n=1, chunk=2,
                                         seek=28+blockOffset,
                                         add=updateBlock)[0]
            depth = bytes_to_float(file, n=1, chunk=4,
                                   seek=62+blockOffset,
                                   add=updateBlock)[0] #water depth [ft]
            lon = bytes_to_int(file, n=1, chunk=4,
                               seek=106+blockOffset,
                               add=updateBlock) #longitude
            lat = bytes_to_int(file, n=1, chunk=4,
                               seek=110+blockOffset,
                               add=updateBlock) #latitude
            
            #avoid writing duplicate data or data with depth=0
            if (prevLat != lat[0] and prevLon != lon[0]) or depth != 0:  
                #lat and lon conversions
                longitude = convertLon(lon, MAX_UINT4=MAX_UINT4,
                                       POLAR_EARTH_RADIUS=POLAR_EARTH_RADIUS)
                latitude = convertLat(lat,
                                      POLAR_EARTH_RADIUS=POLAR_EARTH_RADIUS)

                #write to csv file
                attr = [latitude, longitude, depth*DEPTH_UNIT_CONVERSION]
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

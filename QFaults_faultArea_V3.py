##############################################
## Generate surface-projected fault planes  ##
## from Quaternary Fault and Fold Database  ##
## mtan@usgs.gov, edited: 23 May 2024       ##
## to the user: run code in QGIS Python w/  ##
##              fault .shp layer selected   ##
##############################################

# Import packages (MT running on Python 3.11.3 default kernel)

import matplotlib.pyplot as plt
import numpy as np
import re
import pandas as pd
import csv
import geopandas as gpd
import shapely
from shapely.geometry import Polygon, MultiPolygon

# Functions

def getDip(slip):
    """
    determines dip degree of 'slip_sense' attribute
    input: 
        slip (str)
    output: 
        standard dip [degrees]
    """
    slip = str(slip) #6/30/23: TypeError: expected string or bytes-like object error
    if slip == 'Normal':
        return 50 
    elif slip == 'Reverse*':
        return 60 #05/08/24 updated from 45
    elif slip == 'Thrust*':
        return 50 #05/08/24 updated from 15
    elif slip == 'Strike slip' or re.match('lateral', slip):
        return 90 
        #nest this statement below 'Reverse*', because there are 'Reverse; Left lateral' strings
        #this statement catches 'Left lateral; Normal', 'Left lateral', 'Right lateral'
    elif slip == 'Unspecified*' or 'NULL' or 'No data':
        return 90 #maybe set as really high number so that later functions can flag it.....
    else:
        return 90 #TODO?: create some sort of flag for things that have no slip sense

def getLength(depth, dip):
    """
    calculates how far the fault plane extends at the surface 
    input: 
        depth [km]
        dip [degrees]
    output:
        horizontal extent of fault plane, to be projected to surface [degrees]
    """
    if type(dip) == str: 
        return 90
    else:
        dip = 90-dip
        return depth * np.tan(np.radians(dip))
    
def getVert(geom):
    """
    converts QGIS geom line feature to individual points
    input:
        QGIS feature of lines?
    output:
        x, latitude coordinates of points that make up the feature line
        y, longitude coordinates of points that make up the feature line
    """
    x = []
    y = []
    if QgsWkbTypes.isSingleType(geom.wkbType()):
        # single
        for pnt in geom.asPolyline():
            return pnt.x(), pnt.y()
    else:
        # multipart
        #all qfaults will likely be MultiPolylines, including cases to be safe
        for part in geom.asMultiPolyline():
            for pnt in part:
                y = np.append(pnt.y(), y)
                x = np.append(pnt.x(), x)
    return x, y

def getShift(dir, ll):
    """
    determines how to shift points in lat and lon at the surface to represent a dip direction 
    input:
        dir, direction of dip (str)
        ll, horizontal extent of fault plane projected to surface
    output:
        lat, how much to shift original latitude point
        lon, how much to shift original longitude point
    """
    dir = str(dir)
    if dir == '?nspecified' or dir == 'NULL':
        return 0, 0
    elif re.match("V+", dir):
        return 0, 0 
    elif re.match("N+", dir):
        if dir == 'NE*':
            return ll, ll * np.tan(np.radians(45))
        elif dir == 'NW*':
            return -ll, ll * np.tan(np.radians(45))
        elif dir == 'N':
            return 0, ll #for N
        else:
            return 0, 0
    elif re.match("S+", dir):
        if dir == 'SE':
            return ll * np.tan(np.radians(45)), -ll
        elif dir == 'SW*':
            return -ll * np.tan(np.radians(45)), -ll
        else:
            return 0, -ll #for S
    elif dir == 'E*':
        return ll, 0
    elif dir == 'W':
        return -ll, 0
    else: #this catches the singular 'Center' feature
        return 0, 0 
    
def draw(list1, list2):
    """
    sorts elements in lists so that polygon is drawn in the correct order
    input:
        list1, list of original coordinates (getVert), which are ordered by QGIS
        list2, list of shifted coordinates (getVert + getShift)
    output:
        list of sorted values to draw polygon around points
        i.e.: 
            list1 = [0, 1, 2]
            list2 = [3, 4, 5]
            lt1 = [0]
            ltR = [2, 1, 0]
            returns [0, 3, 4, 5, 2, 1, 0]

            note: first and last element of returned list must be the same to close the polygon
    """
    lt1 = list1[0]
    ltR = list1[::-1]
    #uncomment below for arrays
    #return np.insert((np.append(list2, ltR)), 0, lt1)
    return [lt1] + (list2 + ltR)

def combine(a, b):
    """
    combines latitude and longitude coordinates from two lists
    input:
        a, latitude array
        b, longitude array
    output:
        list of combined lat/lon coordinates
    """
    coord = []
    for i, j in zip(a, b):
        coord.append( [i, j] )
    return coord

def locFormat(s):
    """
    reformats elements within a dataframe column
    input:
        elements in a dataframe column
    output:
        elements as floats without double brackets, splits lists
    """
    s = s.replace('[[', '[').replace('[', '').replace(']', '').replace(']]', ']').split(',')
    s = [float(x) for x in s]
    return s

#initialize arrays
coorID = []
cc = []
count = 0

#operate on QGIS active layer
layer = iface.activeLayer()

#how to name the .csv and .shp files
project = '050924-QF'

#create or append to csv files
with open("/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}.csv".format(project),'a+') as f, open("/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}_buffer.csv".format(project), 'a+') as f2:
    fieldnames = ['fault_id', 'loc']
    #fieldnames2 = ['fid', 'fault_name', 'fault_id']

    writeHeader = csv.DictWriter(f, fieldnames=fieldnames)
    writeHeader.writeheader()
    writer = csv.writer(f)

    writeHeader2 = csv.DictWriter(f2, fieldnames=fieldnames)
    writeHeader2.writeheader()
    writer2 = csv.writer(f2)

    #loop through features in the active layer
    buffCount = -1
    regCount = -1
    for idx, feature in enumerate(layer.getFeatures()):
        
        #gets points that make up fault lines
        geom = feature.geometry()

        #gets dip direction, slip sense from feature attribute
        dir = feature['dip_direct']
        slip = feature['slip_sense']
        faultName = feature['fault_name']
        faultId = feature['fault_id']
        #slipRate = feature['slip_rate']
        #faultLength = feature['Shape_Leng']
        
        dip = getDip(slip)
        ll = 0.008 * getLength(15, dip)

        #coordinates, shifts for each feature
        xs, ys = getVert(geom)
        x_shift, y_shift = getShift(dir, ll)

        if x_shift == 0 and y_shift == 0:
            buffCount += 1
            xx = []
            yy = []
            buffer = geom.buffer(0.008*15, 2) #default 15km buffer on each side
            try:
                for part in buffer.asPolygon():
                    for pnt in part:
                        yy = np.append(pnt.y(), yy)
                        xx = np.append(pnt.x(), xx)
                        
                list1 = xx.ravel().tolist()
                list2 = yy.ravel().tolist()

                xyPoly = combine(list1, list2)
                
                try:
                    writer2.writerow([faultId, xyPoly])

                    #attr = [buffCount, faultName, faultId]
                    #writer4.writerow(attr)
                except ValueError:
                    pass

            #FOR SOME REASON.... some multilinestrings are converted to multipolygons?
            # most are converted to polygons, which avoids this TypeError
            # can't figure out what determines the converted geometry   
            except TypeError:
                count += 1 #bundle 2 ~5; bundle 3 =57; bundle 4 =10
            
        else:
            regCount += 1

            indivX = np.array(xs) + x_shift
            indivY = np.array(ys) + y_shift

            #link points in original and in shifted faults
            #indivID = np.full((len(xs),), idx)
            #coorID.append(indivID)
                    
            #sort elements in list to draw polygon correctly
            xPoly = draw(list1=xs.ravel().tolist(), list2=indivX.ravel().tolist())
            yPoly = draw(list1=ys.ravel().tolist(), list2=indivY.ravel().tolist())
                    
            #combine latitude and longitude elements
            xyPoly = combine(xPoly, yPoly)
        
            #save each polygon coordinates to row in csv file
            try:
                writer.writerow([faultId, xyPoly])

                #attr = [regCount, faultName, faultId]
                #writer3.writerow(attr)
            except ValueError:
                pass

    del geom #I don't think I need this, but may need to figure out memory optimization later......

#read in dataframe to create polygons
df = pd.read_csv('/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}.csv'.format(project), index_col=None)

#check if dataframe is empty, if not make shapefile
if not df.empty:
    #Courtesy of Stack Overflow
    df['loc'] = df['loc'].apply(locFormat)
    df['LAT'] = df['loc'].apply(lambda x: x[1::2])
    df['LON'] = df['loc'].apply(lambda x: x[::2])
    geom_list = [(x, y) for x, y  in zip(df['LON'],df['LAT'])]
    geom_list_2 = [Polygon(tuple(zip(x, y))) for x, y in geom_list]
    polygon_gdf =  gpd.GeoDataFrame(geometry=geom_list_2, crs="EPSG:4326")
    polygon_gdf['fault_id'] = df['fault_id']

    #save shapefile to computer
    polygon_gdf.to_file('/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}.shp'.format(project))

#read in dataframe to create buffer polygon
df2 = pd.read_csv('/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}_buffer.csv'.format(project), index_col=None)

#check if dataframe is empty, if not make shapefile
if not df2.empty:
    #Courtesy of Stack Overflow
    df2['loc'] = df2['loc'].apply(locFormat)
    df2['LAT'] = df2['loc'].apply(lambda x: x[1::2])
    df2['LON'] = df2['loc'].apply(lambda x: x[::2])
    geom_list2 = [(x, y) for x, y  in zip(df2['LON'],df2['LAT'])]
    geom_list_3 = [Polygon(tuple(zip(x, y))) for x, y in geom_list2]
    polygon_gdf2 =  gpd.GeoDataFrame(geometry=geom_list_3, crs="EPSG:4326")
    polygon_gdf2['fault_id'] = df['fault_id']

    polygon_gdf2.to_file('/Users/mtan/Library/CloudStorage/OneDrive-DOI/data/Mw6_surfRup_onOffFault/output/QF/AZ/{}_buffer.shp'.format(project))
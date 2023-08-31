import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from decimal import Decimal

"""
Calculate moment and moment rate for magnitude groups
"""

project = 'NSHM23Catalog_onFault_NSHMPolys_WardPolys' #csv filename (on fault events csv)
refProject = 'NSHM23Catalog_onShore_WardPolys' #csv filename (all events csv)
region = 'SC' #regions in Ward polygon

print("Calculating moment for Ward 1998 {} region...".format(region))

#read in dataframe
df = pd.read_csv('/Users/mtan/Documents/Data/WardPolygons/{}.csv'.format(project), index_col=None)
df2 = pd.read_csv('/Users/mtan/Documents/Data/WardPolygons/{}.csv'.format(refProject), index_col=None)

#filter dataframe to only include region of interest
df_Ward = df[df['region'].str.contains(region)]
df2_Ward = df2[df2['region'].str.contains(region)]

beginYear = df_Ward['Year'].min()
endYear = df_Ward['Year'].max()
totYear = endYear - beginYear

#convert dataframe rows to lists of moment, grouped by magnitude
#currently hardcoded
class cls:
    x3 = list(df_Ward[df_Ward['Mw'].between(3,4, inclusive='left')]['M0'])
    x4 = list(df_Ward[df_Ward['Mw'].between(4,5, inclusive='left')]['M0'])
    x5 = list(df_Ward[df_Ward['Mw'].between(5,6, inclusive='left')]['M0'])
    x6 = list(df_Ward[df_Ward['Mw'].between(6,7, inclusive='left')]['M0'])
    x7 = list(df_Ward[df_Ward['Mw'].between(7,8, inclusive='left')]['M0'])
    x8 = list(df_Ward[df_Ward['Mw'].between(8,9, inclusive='left')]['M0'])
    x9 = list(df_Ward[df_Ward['Mw'].between(9,10, inclusive='left')]['M0'])

class cls2:
    y3 = list(df2_Ward[df2_Ward['Mw'].between(3,4, inclusive='left')]['M0'])
    y4 = list(df2_Ward[df2_Ward['Mw'].between(4,5, inclusive='left')]['M0'])
    y5 = list(df2_Ward[df2_Ward['Mw'].between(5,6, inclusive='left')]['M0'])
    y6 = list(df2_Ward[df2_Ward['Mw'].between(6,7, inclusive='left')]['M0'])
    y7 = list(df2_Ward[df2_Ward['Mw'].between(7,8, inclusive='left')]['M0'])
    y8 = list(df2_Ward[df2_Ward['Mw'].between(8,9, inclusive='left')]['M0'])
    y9 = list(df2_Ward[df2_Ward['Mw'].between(9,10, inclusive='left')]['M0'])

for i in range(3, 10):
    onFaultEvents = getattr(cls, 'x{}'.format(i)) #call list variable using class
    allEvents = getattr(cls2, 'y{}'.format(i)) #call list variable using class
    onFaultMoment = sum(onFaultEvents)
    allMoment = sum(allEvents)
    offFaultMoment = allMoment - onFaultMoment
    print("M", i, "=")
    try: #use f-string and Decimal package to set sig figs
        #print("on fault =", f"{Decimal(onFaultMoment):.2E}", "Nm |", f"{Decimal(onFaultMoment/totYear):.2E}", "Nm/yr |", f"{Decimal((onFaultMoment/allMoment)*100):.2E}")
        print(f"{Decimal(onFaultMoment/totYear):.2E}")
    except ZeroDivisionError: #avoids error that stops code
        print("on fault = 0")
    try: #use f-string and Decimal package to set sig figs
        #print("off fault =", f"{Decimal(offFaultMoment):.2E}", "Nm |", f"{Decimal(offFaultMoment/totYear):.2E}", "Nm/yr |", f"{Decimal((offFaultMoment/allMoment)*100):.2E}")
        print(f"{Decimal(offFaultMoment/totYear):.2E}")
    except ZeroDivisionError: #avoids error that stops code
        print("off fault = 0")
    print("--")
    
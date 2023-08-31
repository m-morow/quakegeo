import pandas as pd

"""
Update length or quality of offset measurements in QGIS attribute table
change lines 18, 32, 34 based on what kind of length you want (i.e. min, max, avg/preferred)
if there is no value, default to '1000'
"""

attributes = []
layer = iface.activeLayer() #get attribute values from active layer

for idx, feature in enumerate(layer.getFeatures()):
    attrs = feature.attributes()
    strValue = str(attrs[0])
    attributes.append([strValue, attrs[1]])

df = pd.DataFrame(attributes, columns=['id', 'length'])
df = df[df["id"].str.contains("max")] #max, min, avg

layer = QgsProject.instance().mapLayersByName("optical_offset_measurements_b_UTM")[0]
iface.setActiveLayer(layer)

NULL_VALUE = 1000

with edit(layer):
    for idx, feature in enumerate(layer.getFeatures()):
        featureID = str(feature["id_2"])
        try:
            feature["max"] = float(df[df['id'].str.contains(featureID)]["length"].iloc[0]) #max, min, preferred
        except IndexError or ValueError:
            feature["max"] = NULL_VALUE #max, min, preferred

        layer.updateFeature(feature)


######################################################

"""

#update quality
#if there is no value, default to '0'

attributes = []
layer = iface.activeLayer() #get attribute values from active layer

for idx, feature in enumerate(layer.getFeatures()):
    attrs = feature.attributes()
    strValue = str(attrs[3])
    attributes.append([strValue, attrs[1]])

df = pd.DataFrame(attributes, columns=['id', 'quality'])

layer = QgsProject.instance().mapLayersByName("optical_offset_measurements_b")[0]
iface.setActiveLayer(layer)

NULL_VALUE = 0

with edit(layer):
    for idx, feature in enumerate(layer.getFeatures()):
        featureID = str(feature["id_2"])
        try:
            feature["quality"] = df[df['id'].str.contains(featureID)]["quality"].iloc[0]
        except IndexError or ValueError:
            feature["quality"] = NULL_VALUE

        layer.updateFeature(feature)

"""
import ee
ee.Initialize()
from ee.batch import Export
import geopandas as gpd
import sys

MAX_TASKS = 3000

def help():
    print(f'''
Extract zonal LTS from the Google Earth Engine (v0.1)
-----------------------------------------------------

Usage:
------
    python {sys.argv[0]} -i [input_shapfile] -o [output_file_prefix] {{-s [start_feature=0]}} {{-e [end_feature=3000]}} {{-f [target_FID])}} {{-d [description="LTS Export"]}} {{-F [folder=""]}}

Required arguments:
-------------------
    -i/--input        Input filename (shapfile). The shapefile must be comprised of single or multi-polygon geometries and have an attribute with the name "TARGET_FID".
    -o/--output       Output file prefix. The TARGET_FID and ".csv" will be appended automatically to each output file, so you do not need to give this information here.
    
Optional arguments:
-------------------
    -h/--help         Show this help message and exit.
    -s/--start        The feature in the input shapefile to start with. Numbering starts at 0, so the first feature will be 0.
    -e/--end          The feature in the input shapefile to end with. All features from 'start' to 'end' will be tasked for processing (1 task per feature).
    -c/--cloud        Cloud cover threshold (%). Only select images with cloud cover less than this (default = 100%; all images used)
    -f/--fid          A single FID to process. `start` and `end` will be ignored if this is given, and only one task will be created. The FID must be typed exactly as it is written in the TARGET_FID column of the input shapfile; otherwise, and error will be returned.
    -d/--description  The description to be shown in your GEE task list. The FID for that task will be automatically appended to the description, so you don't need to give that information here.
    -F/--folder       The output folder on your GoogleDrive where output files will be placed. By default, files will be placed directly in the top level of your Drive folder.

Examples:
---------
    To show this help message:
        python {sys.argv[0]} -h

    To extract LTS from the first 3000 features in a shapfile and set the output filename prefix to "exported_LTS" and save to a folder on your Drive named "gee_exports":
        python {sys.argv[0]} -i input.shp -o exported_LTS -F gee_exports

    To extract LTS from the first 10 features in a shapfile:
        python {sys.argv[0]} -i input.shp -o exported_LTS -s 0 -e 10 -F gee_exports

    To extract LTS from features 3000 to 6000 in a shapfile:
        python {sys.argv[0]} -i input.shp -o exported_LTS -s 3000 -e 6000 -F gee_exports

    To extract LTS from a single, specific TARGET_FID (e.g., "XXXX"):
        python {sys.argv[0]} -i input.shp -o exported_LTS -f XXXX -F gee_exports

Notes:
------
    - Try to avoid spaces in folder and filenames. If you can't avoid it, make sure you put them in quotes when entering them in the command line.
    - You can monitor ongoing tasks in your "task" window at https://code.earthengine.google.com/tasks
    - When an individual task is finished, it will output a .csv file to your GoogleDrive. By default, the file will be placed in the top level of your Drive folder. You can specify a folder using the `-F` argument (note: it is case sensitive, as `-f` relates to a target FID).
    - The output csv files will have columns corresponding to each non-thermal spectral band, as well as NDVI and NDMI (more spectral bands can be added to this script and the output csv's if you want). There is also a date column, allowing you to import each csv in R and run the time series segmentation code I shared with you earlier.
    - All spectral reflectance and index (NDVI, NDMI) values have been scaled by 10000. To get the true reflectance or index value (if that matters to you), divide by 10000
    - I believe that the GEE limits the number of concurrent tasks to 3000. If I am wrong about this, or if this policy changes in the future, just limit the # of features to be processed yourself. You may also change the MAX_TASKS variable in this script.
    - Sometimes, individual tasks may fail if there are not enough computing resources available at that time. If that happens, take note of the FID that failed and re-run it separately using the `-f` argument.
    ''')
    sys.exit(0)


## defaults
scan = True
description = ""
start = 0
end = None
infile = ""
folder = ""
outprefix = ""
description = "LTS Export"
cloud = 100
fid = None
folder = ""

## command-line arguments
i = 1
while scan:
    try:
        if sys.argv[i] in ['-h', '--help']:
            help()
        elif sys.argv[i] in ['-i', '--input']:
            infile = sys.argv[i+1]
            i += 2
        elif sys.argv[i] in ['-o', '--output']:
            outprefix = sys.argv[i+1]
            i += 2
        elif sys.argv[i] in ['-f', '--folder']:
            folder = sys.argv[i+1]
            i += 2
        elif sys.argv[i] in ['-d', '--description']:
            description = sys.argv[i+1]
            i += 2
        elif sys.argv[i] in ['-s', '--start']:
            start = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] in ['-e', '--end']:
            end = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] in ['-c', '--cloud']:
            cloud = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] in ['-f', '--fid']:
            fid = sys.argv[i+1]
            i += 2
        elif sys.argv[i] in ['-F', '--folder']:
            folder = sys.argv[i+1]
            i += 2
        else:
            raise ValueError(f"Unrecognized argument: {sys.argv[i]}")
    except IndexError:
        scan = False

if infile == "" or outprefix == "":
    raise ValueError("Input file and/or output file prefix missing.")

shp = gpd.read_file(infile).to_crs("EPSG:4326")
if fid is not None:
    shp = shp.query(f"TARGET_FID == '{fid}'")
    shp.reset_index(drop = True, inplace = True)
else:
    if not end:
        end = len(shp)
    elif end > len(shp):
        print(f"Warning: {len(shp)} features found. `end` set to {len(shp)}.")
        end = len(shp)

    if end - start > MAX_TASKS:
        end = start+MAX_TASKS
        print(f"Warning: Maximum of {MAX_TASKS} tasks can be started; limiting to features {start} to {end}.")

    shp.sort_values("TARGET_FID", inplace = True)
    shp.reset_index(drop = True, inplace = True)
    shp = shp[start:end]

## function to convert an individual features to a GEE geometry object
def gpdToEEGeom(i):
    poly = shp.loc[i,'geometry']
    if poly.geometryType() == "MultiPolygon":    
        poly = poly.geoms
    else:
        poly = [poly]
    
    geoms = []
    for p in poly:
        coordstring = p.wkt
        coordstring = coordstring\
            .replace("POLYGON ((", "")\
            .replace("))", "")\
            .replace("(", "")\
            .replace(")", "")\
            .split(",")
        xy = [c.lstrip(' ').split(' ') for c in coordstring]
        vertices = [(float(z[0]), float(z[1])) for z in xy]
        geoms.append(ee.Geometry.Polygon(ee.List(vertices)))
        
    if len(geoms) > 1:
        geoms = ee.Geometry.MultiPolygon(geoms)
    else:
        geoms = geoms[0]   
    return geoms

shp = shp.assign(ee_geom = [gpdToEEGeom(i) for i in shp.index])

## function to apply cloud mask to a single image
def maskImage(image):   
    qa = image.select("pixel_qa")
   
    cloudBit = ee.Number(2).pow(5).int()
    cloudShadowBit = ee.Number(2).pow(3).int()
    snowBit = ee.Number(2).pow(4).int()
    
    mask = ee.Image(0)\
        .where(qa.bitwiseAnd(cloudBit).neq(0), 1)\
        .where(qa.bitwiseAnd(cloudShadowBit).neq(0), 2)\
        .where(qa.bitwiseAnd(snowBit).neq(0), 3)\
        .updateMask(image.select('pixel_qa').mask())
    
    return image.updateMask(mask.eq(0))


## functions to compute NDVI, NDMI (and maybe other indices in the future) from L5/L7 and L8 images
def calcIndices(image):
    ndvi = image.normalizedDifference(['B4', 'B3']).multiply(10000).rename("NDVI")
    ndmi = image.normalizedDifference(['B4', 'B5']).multiply(10000).rename("NDMI")
    ### We can add more spectral indices here if needed ###
    return image\
        .select(
            ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'],
            ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
        )\
        .addBands([ndvi, ndmi])

def calcIndicesL8(image):
    ndvi = image.normalizedDifference(['B5', 'B4']).multiply(10000).rename("NDVI")
    ndmi = image.normalizedDifference(['B5', 'B6']).multiply(10000).rename("NDMI")
    ### We can add more spectral indices here if needed ###
    return image\
        .select(
            ['B2', 'B3', 'B4', 'B5', 'B6', 'B7'],
            ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
        )\
        .addBands([ndvi, ndmi])

print("Processing FIDs:")
for i in shp.index:
    print(shp.loc[i,'TARGET_FID'], ", ", sep = '', end = '', flush = True)
    
    filters = [
        ee.Filter.intersects('.geo', shp.loc[i,'ee_geom']),
    ]
    
    l8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")\
        .filter(filters)\
        .map(maskImage)\
        .map(calcIndicesL8)

    l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")\
        .filter(filters)\
        .map(maskImage)\
        .map(calcIndices)

    l5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")\
        .filter(filters)\
        .map(maskImage)\
        .map(calcIndices)

    landsat = ee.ImageCollection(l7.merge(l8).merge(l5)).sort('system:time_start')
    
    reducer_args = {
        'reducer': ee.Reducer.mean(),
        'scale': 30,
        'geometry': shp.loc[i,'ee_geom'],
    }

    def getMeanSR(image):
        meanBands = image.reduceRegion(**reducer_args)
        date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd')
        satellite = ee.String(image.get('SATELLITE'))
        meanBands = meanBands.set('date', date).set('satellite', satellite)
        return ee.Feature(None, meanBands)

    ts = landsat.map(getMeanSR)

    export_params = {
        'collection': ts,
        'description': f"{description} - {shp.loc[i,'TARGET_FID']}",
        'folder': folder,
        'fileNamePrefix': f"{outprefix}_{shp.loc[i,'TARGET_FID']}"
    }

    task = Export.table.toDrive(**export_params)

    task.start()
print("done.")
print(f"\n{len(shp)} tasks submitted.")
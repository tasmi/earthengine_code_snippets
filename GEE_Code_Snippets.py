# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 11:16:43 2020

@author: tsmith
"""

import ee
ee.Initialize()

#%% General Helper Functions
def run_export(image, crs, filename, scale, region, maxPixels=1e12):
    '''
    Runs an export function on GEE servers
    '''
    task_config = {'fileNamePrefix': filename,'crs': crs,'scale': scale,'maxPixels': maxPixels,'fileFormat': 'GeoTIFF','region': region,}
    task = ee.batch.Export.image.toDrive(image, filename, **task_config)
    task.start()

def gee_geometry_from_shapely(geom, ty='Polygon', crs='epsg:4326'):
    """ 
    Simple helper function to take a shapely geometry and a coordinate system and convert them to a 
    Google Earth Engine Geometry.
    """
    from shapely.geometry import mapping
    if ty == 'Polygon':
        return ee.Geometry.Polygon(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)
    elif ty == 'Point':
        return ee.Geometry.Point(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)

#%% Time Series Functions
def detrend(collection, return_detrend=False):
    """ 
    Simple linear detrender to center a dataset around mean zero with no linear trend.
    Used as a first step before fitting Harmonic Models.
    
    Original function adapted from:
    https://docs.google.com/document/d/1mNIRB90jwLuASO1JYas1kuOXCLbOoy1Z4NlV1qIXM10/edit
    
    Input is an ImageCollection, default output is a de-trended ImageCollection. Can return a 
    multi-band image with both the original and the detrended data as bands for comparison.
    
    """
    
    def addVariables(image):
        date = ee.Date(image.get('system:time_start'))
        years = date.difference(ee.Date('1970-01-01'), 'year');
        return image.addBands(ee.Image(years).rename('t')).float().addBands(ee.Image.constant(1))

    independents = ee.List(['constant', 't'])
    img = collection.first()
    bn = img.bandNames().get(0)
    dependent = ee.String(bn)
    
    coll = collection.map(addVariables)
    
    #Compute a linear trend.  This will have two bands: 'residuals' and 
    #a 2x1 band called coefficients (columns are for dependent variables).
    trend = coll.select(independents.add(dependent)).reduce(ee.Reducer.linearRegression(independents.length(), 1))
    
    #Flatten the coefficients into a 2-band image
    coefficients = trend.select('coefficients').arrayProject([0]).arrayFlatten([independents])
    
    #Compute a de-trended series.
    def remove_trend(image):
        detrended = image.select(dependent).subtract(image.select(independents).multiply(coefficients).reduce('sum'))\
            .rename('detrend').copyProperties(image, ['system:time_start'])
        return image.addBands(detrended)
    
    detrended = coll.map(remove_trend)
    
    if return_detrend:
        #This has both the original and detrended bands in there
        return detrended.select([dependent, 'detrend']) 
    else:
        #Default returns only the detrended data
        return detrended.select('detrend')

def Rolling_LinearFit(collection, windowSize=17):
    """
    Function to smooth a given ImageCollection by fitting piecewise linear regressions. Original code adapted from:
    https://stackoverflow.com/questions/47134534/how-to-smooth-ndvi-curve-in-google-earth-engine
    
    WindowSize is number of days on either side of date -- here the default is 17 to capture the previous and next 
    Landsat scene from a given date. This can be modified for any data type/smoothing length.
    
    Returns a smoothed ImageCollection with the same number of images/time sampling as the input dataset
    CAUTION: This does not do a good job of smoothing over gaps! Mileage may vary.
    
    """
    
    #First get the band name of the given Input ImageCollection
    img = collection.first()
    bn = img.bandNames().getInfo()[0]
    
    #Add a time band
    def add_time(img):
        return img.addBands(img.metadata('system:time_start').divide(1e18).rename('time'))
        
    #Create the windowed smoother
    def smoother(t):
        def applyFit(img):
            #Helper function to apply linear regression equation to each sub-collection
            out = img.select('time').multiply(fit.select('scale')).add(fit.select('offset'))\
                .set('system:time_start',img.get('system:time_start'))
            return out.rename(bn)
        
        t = ee.Date(t) #Get the date
        #Use the date to create a search window of windowSize in both directions
        window = data.filterDate(t.advance(-windowSize,'day'),t.advance(windowSize,'day'))
        
        #Perform a simple linear regression between time and the given band name
        fit = window.select(['time',bn]).reduce(ee.Reducer.linearFit())
    
        #Use the helper function to apply the piecewise linear fits to the original data
        #Forcing it to a list drops any empty images, the number of elements kept helps keep memory usage down
        out = window.map(applyFit).toList(5)
        return out

    #function to reduce time stacked linear regression results
    def reduceFits(t):
        t = ee.Date(t)
        return fitIC.filterDate(t.advance(-windowSize,'day'),t.advance(windowSize,'day')).mean()\
            .set('system:time_start',t.millis()).rename(bn)
    
    data = collection.map(add_time)
    dates = ee.List(data.aggregate_array('system:time_start'))
    fitIC = ee.ImageCollection(dates.map(smoother).flatten())
    smoothed = ee.ImageCollection(dates.map(reduceFits))
    return smoothed

#%% Harmonic Fitting Functions
def fit_harmonic(collection):
    """
    Function to fit a simple harmonic model to an ImageCollection. Uses only one harmonic with a frequency of 
    one calendar year. Original function adapted from:
    https://docs.google.com/document/d/1mNIRB90jwLuASO1JYas1kuOXCLbOoy1Z4NlV1qIXM10/edit
    
    Returns the harmonic-smoothed dataset, the phase of the harmonic, the amplitude of the harmonic
    and a simple R squared value for the fit. 

    """
    import numpy as np

    #Use these independent variables in the harmonic regression.
    harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin'])

    def addVariables(image):
        date = ee.Date(image.get('system:time_start'))
        #THE HARMONICS ARE FIT IN UNITS OF TIME (in this case, one year)
        years = date.difference(ee.Date('1970-01-01'), 'year')
        return image.addBands(ee.Image(years).rename('t')).float().addBands(ee.Image.constant(1))

    #Get the name of the image band
    img = collection.first()
    bn = img.bandNames().getInfo()[0]
    def add_harm_terms(image):
        #Add harmonic terms as new image bands.
        timeRadians = image.select('t').multiply(2 * np.pi)
        return image.addBands(timeRadians.cos().rename('cos')).addBands(timeRadians.sin().rename('sin'))
    
    #Add the new variables to the collection
    coll = collection.map(addVariables)
    harm_coll = coll.map(add_harm_terms)
    
    #The output of this reducer is a 4x1 array image.
    harmonicTrend = harm_coll.select(harmonicIndependents.add(bn))\
        .reduce(ee.Reducer.linearRegression(ee.Number(harmonicIndependents.length())))
    
    #Turn the array image into a multi-band image of coefficients.
    harmonicTrendCoefficients = harmonicTrend.select('coefficients').arrayProject([0]).arrayFlatten([harmonicIndependents])

    #Compute fitted values.
    def add_harm_fitted(image):
        return image.addBands(image.select(harmonicIndependents).multiply(harmonicTrendCoefficients).reduce('sum').rename('fitted'))
    fittedHarmonic = harm_coll.map(add_harm_fitted)
    
    #Compute distance between fitted and original
    def distance_to_orig(image):
        resids = image.select('fitted').subtract(image.select(bn))
        ss_res = resids.pow(ee.Number(2)).rename('SumSQ')
        dist_mean = image.select(bn).subtract(mn)
        ss_tot = dist_mean.pow(ee.Number(2)).rename('DistMean')
        #rsq = 1 - (ss_res / ss_tot)
        
        return image.addBands(ss_res).addBands(ss_tot)
    
    #Get the collection long-term mean
    mn = collection.reduce(ee.Reducer.mean())
    #Find the distance between the harmonic-smoothed data and the original data
    distfit = fittedHarmonic.map(distance_to_orig)
    
    #Sum of Resids
    sum_resids = distfit.select('SumSQ').reduce(ee.Reducer.sum())
    sum_sqtot = distfit.select('DistMean').reduce(ee.Reducer.sum())
    
    #Divide the summed images and subtract from 1 for a classic RSQ value
    rsq = ee.Image(1).subtract(sum_resids.divide(sum_sqtot)).rename('RSQ')
    
    def phase_from_harmonic(harmonicTrendCoefficients):
        #Compute phase and amplitude.
        phase = harmonicTrendCoefficients.select('sin').atan2(harmonicTrendCoefficients.select('cos')).rename('phase')
        amplitude = harmonicTrendCoefficients.select('sin').hypot(harmonicTrendCoefficients.select('cos')).rename('amplitude')
        
        return ee.Image(phase), ee.Image(amplitude)
    
    #Return the phase and amplitude values    
    phase, amplitude = phase_from_harmonic(harmonicTrendCoefficients)
    
    return [fittedHarmonic.select('fitted'), phase.select('phase').toFloat(), \
            amplitude.select('amplitude').toFloat(), rsq.select('RSQ').toFloat()]

def fit_multi_harmonic(collection, harmonics=3):
    """
    Function to fit a complex harmonic model to an ImageCollection. Uses any number of desired harmonics. 
    Original function adapted from:
    https://code.earthengine.google.com/2669122497313113fc4bb81bc8352828
    
    Returns the harmonic-smoothed dataset, the phase of the harmonic, the amplitude of the harmonic
    and a simple R squared value for the fit. 

    """
    import numpy as np

    #Get the name of the image band
    img = collection.first()
    bn = img.bandNames().getInfo()[0]

    #Get list of harmonic terms
    harmonicFrequencies = list(range(1, harmonics+1))#ee.List.sequence(1, harmonics)
    
    def getNames(base_bn, l):
        out = []
        for i in l:
            out.append(base_bn + str(i))
        return ee.List(out)
    
    #Construct lists of names for the harmonic terms.
    cosNames = getNames('cos_', harmonicFrequencies)
    sinNames = getNames('sin_', harmonicFrequencies)
    
    #Add the list of all sin/cos terms to the independent variables list
    independents = ee.List(['constant', 't']).cat(cosNames).cat(sinNames)
    
    def addVariables(image):
        date = ee.Date(image.get('system:time_start'))
        years = date.difference(ee.Date('1970-01-01'), 'year')
        #THE HARMONICS ARE FIT IN UNITS OF TIME (in this case, years)
        timeRadians = ee.Image(years.multiply(2 * np.pi))
        return image.addBands(timeRadians.rename('t').float()).addBands(ee.Image.constant(1))
        
    #Function to compute the specified number of harmonics and add them as bands.  
    #Assumes the time band is present.
    def addHarmonics(freqs):
        #This returns a function to be mapped over the ImageCollection, adding all required bands
        def fx(image):
            #Make an image of frequencies.
            frequencies = ee.Image.constant(freqs)
            #This band should represent time in radians.
            time = ee.Image(image).select('t')
            #Get the cosine terms.
            cosines = time.multiply(frequencies).cos().rename(cosNames)
            #Get the sine terms
            sines = time.multiply(frequencies).sin().rename(sinNames)
            return image.addBands(cosines).addBands(sines)
        return fx
    
    #Add the time and constant bands
    coll = collection.map(addVariables)
    #Add the arbitrary number of sin/cos bands
    harmonic_coll = coll.map(addHarmonics(harmonicFrequencies))
    
    #The output of the regression reduction is a 4x1 array image.
    harmonicTrend = harmonic_coll.select(independents.add(bn)).reduce(ee.Reducer.linearRegression(independents.length(), 1))
    
    #Turn the array image into a multi-band image of coefficients.
    harmonicTrendCoefficients = harmonicTrend.select('coefficients').arrayProject([0]).arrayFlatten([independents])
    #This has individual coefficents for each frequency (plus t and constant)

    #Compute fitted values.
    def add_harm_fitted(image):
        return image.addBands(image.select(independents).multiply(harmonicTrendCoefficients).reduce('sum').rename('fitted'))
    
    #This will return one fitted line -- the sum of all harmonic terms
    fittedHarmonic = harmonic_coll.map(add_harm_fitted)

    def phase_from_harmonic(harmonicTrendCoefficients, number):
        #Compute phase and amplitude.
        phase = harmonicTrendCoefficients.select('sin_' + str(number)).atan2(harmonicTrendCoefficients.select('cos_' + str(number))).rename('phase_' + str(number)) #.unitScale(-np.pi, np.pi)
        amplitude = harmonicTrendCoefficients.select('sin_' + str(number)).hypot(harmonicTrendCoefficients.select('cos_' + str(number))).rename('amplitude_' + str(number))
        
        return ee.Image(phase), ee.Image(amplitude)
    
    #To return all phase/amplitude values, loop over the number of sin/cos terms
    multiphase, multiamp = [], []
    for i in harmonicFrequencies:
        phase, amplitude = phase_from_harmonic(harmonicTrendCoefficients, i)
        multiphase.append(phase.toFloat())
        multiamp.append(amplitude.toFloat())
        
    #Convert these lists into multiband images for easier processing/exporting
    multiphase = ee.ImageCollection.fromImages(multiphase).toBands()
    multiamp = ee.ImageCollection.fromImages(multiamp).toBands()

    #Get Model fit
    def distance_to_orig(image):
        resids = image.select('fitted').subtract(image.select(bn))
        ss_res = resids.pow(ee.Number(2)).rename('SumSQ')
        dist_mean = image.select(bn).subtract(mn)
        ss_tot = dist_mean.pow(ee.Number(2)).rename('DistMean')
        #rsq = 1 - (ss_res / ss_tot)
        
        return image.addBands(ss_res).addBands(ss_tot)
    
    mn = collection.reduce(ee.Reducer.mean())
    distfit = fittedHarmonic.map(distance_to_orig)
    
    #Sum of Resids
    sum_resids = distfit.select('SumSQ').reduce(ee.Reducer.sum())
    sum_sqtot = distfit.select('DistMean').reduce(ee.Reducer.sum())
    
    #Divide the summed images and subtract from 1 for a classic RSQ value
    rsq = ee.Image(1).subtract(sum_resids.divide(sum_sqtot)).rename('RSQ')

    return [fittedHarmonic.select('fitted'), multiphase, multiamp, rsq]

#%% Data Aggregation Functions
def aggregate_to_yearly(collection, ds, de, agg_fx='sum'):
    '''
    Take an ImageCollection and convert it into a summed or average yearly value
    '''
    start, end = ee.Date(ds), ee.Date(de)
    #Generate list of years
    difdate = end.difference(start, 'year')
    length = ee.List.sequence(0, difdate.subtract(1))
    
    def gen_datelist(yr):
        return start.advance(yr, 'year')
    
    dates = length.map(gen_datelist)

    #Get band name
    bn = collection.first().bandNames().getInfo()[0]
    
    def reduceSum(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t, t.advance(1, 'year'))
        daysum = filt_coll.reduce(ee.Reducer.sum()).set('system:time_start', t.millis()).rename(bn)
        return daysum
    
    def reduceMean(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t, t.advance(1, 'month'))
        daymn = filt_coll.reduce(ee.Reducer.mean()).set('system:time_start', t.millis()).rename(bn)
        return daymn
    
    if agg_fx == 'sum':
        yr_agg = dates.map(reduceSum)
    elif agg_fx == 'mean':
        yr_agg = dates.map(reduceMean)
        
    #Convert back into an image collection
    yearly = ee.ImageCollection.fromImages(yr_agg)
    
    return yearly

def aggregate_to_monthly(collection, ds, de, agg_fx='sum'):
    '''
    Take an ImageCollection and convert it into a summed or average monthly value
    '''
    start, end = ee.Date(ds), ee.Date(de)
    #Generate length of months to look through
    difdate = end.difference(start, 'month')
    length = ee.List.sequence(0, difdate.subtract(1))
    
    def gen_datelist(mo):
        return start.advance(mo, 'month')
    
    dates = length.map(gen_datelist)

    #Get band name
    bn = collection.first().bandNames().getInfo()[0]
    
    def reduceSum(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t, t.advance(1, 'month'))
        daysum = filt_coll.reduce(ee.Reducer.sum()).set('system:time_start', t.millis()).rename(bn)
        return daysum
    
    def reduceMean(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t, t.advance(1, 'month'))
        daymn = filt_coll.reduce(ee.Reducer.mean()).set('system:time_start', t.millis()).rename(bn)
        return daymn
    
    #Map over the list of months, return either a mean or a sum of those values
    if agg_fx == 'sum':
        mo_agg = dates.map(reduceSum)
    elif agg_fx == 'mean':
        mo_agg = dates.map(reduceMean)
        
    #Convert back into an image collection
    monthly = ee.ImageCollection.fromImages(mo_agg)
    
    return monthly

def monthly_averages(collection, agg_fx='mean', years=None):
    '''
    Get the long-term average value for each month of the year over an ImageCollection
    Optionally select only certain years to use for the aggregation. 
    '''
    months = ee.List.sequence(1, 12)
    if agg_fx == 'mean':
        def fx(m):
            return collection.filter(ee.Filter.calendarRange(m, m, 'month')).mean().set('month', m)
    elif agg_fx == 'sum':
        def fx(m):
            return collection.filter(ee.Filter.calendarRange(m, m, 'month')).sum().set('month', m)
    
    if years:
        filtlist = []
        #Create a list of yearly filters
        for yr in years:
            #Create a set of filters to separate out the chosen years
            filt = ee.Filter.calendarRange(yr, yr, 'year')
            filtlist.append(filt)
            
        #Merge the filters together with Or
        multifilt = ee.Filter.Or(filtlist)
        #Filter the collection to only contain the chosen years
        collection = collection.filter(multifilt)
    
    monthly = ee.ImageCollection.fromImages(months.map(fx).flatten())
    return monthly

def create_anomalies(collection, ds, de):
    '''
    Calculate anomalies from long-term mean   
    '''
    coll_yearly = aggregate_to_yearly(collection, ds, de)
    lt_mn = coll_yearly.reduce(ee.Reducer.mean())
    
    def mapped_subtract(image):
        return image.subtract(lt_mn).set('system:time_start', image.get('system:time_start'))

    #Subtract the lont-term mean from each value
    anoms = coll_yearly.map(mapped_subtract)
    return anoms

#%% Data Import and Cleaning Functions
### GPM
def mask_GPM(image):
    mask = image.gt(0.1)
    precip_filt = image.updateMask(mask).set('system:time_start', image.get('system:time_start'))
    return precip_filt

### Landsat
def NDVI_L7(image):
    ndvi = image.normalizedDifference(['B4', 'B3']).rename('NDVI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def NDVI_L8(image):
    ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def NDWI_L8(image):
    ndwi = image.normalizedDifference(['B5', 'B6']).rename('NDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def NDWI_L7(image):
    ndwi = image.normalizedDifference(['B4', 'B5']).rename('NDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def maskL8sr(image):
    #Bits 3 and 5 are cloud shadow and cloud, respectively.
    cloudShadowBitMask = (1 << 3)
    cloudsBitMask = (1 << 5)
    #Get the pixel QA band.
    qa = image.select('pixel_qa')
    #Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))
    return image.updateMask(mask)

def cloudMaskL457(image):
    qa = image.select('pixel_qa')
    #If the cloud bit (5) is set and the cloud confidence (7) is high
    #or the cloud shadow bit is set (3), then it's a bad pixel.
    cloud = qa.bitwiseAnd(1 << 5).And(qa.bitwiseAnd(1 << 7)).Or(qa.bitwiseAnd(1 << 3))
    #Remove edge pixels that don't occur in all bands
    mask2 = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(cloud.Not()).updateMask(mask2)

### Sentinel-2
def NDVI_S2(image):
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def NDWI_S2(image):
    ndwi = image.normalizedDifference(['B8', 'B11']).rename('NDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def maskS2clouds(image):
    qa = image.select('QA60')
    
    #Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = (1 << 10)
    cirrusBitMask = (1 << 11)
    
    #Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    return image.updateMask(mask).divide(10000).set('system:time_start', image.get('system:time_start'))

#%% Conversion to Python Time Series
def export_to_pandas(collection, clipper, aggregation_scale, save=None):
    '''
    Takes an ImageCollection, an Earth Engine Geometry, and an aggregation scale (e.g., 30m for Landsat, 250m for MODIS, etc)
    
    Returns a pandas time series for the mean and standard deviation values over the 
    aggregation area. 
    
    Optionally saves those time series to a CSV file    
    
    '''
    import pandas as pd, numpy as np
    
    def createTS(image):
        date = image.get('system:time_start')
        value = image.reduceRegion(ee.Reducer.median(), clipper, aggregation_scale)
        std = image.reduceRegion(ee.Reducer.stdDev(), clipper, aggregation_scale)
        ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'Mean': value, 'STD': std})
        return ft
    
    TS = collection.map(createTS)
    dump = TS.getInfo()
    fts = dump['features']
    out_vals = np.empty((len(fts)))
    out_dates = []
    out_std = np.empty((len(fts)))
    
    for i, f in enumerate(fts):
        props = f['properties']
        date = props['date']
        try:
            val = list(props['Mean'].values())[0]
        except:
            val = np.nan
        try:
            std = list(props['STD'].values())[0]
        except:
            std = np.nan
        out_vals[i] = val
        out_std[i] = std
        out_dates.append(pd.Timestamp(date))
    
    ser = pd.Series(out_vals, index=out_dates)
    serstd = pd.Series(out_std, index=out_dates)
    if save:
        df = pd.DataFrame({'mean':out_vals, 'std':out_std, 'time':out_dates})
        df.to_csv(save + '.csv', index=False)
        print(save)
    return ser, serstd

#%% Short Example
minx2, maxx2 = 18, 20
miny2, maxy2 = -20, -18
from shapely.geometry import Polygon
geom = Polygon([[minx2, maxy2], [maxx2, maxy2], [maxx2, miny2], [minx2, miny2]])

search_area = gee_geometry_from_shapely(geom)

L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").filter(ee.Filter.date('2013-01-01', '2020-01-01'))\
    .filterBounds(search_area)\
    .filter(ee.Filter.lt('CLOUD_COVER', 15))\
    .map(maskL8sr)\
    .sort('system:time_start')

NDVI = L8.map(NDVI_L8).select('NDVI')

harmonic, phase, amplitude, rsq = fit_multi_harmonic(detrend(NDVI), 5)
run_export(phase, 'epsg:4326', 'Landsat8_HarmonicPhase', scale=30, region=search_area, maxPixels=1e12)
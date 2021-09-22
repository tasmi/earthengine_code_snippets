# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 11:16:43 2020

@author: tsmith
"""

import ee
ee.Initialize()
import numpy as np

#%% General Helper Functions
def run_export(image, crs, filename, scale, region, folder=None, maxPixels=1e12, cloud_optimized=True):
    '''
    Runs an export function on GEE servers
    '''
    task_config = {'fileNamePrefix': filename,'crs': crs,'scale': scale,'maxPixels': maxPixels, 'fileFormat': 'GeoTIFF', 'formatOptions': {'cloudOptimized': cloud_optimized}, 'region': region,}
    if folder:
        task_config['folder'] = folder
    task = ee.batch.Export.image.toDrive(image, filename, **task_config)
    task.start()
    
def export_collection(collection, region, prefix, crs=None, scale=100, start_image=0, max_images=None, folder=None):
    '''
    Exports all images within an image collection for a given region. All files named by a prefix (given)
    and their image date (formated YYYYMMDD). 
    region: area to export
    prefix: file name prefix
    crs: can be provided, or determined automatically
    scale: output image pixel size in meters
    start_image: where to start in the list of images (e.g., if you need to break your job up into pieces)
    max_images: number of images to export (e.g., maximum)
    folder: if you want to store all images in a separate folder in your GDrive
    '''
    if not crs:
        crs = collection.first().projection()
    
    nr_images = int(collection.size().getInfo())    
    image_list = collection.toList(nr_images)
    
    if max_images:
        if max_images < nr_images:
            nr_images = start_image + max_images #Make sure not to export too many if you want to test something
        else:
            #If the number of images to export is less than the max_images, pass
            pass
        
    print('Exporting up to %i Images' % nr_images)
    
    #Run a list from the starting image to the number you want
    for i in range(start_image, nr_images):
        if i >= start_image:
            image = ee.Image(image_list.get(i))
            try:
                #If there are defined start and end dates, add them to the file names
                ds = image.get('sdate')
                de = image.get('edate')
                date_name0 = ee.Date(ds).format('YYYYMMdd').getInfo()
                date_name1 = ee.Date(de).format('YYYYMMdd').getInfo()
                date_name = date_name0 + '-' + date_name1
            except:
                date = image.get('system:time_start')
                date_name = ee.Date(date).format('YYYYMMdd').getInfo()
            output_name = prefix + '_' + date_name + '_' + str(scale) + 'm.tif'
            run_export(image, crs=crs, filename=output_name, scale=scale, region=region, folder=folder)
            print('Started export for image ' + str(i) + '(' + date_name + ')')
            
def build_collection_from_search(ic, searchranges, var, agg_fx):
    '''
    Given a set of date ranges (e.g., [start_date, end_date]), create an image collection that is 
    one image per date range. Each image will have start and end dates as additional metadata.
    The image will be single band if 'var' is one band, or multiband if var is of the form ['band1', 'band2', etc]
    '''
    def get_red(agg_fx):
        '''Get a reducer from a dictionary'''
        agg_dict = {'mn':ee.Reducer.mean(), 'max':ee.Reducer.max(), 'min':ee.Reducer.min(),\
                   'sum':ee.Reducer.sum(), 'median':ee.Reducer.median()}
        return agg_dict[agg_fx]
        
    merged_images = []
    for i in searchranges:
        ds, de = i
        subcol = ic.filterDate(ds, de)
        d = ee.Date(ds)
        e = ee.Date(de)
        agg = subcol.select(var).reduce(get_red(agg_fx))
        agg = ee.Image(agg).set('system:time_start', d)
        agg = ee.Image(agg).set('sdate', d)
        agg = ee.Image(agg).set('edate', e)
        merged_images.append(agg)
    
    #Turn them into a single imagecollection
    output = ee.ImageCollection.fromImages(merged_images)
    return output

def gee_geometry_from_shapely(geom, crs='epsg:4326'):
    """ 
    Simple helper function to take a shapely geometry and a coordinate system and convert them to a 
    Google Earth Engine Geometry.
    """
    from shapely.geometry import mapping
    ty = geom.type
    if ty == 'Polygon':
        return ee.Geometry.Polygon(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)
    elif ty == 'Point':
        return ee.Geometry.Point(ee.List(mapping(geom)['coordinates']), proj=crs)    
    elif ty == 'MultiPolygon':
        return ee.Geometry.MultiPolygon(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)
    elif ty == 'MultiPoint':
        return ee.Geometry.MultiPoint(ee.List(mapping(geom)['coordinates']), proj=crs)
    
def customRemap(image, lowerLimit, upperLimit, newValue):
    mask = image.gt(lowerLimit).And(image.lte(upperLimit))
    return image.where(mask, newValue)

def remove_null(collection):
    '''
    Remove empty images from a collection
    '''
    def flag_null(image):
        return image.set('band_count', image.bandNames().length())
    return collection.map(flag_null).filter(ee.Filter.gt('band_count', 0))

def mask_invalid(collection, minval, maxval, band=None):
    '''
    Mask all images in a collection by some min and max value
    '''
    
    if band:
        collection = collection.select(band)
    
    def apply_mask(image):
        mask1 = image.lt(maxval)
        mask2 = image.gt(minval)
        return image.updateMask(mask1).updateMask(mask2)
    return collection.map(apply_mask)

def apply_mask(collection, mask):
    '''
    Simple function to apply a static mask to all images in a collection
    '''
    
    def apply_fx(image):
        return image.updateMask(mask)
    return collection.map(apply_fx)

def add_doy(image):
    ''' Add a day of year as an image band '''
    #ts = image.get('system:time_start')
    #return image.addBands(ee.Image.constant(ee.Number.parse(image.date().millis())).rename('day').float())
    #return image.addBands(ee.Image.constant(ee.Number.parse(image.date().format("YYYYMMdd"))).rename('day').float())
    doy = ee.Image.constant(ee.Number.parse(image.date().format("D"))).rename('day').toUint32()#.float().set('system:time_start', ts)
    return image.addBands(doy)

def add_date(image):
    ''' Add a date as an image band '''
    date = ee.Image.constant(ee.Number.parse(image.date().format('YYYYMMdd'))).rename('date').toUint32()
    return image.addBands(date)

def invert(image):
    ''' Invert an Image '''
    return image.multiply(-1)

def get_ic_dates(ic):
    '''
    Get all of the dates in a given ImageCollection. Useful for finding unique sets
    '''
    def ret_date(image):
        date = ee.Image(image).get('system:time_start')
        d = ee.Date(date).format('Y/M/d')
        return image.set('date', d)
    
    ic_d = ic.map(ret_date)
    datelist = ee.List(ic_d.aggregate_array('date'))
    return datelist

#%% Time Series Functions
def prevdif(collection):
    ''' 
    Get the simple time difference between sequential images in an image collection.
    NOTE: It is important to do a tight spatial filtering first!
    '''
    def dif(f):
        #sdate = ee.Date(ee.Image(ee.List(f).get(0)).get('system:time_start'))
        f = ee.Number(f)
        simage = ee.Image(ic_list.get(f))
        snext = ee.Image(ic_list.get(f.add(1)))
        sdate = simage.get('system:time_start')
        
        d = snext.subtract(simage).set('system:time_start', sdate)
        
        return d
    
    ic_list = collection.sort('system:time_start').toList(collection.size())
    seq = ee.List.sequence(0, collection.size().subtract(2)) #Drop the final image since we don't have the next to subtract
    ic_diff = seq.map(dif)
    
    return ee.ImageCollection.fromImages(ic_diff)

def windowed_difference(collection, window):
    ''' 
    Get a windowed difference between sequential images in an image collection.
    NOTE: It is important to do a tight spatial filtering first!
    '''
    
    half_window = int(window / 2)
    
    def windif_fx(f):
        f = ee.Number(f)
        edate = ee.Image(image_list.get(f)).get('system:time_start')
        
        forward_slice = image_list.slice(f, f.add(half_window))
        backward_slice = image_list.slice(f.subtract(half_window), f)
        
        forward_mean = ee.ImageCollection.fromImages(forward_slice).reduce(ee.Reducer.mean())
        backward_mean = ee.ImageCollection.fromImages(backward_slice).reduce(ee.Reducer.mean())
        
        return ee.Image(forward_mean.subtract(backward_mean)).set('system:time_start', edate)
            
    ln = ee.Number(collection.size()) #Length of time series
    image_list = collection.sort('system:time_start').toList(ln)
    
    seq = ee.List.sequence(half_window, ln.subtract(half_window)) #Drop the final image since we don't have the next to subtract
    windif = seq.map(windif_fx)
    
    return ee.ImageCollection.fromImages(windif)
    
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
        years = date.difference(ee.Date('1970-01-01'), 'year')
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
        #Forcing it to a list drops any empty images
        out = window.map(applyFit).toList(data.size())
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
def match_time_sampling(collection1, collection2, ds, de, reducer):
    ''' 
    Take an IC and aggregate to the time sampling of another IC. 
    Referenced to the first collection.
    ds, de = start and end date
    '''
    
    def make_filter(item):
        image = ee.Image(item)
        sd = image.get('system:time_start')
        ed = image.get('system:time_end')
        filt = ee.Filter.date(sd, ed)
        return filt
    
    c1 = collection1.toList(collection1.size())
    filts = c1.map(make_filter)
    
    def create_sub_collections(filt):
        filt_coll = collection2.filter(filt)
        return filt_coll.set('bandcount', ee.Number(filt_coll.size()))
    
    mc = filts.map(create_sub_collections)
    mc_filt = mc.filter(ee.Filter.gt('bandcount',0))
    
    def map_red(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        filt_date = filt_coll.first().get('system:time_start')
        return filt_coll.reduce(reducer).set('system:time_start', filt_date)
    
    agged = mc_filt.map(map_red)
    
    return ee.ImageCollection(agged)  

def aggregate_to(collection, ds, de, timeframe='month', skip=1, agg_fx='sum', agg_var=None):
    '''
    Take an ImageCollection and convert it into an aggregated value based on an arbitrary function.
    Several useful functions are included with keywords (mean, sum, etc), but an arbitrary reducer can be supplied
    
    day, month, year are the typical arguments
    
    skip will allow you to make larger windows (e.g., 5 days)
    
    '''
    start, end = ee.Date(ds), ee.Date(de)
    #Generate length of months to look through
    difdate = end.difference(start, timeframe)
    length = ee.List.sequence(0, difdate.subtract(1))
    
    if not skip == 1:
        length_py = length.getInfo()
        length_skip = length_py[::skip]
        length = ee.List(length_skip)
    
    def gen_datelist(t):
        return start.advance(t, timeframe)
    
    dates = length.map(gen_datelist)
    
    def create_sub_collections(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t, t.advance(skip, timeframe)) #Move forward the 'skip' amount of time units
        return filt_coll.set('bandcount', ee.Number(filt_coll.size())).set('system:time_start', t.millis())
    
    mc = dates.map(create_sub_collections)
    mc_filt = mc.filter(ee.Filter.gt('bandcount',0))
    
    #Get band name -- this covers for possible empty parts of the collection
    bn = ee.Image(collection.reduce(ee.Reducer.mean())).bandNames().getInfo()
    try:
        bn = bn[0]
    except:
        bn = 'empty'
    #print(bn)

    def reduceSum(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daysum = filt_coll.reduce(ee.Reducer.sum()).set('system:time_start', t).rename(bn)
        return daysum
    
    def reduceMean(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.mean()).set('system:time_start', t).rename(bn)
        return daymn
    
    def reduceMin(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.min()).set('system:time_start', t).rename(bn)
        return daymn
    
    def reduceMax(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.max()).set('system:time_start', t).rename(bn)
        return daymn

    def reduceMedian(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.median()).set('system:time_start', t).rename(bn)
        return daymn
    
    def reduceSTD(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.stdDev()).set('system:time_start', t).rename(bn)
        return daymn
    
    def reduceIQR(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        pcts = filt_coll.reduce(ee.Reducer.percentile([25,75]))
        iqr = pcts.select(bn + '_p75').subtract(pcts.select(bn + '_p25')).toFloat().set('system:time_start', t).rename(bn)
        return iqr
    
    def reduce9010(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)  
        t = filt_coll.get('system:time_start')
        pcts = filt_coll.reduce(ee.Reducer.percentile([10,90]))
        iqr = pcts.select(bn + '_p90').subtract(pcts.select(bn + '_p10')).toFloat().set('system:time_start', t).rename(bn)
        return iqr
    
    def reduce955(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)  
        pcts = filt_coll.reduce(ee.Reducer.percentile([5,95]))
        t = filt_coll.get('system:time_start')
        iqr = pcts.select(bn + '_p95').subtract(pcts.select(bn + '_p5')).toFloat().set('system:time_start', t).rename(bn)
        return iqr
    
    def binary_mask(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        if agg_var == 'NDVI':
            def reclass(image):
                #No veg == 0, veg == 1
                remapped = customRemap(image, -1, 0.4, 0)
                remapped2 = customRemap(remapped, 0.4, 1, 1)
                return remapped2
        if agg_var == 'NDVI2':
            def reclass(image):
                #No veg == 0, veg == 1
                remapped = customRemap(image, -1, 0.5, 0)
                remapped2 = customRemap(remapped, 0.5, 1, 1)
                return remapped2
        if agg_var == 'NDWI':
            #Map all positive values to zero, all negative values to 1
            def reclass(image):
                remapped = customRemap(image, -1, 0.45, 0)
                remapped2 = customRemap(remapped, 0.45, 1, 1)
                return remapped2
        if agg_var == 'MNDWI':
            #Map all positive values to zero, all negative values to 1
            def reclass(image):
                remapped = customRemap(image, -1, 0.5, 0)
                remapped2 = customRemap(remapped, 0.5, 1, 1)
                return remapped2
        if agg_var == 'S1':
            def reclass(image):
                remapped = customRemap(image, -25, 0, 0)
                remapped2 = customRemap(remapped, -60, -25, 1)
                return remapped2

        rm = filt_coll.map(reclass)
        occur = rm.reduce(ee.Reducer.sum()).toUint8().set('system:time_start', t).rename(bn)
        return occur
    
    #Map over the list of months, return either a mean or a sum of those values
    if agg_fx == 'sum':
        mo_agg = mc_filt.map(reduceSum)
    elif agg_fx == 'min':
        mo_agg = mc_filt.map(reduceMin)
    elif agg_fx == 'max':
        mo_agg = mc_filt.map(reduceMax)
    elif agg_fx == 'mean':
        mo_agg = mc_filt.map(reduceMean)
    elif agg_fx == 'median':
        mo_agg = mc_filt.map(reduceMedian)
    elif agg_fx == 'std':
        mo_agg = mc_filt.map(reduceSTD)
    elif agg_fx == 'binary_mask':
        mo_agg = mc_filt.map(binary_mask)
    elif agg_fx == 'iqr':
        mo_agg = mc_filt.map(reduceIQR)
    elif agg_fx == '9010':
        mo_agg = mc_filt.map(reduce9010)
    else:
        mo_agg = mc_filt.map(agg_fx)
        
    #Convert back into an image collection
    agged = ee.ImageCollection.fromImages(mo_agg)
    agged = agged.filter(ee.Filter.listContains('system:band_names', bn))
    
    return agged

def windowed_difference_date(collection, ds, de, window_size, window_unit):
    ''' 
    Get a windowed difference between sequential images in an image collection.
    NOTE: It is important to do a tight spatial filtering first!
    '''
    
    start, end = ee.Date(ds), ee.Date(de)
    difdate = end.difference(start, window_unit)
    ln = ee.List.sequence(0, difdate.subtract(1))
    
    half_window = int(window_size / 2)
        
    def gen_datelist(t):
        return start.advance(t, window_unit)
    
    dates = ln.map(gen_datelist)
    
    def create_sub_collections(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t.advance(-1*half_window, window_unit), t.advance(half_window, window_unit))
        return filt_coll.set('bandcount', ee.Number(filt_coll.size())).set('date', t)
    
    mc = dates.map(create_sub_collections)
    mc_filt = mc.filter(ee.Filter.gt('bandcount',0))
    
    def windif_fx(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        edate = ee.Date(filt_coll.get('date'))
        
        forward_slice = filt_coll.filterDate(edate.advance(-1*half_window, window_unit), edate)
        backward_slice = filt_coll.filterDate(edate, edate.advance(half_window, window_unit))
        
        #forward_slice = image_list.slice(f, f.add(half_window))
        #backward_slice = image_list.slice(f.subtract(half_window), f)
        
        forward_mean = forward_slice.reduce(ee.Reducer.mean())
        backward_mean = backward_slice.reduce(ee.Reducer.mean())
        
        return ee.Image(forward_mean.subtract(backward_mean)).set('system:time_start', edate)
            
    windif = mc_filt.map(windif_fx)
    
    return ee.ImageCollection(windif)

def reduceFit(collection):
    t = collection.get('system:time_start')
    bn = collection.first().bandNames().getInfo()[0]
    def createTimeBand(image):
        date = ee.Date(image.get('system:time_start'))
        years = date.difference(ee.Date('1970-01-01'), 'month')
        return image.addBands(ee.Image(years).rename('t')).float().addBands(ee.Image.constant(1))
    
    c = collection.map(createTimeBand)
    fit = c.select(['t',bn]).reduce(ee.Reducer.linearFit()).select('scale')
    return fit.set('system:time_start', t)

def windowed_monthly_agg(collection, ds, de, window=1, agg_fx='sum'):
    '''
    Centered window aggregation (pm window size)
    '''
    start, end = ee.Date(ds), ee.Date(de)
    difdate = end.difference(start, 'month')
    length = ee.List.sequence(0, difdate.subtract(1))
    
    def gen_datelist(mo):
        return start.advance(mo, 'month')
    
    dates = length.map(gen_datelist)

    bn = collection.first().bandNames().getInfo()[0]
    
    def create_sub_collections(t):
        t = ee.Date(t)
        filt_coll = collection.filterDate(t.advance(-1*window, 'month'), t.advance(window, 'month'))
        return filt_coll.set('bandcount', ee.Number(filt_coll.size())).set('system:time_start', t.millis())
    
    mc = dates.map(create_sub_collections)
    mc_filt = mc.filter(ee.Filter.gt('bandcount',0))
    
    def reduceFit(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        def createTimeBand(image):
            date = ee.Date(image.get('system:time_start'))
            years = date.difference(ee.Date('1970-01-01'), 'month')
            return image.addBands(ee.Image(years).rename('t')).float().addBands(ee.Image.constant(1))
        
        c = filt_coll.map(createTimeBand)
        fit = c.select(['t',bn]).reduce(ee.Reducer.linearFit()).select('scale')
        return fit.set('system:time_start', t)
    
    def reduceSum(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daysum = filt_coll.reduce(ee.Reducer.sum()).set('system:time_start', t).rename(bn)
        return daysum
    
    def reduceMean(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.mean()).set('system:time_start', t).rename(bn)
        return daymn

    def reduceMedian(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.median()).set('system:time_start', t).rename(bn)
        return daymn
    
    def reduceSTD(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        daymn = filt_coll.reduce(ee.Reducer.stdDev()).set('system:time_start', t).rename(bn)
        return daymn
        
    def reduceIQR(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        pcts = filt_coll.reduce(ee.Reducer.percentile([25,75]))
        iqr = pcts.select(bn + '_p75').subtract(pcts.select(bn + '_p25')).toFloat().set('system:time_start', t).rename(bn)
        return iqr
    
    def reduce9010(filt_coll):
        filt_coll = ee.ImageCollection(filt_coll)
        t = filt_coll.get('system:time_start')
        pcts = filt_coll.reduce(ee.Reducer.percentile([10,90]))
        iqr = pcts.select(bn + '_p90').subtract(pcts.select(bn + '_p10')).toFloat().set('system:time_start', t).rename(bn)
        return iqr
        
    #Map over the list of months, return either a mean or a sum of those values
    if agg_fx == 'sum':
        mo_agg = mc_filt.map(reduceSum)
    elif agg_fx == 'mean':
        mo_agg = mc_filt.map(reduceMean)
    elif agg_fx == 'median':
        mo_agg = mc_filt.map(reduceMedian)
    elif agg_fx == 'std':
        mo_agg = mc_filt.map(reduceSTD)
    elif agg_fx == 'fit':
        mo_agg = mc_filt.map(reduceFit)
    elif agg_fx == 'iqr':
        mo_agg = mc_filt.map(reduceIQR)
    elif agg_fx == '9010':
        mo_agg = mc_filt.map(reduce9010)
        
    monthly = ee.ImageCollection.fromImages(mo_agg)
    #monthly = monthly.filter(ee.Filter.listContains('system:band_names', bn))
    
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

def seasonal_composite(monthly, season):
    '''
    Break out seasonal slices of a collection
    '''
    allfilts = []
    if season == 'DJF':
        for m in [12, 1, 2]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
        filt = ee.Filter.Or(allfilts)
    if season == 'MAM':
        for m in [3, 4, 5]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
        filt = ee.Filter.Or(allfilts)
    if season == 'JJA':
        for m in [6, 7, 8]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
        filt = ee.Filter.Or(allfilts)
    if season == 'SON':
        for m in [9, 10, 11]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
        filt = ee.Filter.Or(allfilts)    
    return monthly.filter(filt)

def to_percent(collection, threshold_min, threshold_max):
    '''
    Finds the amount of time (percentage) a collection is within a given threshold.
    '''
    def newbounds(image):
        #remap here to max of one
        rc = ee.Image(0).where(image.gt(threshold_min).And(image.lte(threshold_max)), 1)
        return rc
    
    flt = collection.map(newbounds)
    su = flt.reduce(ee.Reducer.sum()).toFloat()
    ct = flt.reduce(ee.Reducer.count()).toFloat()
    pct = su.divide(ct).multiply(100).rename('PCT')
    return pct

#%% Join Collections
def join_timeagg_collections(c1, c2):
    filt = ee.Filter.equals(leftField='system:time_start', rightField='system:time_start') 
    #filt = ee.Filter.equals(leftField='system:index', rightField='system:index') 
    innerJoin = ee.Join.inner() #initialize the join
    innerJoined = innerJoin.apply(c1, c2, filt) #This is a FEATURE COLLECTION
    def combine_joined(feature):
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
    
    joined_collect = ee.ImageCollection(innerJoined.map(combine_joined))
    return joined_collect

def two_band_reg(c1, c2, crs, name, polygon, scale=30):
    #https://developers.google.com/earth-engine/joins_simple
    filt = ee.Filter.equals(leftField='system:time_start', rightField='system:time_start') 
    #filt = ee.Filter.equals(leftField='system:index', rightField='system:index')     
    innerJoin = ee.Join.inner() #initialize the join
    innerJoined = innerJoin.apply(c1, c2, filt) #This is a FEATURE COLLECTION
    def combine_joined(feature):
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
    
    joined_collect = ee.ImageCollection(innerJoined.map(combine_joined))
    
    #Get the band names
    bn1 = c1.first().bandNames().get(0)
    bn2 = c2.first().bandNames().get(0)
        
    #Add a constant band
    def createConstantBand(image):
        return ee.Image(1).addBands(image)
    
    prepped = joined_collect.map(createConstantBand)
    
    #Filter to make sure only images with both bands are used
    #f1 = ee.Filter.listContains('bands', bn1)
    #f2 = ee.Filter.listContains('bands', bn2)
    #filt = prepped.filter(ee.Filter.And([f1, f2]))
    #filt = prepped.filterMetadata('bandNames','contains',bn1)
    
    #filt = ee.ImageCollection.fromImages(prepped.toList(prepped.size()))
        
    #Now using that joined collection, do a regression
    #fit = prepped.select(['constant', bn1, bn2]).reduce(ee.Reducer.linearRegression(numX=2, numY=1))
    fit = prepped.select(['constant', bn1, bn2]).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1))
    #fit = filt.select(['constant', bn1, bn2]).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1))
    lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['constant', 'trend']]).select('trend')
    run_export(lrImage, crs, name + '_RobustLinReg', scale, polygon)

def same_inst_twoband_reg(collection, crs, name, polygon, scale=30, export='Slope', rmse=None, constant=True):
    '''
    Export can be 'Slope', 'Intercept', or 'Both'
    '''
    bn = collection.first().bandNames()
    
    def createConstantBand(image):
        return ee.Image(1).addBands(image)
    
    if constant:
        prepped = collection.map(createConstantBand)
        
        var = ee.List(['constant']).cat(bn)
        
        fit = prepped.select(var).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1))
        lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['constant', 'trend']])
    else:
        fit = collection.select(bn).reduce(ee.Reducer.robustLinearRegression(numX=1, numY=1))
        lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['trend']])
    if export in ['Slope', 'Both']:
        run_export(lrImage.select('trend'), crs, name + '_RobustLinReg_slope', scale, polygon)
    if export in ['Intercept', 'Both']:
        run_export(lrImage.select('constant'), crs, name + '_RobustLinReg_intercept', scale, polygon)
        
    #Get the RMSE
    if rmse:
        lrImage = fit.select(['residuals']).arrayFlatten([['residuals']])
        run_export(lrImage, crs, name + '_RobustLinReg_rmse', scale, polygon)
        
def bootstrap_slope(collection, crs, name, polygon, scale=30, export='Slope', constant=True, n_iter=100, sample_size=0.75):
    '''
    Export can be 'Slope', 'Intercept', or 'Both'. Does a set n_iter using sample_size percent of the 
    available DATES within a collection. This will work only if neighboring scenes aren't the same
    date (e.g., you have filtered by path/row already). Otherwise the sampling percentages will be 
    off due to multiple geographic regions having the same date!
    '''
    import pandas as pd
    bn = collection.first().bandNames()
    
    def createConstantBand(image):
        return ee.Image(1).addBands(image).set('system:time_start', image.get('system:time_start'))
    
    prepped = collection.map(createConstantBand)
    var = ee.List(['constant']).cat(bn)
    
    def subsample(ic, iter_nr):
        def fmt(d):
            #Reformat a date to EarthEngine standard
            return d.strftime('%Y-%m-%d')
            
        #Get unique dates to use for the random sampling
        UD = get_ic_dates(ic).getInfo()
        UDS = pd.Series(range(len(UD)), index=UD)
        ud = UDS.index.unique()
    
        #Get random selection of unique dates
        r_sample = pd.Series(ud).sample(frac=sample_size,random_state=iter_nr)
        
        #Go through each date and add it to a list of filters
        flist = []
        for date in r_sample.values:
            date = pd.Timestamp(date)
            sd, ed = fmt(date), fmt(date + pd.Timedelta('1 day'))
            f = ee.Filter.date(sd, ed)
            flist.append(f)
            
        #Merge all the filters with an 'or' to make the random subset
        filt = ic.filter(ee.Filter.Or(flist))
        
        return filt
    
    def run_reg(ic):   
        if constant:
            fit = ic.select(var).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1))
            lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['constant', 'trend']])
        else:
            fit = ic.select(bn).reduce(ee.Reducer.robustLinearRegression(numX=1, numY=1))
            lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['trend']])
            
        return lrImage
    
    slopes, intercepts = [], []
    for i in range(n_iter):
        filt = subsample(prepped, i)
        lrImage = run_reg(filt)
        slopes.append(lrImage.select('trend'))
        intercepts.append(lrImage.select('constant'))
        
    slope_collection = ee.ImageCollection.fromImages(slopes)
    intercept_collection = ee.ImageCollection.fromImages(intercepts)
    
    #Now get some statistics
    def return_stats(ic):
        std = ic.reduce(ee.Reducer.stdDev())
        mn = ic.reduce(ee.Reducer.mean())
        pctiles = []
        for p in plist:
            pct = ic.reduce(ee.Reducer.percentile([p]))
            pctiles.append(pct)
        return std, mn, pctiles
    
    plist = [5, 25, 50, 75, 95]
    outname = name + '-n_iter' + str(n_iter) + '-pct' + str(sample_size).replace('.','')

    if export in ['Slope', 'Both']:
        s, m, p = return_stats(slope_collection)
        run_export(s, crs, outname + 'slopeSTD', scale, polygon)
        run_export(m, crs, outname + 'slopeMN', scale, polygon)
        for i, pct in enumerate(plist):
            run_export(p[i], crs, outname + 'slopePCT_p' + str(pct), scale, polygon)
    if export in ['Intercept', 'Both']:
        s, m, p = return_stats(intercept_collection)
        run_export(s, crs, outname + 'intSTD', scale, polygon)
        run_export(m, crs, outname + 'intMN', scale, polygon)
        for i, pct in enumerate(plist):
            run_export(p[i], crs, outname + 'intPCT_p' + str(pct), scale, polygon)

#%% Data Import and Cleaning Functions
### GPM
def mask_GPM(image):
    mask = image.gt(0.1)
    precip_filt = image.updateMask(mask).set('system:time_start', image.get('system:time_start'))
    return precip_filt

### Landsat
def NDSI_L7(image):
    ndvi = image.normalizedDifference(['B2', 'B5']).rename('NDSI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def NDSI_L8(image):
    ndvi = image.normalizedDifference(['B3', 'B6']).rename('NDSI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

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

def MNDWI_L8(image):
    ndwi = image.normalizedDifference(['B3', 'B6']).rename('MNDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def MNDWI_L7(image):
    ndwi = image.normalizedDifference(['B2', 'B5']).rename('MNDWI').set('system:time_start', image.get('system:time_start'))
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

def maskL8toa_simple(image):
    cloudsBitMask = (1 << 4)
    #Get the pixel QA band.
    qa = image.select('BQA')
    #Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudsBitMask).eq(0)
    return image.updateMask(mask)

def maskL8toa(image):
    #via: https://gis.stackexchange.com/questions/274612/apply-a-cloud-mask-to-a-landsat8-collection-in-google-earth-engine-time-series
    #via: https://gis.stackexchange.com/questions/292835/using-cloud-confidence-to-create-cloud-mask-from-landsat-8-bqa
    def extractQABits(qaBand, bitStart, bitEnd):
        numBits = bitEnd - bitStart + 1
        qaBits = qaBand.rightShift(bitStart).mod(2**numBits)
        return qaBits
    
    qa_band = image.select('BQA')
    cloud_shadow = extractQABits(qa_band, 7, 8)
    mask_shadow = cloud_shadow.gte(2)
    cloud = extractQABits(qa_band, 5, 6)
    mask_cloud = cloud.gte(2)
    
    combined_mask = (cloud.Or(cloud_shadow))#.Not()
    return image.updateMask(combined_mask)

def cloudMaskL457(image):
    qa = image.select('pixel_qa')
    #If the cloud bit (5) is set and the cloud confidence (7) is high
    #or the cloud shadow bit is set (3), then it's a bad pixel.
    cloud = qa.bitwiseAnd(1 << 5).And(qa.bitwiseAnd(1 << 7)).Or(qa.bitwiseAnd(1 << 3))
    #Remove edge pixels that don't occur in all bands
    mask2 = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(cloud.Not()).updateMask(mask2)

### Sentinel-2
def MSAVI_S2(image):
    msavi2 = image.select('B8').multiply(2).add(1).subtract(image.select('B8').multiply(2).add(1).pow(2)\
                .subtract(image.select('B8').subtract(image.select('B4')).multiply(8)).sqrt())\
                .divide(2)\
                .set('system:time_start', image.get('system:time_start'))\
                .rename('MSAVI2')
    return image.addBands(msavi2)

def NDVI_S2(image):
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def NDWI_S2(image):
    ndwi = image.normalizedDifference(['B8', 'B11']).rename('NDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def MNDWI_S2(image):
    ndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndwi)

def LAI_S2(image):
    #(5.405 * ((B8A - B5) / (B8A + B5))) - 0.114
    calc = image.expression('(5.405 * ((B8A - B5) / (B8A + B5))) - 0.114', {'B8A': image.select('B8A'),'B5': image.select('B5')})
    lai = calc.rename('LAI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(lai)

def maskS2clouds(image):
    qa = image.select('QA60')
    
    #Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = (1 << 10)
    cirrusBitMask = (1 << 11)
    
    #Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    return image.updateMask(mask).divide(10000).set('system:time_start', image.get('system:time_start'))

def rescale_modis(image):
    return image.multiply(0.0001).set('system:time_start', image.get('system:time_start'))

def L8_Temp(image):
    #Get Fractional Veg
    ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI').select('NDVI')
    #mm = ndvi.reduceRegion(ee.Reducer.minMax(), polygon, 100)
    #nmin = ee.Number(mm.get('NDVI_min'))
    #nmax = ee.Number(mm.get('NDVI_max'))
    nmin, nmax = ee.Number(0.2), ee.Number(0.5)
    fv = (ndvi.subtract(nmin).divide(nmax.subtract(nmin))).pow(ee.Number(2)).rename('FV')
    
    #Emissivity
    a = ee.Number(0.004)
    b = ee.Number(0.986)
    EM = fv.multiply(a).add(b).rename('EMM')
    
    #Calc of LST
    thermal = image.select('B10').multiply(0.1)
    LST = thermal.expression('(Tb/(1 + (0.00115* (Tb / 1.438))*log(Ep)))',\
                             {'Tb': thermal.select('B10'), \
                              'Ep': EM.select('EMM')}).rename('LST')
    return LST.select('LST').set('system:time_start', image.get('system:time_start')) #-273.15

def focal_med_filt(collection, radius=100):
    ''' 
    Apply a focal median filter to a selected band, with flexible radius
    '''
    bn = collection.first().bandNames().getInfo()
    
    def applyfx(image):
        for b in bn:
            sel = image.select(b)
            smoothed = sel.focal_median(radius, 'circle', 'meters')
            image = image.addBands(smoothed.rename(b + '_filt'))
        return image
    return collection.map(applyfx)
    
def focal_range_filt(collection, radius=1.5, radiustype='pixels', kernel='square'):
    ''' 
    Apply a focal range (max - min) filter to a selected band, with flexible radius
    '''
    bn = collection.first().bandNames().getInfo()
    
    def applyfx(image):
        for b in bn:
            sel = image.select(b)
            ma = sel.focal_max(radius=radius, kernelType=kernel, units=radiustype)
            mi = sel.focal_min(radius=radius, kernelType=kernel, units=radiustype)
            dif = ma.subtract(mi)
            image = image.addBands(dif.rename(b + '_range'))
        return image
    return collection.map(applyfx)

#%% Sentinel 1 Specific Functions
def slope_correction(collection, model, buffer=0):
    #Via https://github.com/ESA-PhiLab/radiometric-slope-correction
    '''This function applies the slope correction on a collection of Sentinel-1 data
       
       :param collection: ee.Collection of Sentinel-1
       :param elevation: ee.Image of DEM
       :param model: model to be applied (volume/surface)
       :param buffer: buffer in meters for layover/shadow amsk
        
        :returns: ee.Image
    '''
    
    elevation = ee.Image('USGS/SRTMGL1_003')
    
    def _volumetric_model_SCF(theta_iRad, alpha_rRad):
        '''Code for calculation of volumetric model SCF
        
        :param theta_iRad: ee.Image of incidence angle in radians
        :param alpha_rRad: ee.Image of slope steepness in range
        
        :returns: ee.Image
        '''
        
        # create a 90 degree image in radians
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        
        # model
        nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan()
        denominator = (ninetyRad.subtract(theta_iRad)).tan()
        return nominator.divide(denominator) 
    
    
    def _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad):
        '''Code for calculation of direct model SCF
        
        :param theta_iRad: ee.Image of incidence angle in radians
        :param alpha_rRad: ee.Image of slope steepness in range
        :param alpha_azRad: ee.Image of slope steepness in azimuth
        
        :returns: ee.Image
        '''
        
        # create a 90 degree image in radians
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        
        # model  
        nominator = (ninetyRad.subtract(theta_iRad)).cos()
        denominator = (alpha_azRad.cos()
          .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos()))

        return nominator.divide(denominator)


    def _erode(image, distance):
      '''Buffer function for raster

      :param image: ee.Image that shoudl be buffered
      :param distance: distance of buffer in meters
        
      :returns: ee.Image
      '''
      
      d = (image.Not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()))
    
      return image.updateMask(d.gt(distance))
    
    
    def _masking(alpha_rRad, theta_iRad, buffer):
        '''Masking of layover and shadow
        
        
        :param alpha_rRad: ee.Image of slope steepness in range
        :param theta_iRad: ee.Image of incidence angle in radians
        :param buffer: buffer in meters
        
        :returns: ee.Image
        '''
        # layover, where slope > radar viewing angle 
        layover = alpha_rRad.lt(theta_iRad).rename('layover')

        # shadow 
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow')
        
        # add buffer to layover and shadow
        if buffer > 0:
            layover = _erode(layover, buffer)   
            shadow = _erode(shadow, buffer)  

        # combine layover and shadow
        no_data_mask = layover.And(shadow).rename('no_data_mask')
        
        return layover.addBands(shadow).addBands(no_data_mask)
                        
        
    def _correct(image):
        '''This function applies the slope correction and adds layover and shadow masks
        
        '''
        
        # get the image geometry and projection
        geom = image.geometry()
        proj = image.select(1).projection()
        
        # calculate the look direction
        heading = (ee.Terrain.aspect(image.select('angle'))
                                     .reduceRegion(ee.Reducer.mean(), geom, 1000)
                                     .get('aspect'))
                   

        # Sigma0 to Power of input image
        sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0))

        # the numbering follows the article chapters
        # 2.1.1 Radar geometry 
        theta_iRad = image.select('angle').multiply(np.pi/180)
        phi_iRad = ee.Image.constant(heading).multiply(np.pi/180)
        
        # 2.1.2 Terrain geometry
        alpha_sRad = ee.Terrain.slope(elevation).select('slope').multiply(np.pi/180).setDefaultProjection(proj).clip(geom)
        phi_sRad = ee.Terrain.aspect(elevation).select('aspect').multiply(np.pi/180).setDefaultProjection(proj).clip(geom)
        
        # we get the height, for export 
        height = elevation.setDefaultProjection(proj).clip(geom)
        
        # 2.1.3 Model geometry
        #reduce to 3 angle
        phi_rRad = phi_iRad.subtract(phi_sRad)

        # slope steepness in range (eq. 2)
        alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

        # slope steepness in azimuth (eq 3)
        alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

        # local incidence angle (eq. 4)
        theta_liaRad = (alpha_azRad.cos().multiply((theta_iRad.subtract(alpha_rRad)).cos())).acos()
        theta_liaDeg = theta_liaRad.multiply(180/np.pi)

        # 2.2 
        # Gamma_nought
        gamma0 = sigma0Pow.divide(theta_iRad.cos())
        gamma0dB = ee.Image.constant(10).multiply(gamma0.log10()).select(['VV', 'VH'], ['VV_gamma0', 'VH_gamma0'])
        ratio_gamma = (gamma0dB.select('VV_gamma0')
                        .subtract(gamma0dB.select('VH_gamma0'))
                        .rename('ratio_gamma0'))

        if model == 'volume':
            scf = _volumetric_model_SCF(theta_iRad, alpha_rRad)

        if model == 'surface':
            scf = _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)

        # apply model for Gamm0_f
        gamma0_flat = gamma0.divide(scf)
        gamma0_flatDB = (ee.Image.constant(10)
                         .multiply(gamma0_flat.log10())
                         .select(['VV', 'VH'],['VV_gamma0flat', 'VH_gamma0flat'])
                        )

        masks = _masking(alpha_rRad, theta_iRad, buffer)

        # calculate the ratio for RGB vis
        ratio_flat = (gamma0_flatDB.select('VV_gamma0flat')
                        .subtract(gamma0_flatDB.select('VH_gamma0flat'))
                        .rename('ratio_gamma0flat')
                     )

        return (image.rename(['VV_sigma0', 'VH_sigma0', 'incAngle'])
                      .addBands(gamma0dB)
                      .addBands(ratio_gamma)
                      .addBands(gamma0_flatDB)
                      .addBands(ratio_flat)
                      .addBands(alpha_rRad.rename('alpha_rRad'))
                      .addBands(alpha_azRad.rename('alpha_azRad'))
                      .addBands(phi_sRad.rename('aspect'))
                      .addBands(alpha_sRad.rename('slope'))
                      .addBands(theta_iRad.rename('theta_iRad'))
                      .addBands(theta_liaRad.rename('theta_liaRad'))
                      .addBands(masks)
                      .addBands(height.rename('elevation'))
                 )    
    
    # run and return correction
    return collection.map(_correct)

def fix_S1(ds, de, polygon, flt=True, orbit=False, gamma=False, direction='Ascending', platform='both'):
    if flt:
        #This is not log-scaled (raw power)
        S1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
    else:
        #This is log scaled (decibels)
        S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    
    S1 = S1.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
    .filter(ee.Filter.eq('instrumentMode', 'IW'))\
    .filterBounds(polygon)\
    .filterDate(ds, de)
    
    if orbit:
        S1 = S1.filter(ee.Filter.eq('relativeOrbitNumber_start', orbit))
    
    if direction == 'Ascending':
        data = S1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
    else:
        data = S1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        
    if not platform == 'both':
        data = data.filter(ee.Filter.eq('platform_number', platform))
    
    #Apply angle masking
    data = data.map(maskAngGT30)
    data = data.map(maskAngLT45)
    
    #Apply terrain correction
    if gamma:
        data = slope_correction(data, 'surface', buffer=0)
        #Choose the gamma bands and rename
        def rename(collection, which):
            def rnfx(image):
                return image.rename(['VV', 'VH']).set('system:time_start', image.get('system:time_start'))
            return collection.select(which).map(rnfx)            
                
        data = rename(data, ['VV_gamma0', 'VH_gamma0'])
    
    s1_crs = data.select('VV').first().projection()
    
    return data, s1_crs

def filter_s1(Ascending):
    def make_rat(image):
        rat = image.select('VV').divide(image.select('VH'))
        return rat.rename('VVdVH').set('system:time_start', image.get('system:time_start'))
    
    def make_rat_filt(image):
        rat = image.select('VV_filt').divide(image.select('VH_filt'))
        return rat.rename('VVdVH').set('system:time_start', image.get('system:time_start'))
    
    def make_dif(image):
        rat = image.select('VV').subtract(image.select('VH'))
        return rat.rename('VVminVH').set('system:time_start', image.get('system:time_start'))
                                       
    S1A_both = Ascending.select(['VV', 'VH']).sort('system:time_start')
    S1A_both_filt = apply_speckle_filt(S1A_both)
    
    S1A_both_focal = focal_med_filt(S1A_both)
    S1A_ratio_focal = S1A_both_focal.map(make_rat_filt)
    S1A_ratio_focal = mask_invalid(S1A_ratio_focal, -5, 5)
        
    S1A_ratio = S1A_both.map(make_rat)
    S1A_ratio_filt = S1A_both_filt.map(make_rat_filt)
    S1A_ratio_filt = mask_invalid(S1A_ratio_filt, -5, 5)
    S1A_dif = S1A_both.map(make_dif)
    
    return S1A_both, S1A_both_focal, S1A_both_filt, S1A_ratio, S1A_ratio_filt, S1A_ratio_focal

#Translate to Gamma0
def make_gamma0(image):
    angle = image.select('angle').resample('bicubic')
    return image.select('..').divide(angle.multiply(np.pi/180.0).cos()).copyProperties(image, ['system:time_start'])

#Edge masking with high/low angle
def maskAngGT30(image):
    ang = image.select(['angle'])
    return image.updateMask(ang.gt(30.63993))

def maskAngLT45(image):
    ang = image.select(['angle'])
    return image.updateMask(ang.lt(45.53993)) 

def maskAngleGT40(image):
    ang = image.select(['angle'])
    return image.updateMask(ang.gt(40))

# Some modified from here: https://code.earthengine.google.com/2ef38463ebaf5ae133a478f173fd0ab5 [Originally by Guido Lemoine]
def toNatural(img):
    return ee.Image(10.0).pow(img.select(0).divide(10.0))

def toDB(img):
    return ee.Image(img).log10().multiply(10.0)

def RefinedLee(img):
    '''
    Refined Lee Speckle Filter
    NOTE: img must be in natural units, i.e. not in dB!
    '''
    #Set up 3x3 kernels 
    weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
    kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

    mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
    variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

    #Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

    sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False)

    #Calculate mean and variance for the sampled windows and store as 9 bands
    sample_mean = mean3.neighborhoodToBands(sample_kernel)
    sample_var = variance3.neighborhoodToBands(sample_kernel)

    #Determine the 4 gradients for the sampled windows
    gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

    #And find the maximum gradient amongst gradient bands
    max_gradient = gradients.reduce(ee.Reducer.max())

    #Create a mask for band pixels that are the maximum gradient
    gradmask = gradients.eq(max_gradient)

    #duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask)

    #Determine the 8 directions
    directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
  
    #The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).Not().multiply(5))
    directions = directions.addBands(directions.select(1).Not().multiply(6))
    directions = directions.addBands(directions.select(2).Not().multiply(7))
    directions = directions.addBands(directions.select(3).Not().multiply(8))

    #Mask all values that are not 1-8
    directions = directions.updateMask(gradmask)

    #"collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum()) 

    #var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
    #Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

    sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

    #Calculate localNoiseVariance
    sigmaV = ee.Image(sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]))

    #Set up the 7*7 kernels for directional statistics
    rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))
    
    diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]])

    rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
    diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

    #Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
    dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

    #and add the bands for rotated kernels
    #for (var i=1; i<4; i++) {
    for i in range(1,4):
        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
        dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
        dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))

    #"collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum())
    dir_var = dir_var.reduce(ee.Reducer.sum())

    #And finally generate the filtered value
    varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))
    b = varX.divide(dir_var)

    result = ee.Image(dir_mean.add(b.multiply(img.subtract(dir_mean))))
    return result

def apply_speckle_filt(collection):
    bn = collection.first().bandNames().getInfo()
    def applyfx(image):
        for b in bn:
            nat = toNatural(image.select(b)) #Convert to log scale
            filt = RefinedLee(nat) #Speckle Filter
            updated = toDB(filt) #Convert back to decibels
            image = image.addBands(updated.rename(b + '_filt'))
        return ee.Image(image)
    return collection.map(applyfx)

def terrainCorrection(image):
    import numpy as np
    #Implementation by Andreas Vollrath (ESA), inspired by Johannes Reiche (Wageningen)
    #Modified from: https://gis.stackexchange.com/questions/352602/getting-local-incidence-angle-from-sentinel-1-grd-image-collection-in-google-ear
    imgGeom = image.geometry()
    srtm = ee.Image('USGS/SRTMGL1_003').clip(imgGeom) # 30m srtm 
    sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0))

    #Article ( numbers relate to chapters)
    #2.1.1 Radar geometry 
    theta_i = image.select('angle')
    phi_i = ee.Terrain.aspect(theta_i)\
        .reduceRegion(ee.Reducer.mean(), theta_i.get('system:footprint'), 1000)\
        .get('aspect')

    #2.1.2 Terrain geometry
    alpha_s = ee.Terrain.slope(srtm).select('slope')
    phi_s = ee.Terrain.aspect(srtm).select('aspect')

    #2.1.3 Model geometry
    #reduce to 3 angle
    phi_r = ee.Image.constant(phi_i).subtract(phi_s)

    #convert all to radians
    phi_rRad = phi_r.multiply(np.pi / 180)
    alpha_sRad = alpha_s.multiply(np.pi / 180)
    theta_iRad = theta_i.multiply(np.pi / 180)
    ninetyRad = ee.Image.constant(90).multiply(np.pi / 180)

    #slope steepness in range (eq. 2)
    alpha_r = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

    #slope steepness in azimuth (eq 3)
    alpha_az = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

    #local incidence angle (eq. 4)
    theta_lia = (alpha_az.cos().multiply((theta_iRad.subtract(alpha_r)).cos())).acos()
    theta_liaDeg = theta_lia.multiply(180 / np.pi)
    
    #2.2 Gamma_nought_flat
    gamma0 = sigma0Pow.divide(theta_iRad.cos())
    gamma0dB = ee.Image.constant(10).multiply(gamma0.log10())
    ratio_1 = gamma0dB.select('VV').subtract(gamma0dB.select('VH'))

    #Volumetric Model
    nominator = (ninetyRad.subtract(theta_iRad).add(alpha_r)).tan()
    denominator = (ninetyRad.subtract(theta_iRad)).tan()
    volModel = (nominator.divide(denominator)).abs()

    #apply model
    gamma0_Volume = gamma0.divide(volModel)
    gamma0_VolumeDB = ee.Image.constant(10).multiply(gamma0_Volume.log10())

    #we add a layover/shadow maskto the original implmentation
    #layover, where slope > radar viewing angle 
    alpha_rDeg = alpha_r.multiply(180 / np.pi)
    layover = alpha_rDeg.lt(theta_i)

    #shadow where LIA > 90
    shadow = theta_liaDeg.lt(85)

    #calculate the ratio for RGB vis
    ratio = gamma0_VolumeDB.select('VV').subtract(gamma0_VolumeDB.select('VH'))

    output = gamma0_VolumeDB.addBands(ratio).addBands(alpha_r).addBands(phi_s).addBands(theta_iRad)\
    .addBands(layover).addBands(shadow).addBands(gamma0dB).addBands(ratio_1)

    return image.addBands(output.select(['VV', 'VH', 'slope_1', 'slope_2'], ['VV', 'VH', 'layover', 'shadow']),None, True)

def NBMI(collection):
    ''' 
    Calculate the NBMI based on Shoshaney et al. 2000 Int. J. Remote Sensing
    NBMI = A + B / A - B, where A and B are backscatter at two dates. 
    
    NOTE: It is important to do a tight spatial filtering first!
    '''
    def dif(f):
        #sdate = ee.Date(ee.Image(ee.List(f).get(0)).get('system:time_start'))
        f = ee.Number(f)
        simage = ee.Image(ic_list.get(f))
        snext = ee.Image(ic_list.get(f.add(1)))
        sdate = simage.get('system:time_start')
        
        su = simage.add(snext)
        dif = simage.subtract(snext)
        ret = su.divide(dif).set('system:time_start', sdate)
        
        return ret
    
    ic_list = collection.sort('system:time_start').toList(collection.size())
    seq = ee.List.sequence(0, collection.size().subtract(2)) #Drop the final image since we don't have the next to subtract
    ic_diff = seq.map(dif)
    
    return ee.ImageCollection.fromImages(ic_diff)

#%% Online Filtering
def apply_SG(collect, geom, window_size=7, imageAxis=0, bandAxis=1, order=3):
    '''
    Apply a Savitzky-Golay filter to a time series collection (pixelwise)

    '''
    def prep_SG(img):
        #Add predictors for SG fitting, using date difference
        #We prepare for order 3 fitting, but can be adapted to lower order fitting later on
        dstamp = ee.Date(img.get('system:time_start'))
        ddiff = dstamp.difference(ee.Date(ds), 'hour')
        return img.addBands(ee.Image(1).toFloat().rename('constant'))\
        .addBands(ee.Image(ddiff).toFloat().rename('t'))\
        .addBands(ee.Image(ddiff).pow(ee.Image(2)).toFloat().rename('t2'))\
        .addBands(ee.Image(ddiff).pow(ee.Image(3)).toFloat().rename('t3'))\
        .set('date', dstamp)
    
    def getLocalFit(i):
        #Get a slice corresponding to the window_size of the SG smoother
        subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window_size).int())
        #predictors = subarray.arraySlice(bandAxis, 2, 2 + order + 1)
        predictors = subarray.arraySlice(bandAxis, 1, 1 + order + 1) #Here for a one-variable case
        response = subarray.arraySlice(bandAxis, 0, 1)
        coeff = predictors.matrixSolve(response)
        
        coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener)
        return coeff 
    
    def apply_SG_sub(i):
        ref = ee.Image(c.get(ee.Number(i).add(ee.Number(half_window))))
        return getLocalFit(i).multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum()).copyProperties(ref)
    
    half_window = (window_size - 1)/2
    if order == 3:
        coeffFlattener = [['constant', 'x', 'x2', 'x3']]
        indepSelectors = ['constant', 't', 't2', 't3']
    elif order == 2:
        coeffFlattener = [['constant', 'x', 'x2']]
        indepSelectors = ['constant', 't', 't2']
        
    collect_coeffs = collect.map(prep_SG)
    array = collect_coeffs.toArray() #THIS STEP IS EXPENSIVE
    c = collect_coeffs.toList(collect_coeffs.size())
    runLength = ee.List.sequence(0, c.size().subtract(window_size))
    
    sg_series = runLength.map(apply_SG_sub)
    
    #Drop the null values
    sg_sliced = sg_series.slice(half_window, sg_series.size().subtract(half_window))
    sg_series = ee.ImageCollection.fromImages(sg_sliced)
    
    def minmax(img):
        #Map reducer to get global min/max NDVI (for filtering purposes)
        bn = img.bandNames().get(0)
        minMax =  img.reduceRegion(ee.Reducer.minMax(), geom, 1000)
        return img.set({'roi_min': minMax.get(ee.String(bn).cat('_min')),'roi_max': minMax.get(ee.String(bn).cat('_max'))}).set('date', img.get('date'))
    
    sg = sg_series.map(minmax)                    
    
    return sg.filterMetadata('roi_min', 'not_equals', None).filterMetadata('roi_max', 'not_equals', None)

#%% Conversion to Python Time Series
def export_to_pandas(collection, clipper, aggregation_scale, save=None, med='median'):
    '''
    Takes an ImageCollection, an Earth Engine Geometry, and an aggregation scale (e.g., 30m for Landsat, 250m for MODIS, etc)
    
    Returns a pandas time series for the mean/median and standard deviation values over the 
    aggregation area. 
    
    Optionally saves those time series to a CSV file    
    
    '''
    import pandas as pd, numpy as np
    
    def createTS(image):
        date = image.get('system:time_start')
        if med == 'median':
            value = image.reduceRegion(ee.Reducer.median(), clipper, aggregation_scale)
        elif med == 'mean':
            value = image.reduceRegion(ee.Reducer.mean(), clipper, aggregation_scale)
        std = image.reduceRegion(ee.Reducer.stdDev(), clipper, aggregation_scale)
        ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'Mn': value, 'STD': std})
        return ft
    
    TS = collection.filterBounds(clipper).map(createTS)
    dump = TS.getInfo()
    fts = dump['features']
    out_vals = np.empty((len(fts)))
    out_dates = []
    out_std = np.empty((len(fts)))
    
    for i, f in enumerate(fts):
        props = f['properties']
        date = props['date']
        try:
            val = list(props['Mn'].values())[0]
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
        df = pd.DataFrame({'mn':out_vals, 'std':out_std, 'time':out_dates})
        df.to_csv(save, index=False)
        print(save)
    return ser, serstd
    
def percentile_export(collection, percentile, clipper, aggregation_scale=30, save=None):
    '''
    Get a time series at a certain percentile
    '''
    
    import pandas as pd, numpy as np
    
    def createTS(image):
        date = image.get('system:time_start')
        value = image.reduceRegion(ee.Reducer.percentile([percentile]), clipper, aggregation_scale)
        ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'PCT': value})
        return ft
    
    TS = collection.filterBounds(clipper).map(createTS)
    dump = TS.getInfo()
    fts = dump['features']
    out_vals = np.empty((len(fts)))
    out_dates = []
    
    for i, f in enumerate(fts):
        props = f['properties']
        date = props['date']
        try:
            val = list(props['PCT'].values())[0]
        except:
            val = np.nan
        out_vals[i] = val
        out_dates.append(pd.Timestamp(date))
    
    ser = pd.Series(out_vals, index=out_dates)
    if save:
        df = pd.DataFrame({'p' + str(percentile):out_vals, 'time':out_dates})
        df.to_csv(save + '.csv', index=False)
        print(save)
    return ser


#%% Short Example
# minx2, maxx2 = 18, 20
# miny2, maxy2 = -20, -18
# from shapely.geometry import Polygon
# geom = Polygon([[minx2, maxy2], [maxx2, maxy2], [maxx2, miny2], [minx2, miny2]])

# search_area = gee_geometry_from_shapely(geom)

# L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").filter(ee.Filter.date('2013-01-01', '2020-01-01'))\
#     .filterBounds(search_area)\
#     .filter(ee.Filter.lt('CLOUD_COVER', 15))\
#     .map(maskL8sr)\
#     .sort('system:time_start')

# NDVI = L8.map(NDVI_L8).select('NDVI')

# harmonic, phase, amplitude, rsq = fit_multi_harmonic(detrend(NDVI), 5)
# run_export(phase, 'epsg:4326', 'Landsat8_HarmonicPhase', scale=30, region=search_area, maxPixels=1e12)
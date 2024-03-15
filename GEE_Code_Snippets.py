# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 11:16:43 2020

@author: tsmith
"""

import ee
ee.Initialize()
import numpy as np
import os

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
    
def export_timeseries_todrive(collection, filename, scale, region, folder=None, agg_fx=ee.Reducer.median()):
    '''
    Runs an export function on GEE servers
    '''
    #Create the time series    
    def createTS(image):
        date = image.get('system:time_start')
        value = image.reduceRegion(agg_fx, region, scale)
        ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'val': value})
        return ft
    
    TS = collection.filterBounds(region).map(createTS)
    
    task_config = {'description': filename, 'fileFormat': 'CSV'}
    if folder:
        task_config['folder'] = folder
    task = ee.batch.Export.table.toDrive(TS, **task_config)
    task.start()
        
def export_collection(collection, region, prefix, crs=None, scale=100, start_image=0, \
                      max_images=None, folder=None, namelist=None):
    '''
    Exports all images within an image collection for a given region. All files named by a prefix (given)
    and their image date (formatted YYYYMMDD). 
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
    
    if namelist:
        #Run through provided prefixes (e.g., one image for each month or year in a collection)
        for i, name in enumerate(namelist):
            image = ee.Image(image_list.get(i))
            output_name = prefix + '_' + name + '_' + str(scale) + 'm'
            run_export(image, crs=crs, filename=output_name, scale=scale, region=region, folder=folder)
            print('Started export for image ' + str(i) + ' (' + name + ')')
            
    else:
        #Run a list from the starting image to the number you want using the date of the image in the name
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
                    #Otherwise simply name by the image collection date
                    date = image.get('system:time_start')
                    date_name = ee.Date(date).format('YYYYMMdd').getInfo()
                output_name = prefix + '_' + date_name + '_' + str(scale) + 'm'
                run_export(image, crs=crs, filename=output_name, scale=scale, region=region, folder=folder)
                print('Started export for image ' + str(i) + '(' + date_name + ')')
       
def split_export(image, namebase, crs=None, scale=500, minx=-180, maxx=180, miny=-80, maxy=80, step=20, gdrive=None, folder=None, **kwaargs):
    '''
    Split a large (e.g., global) job into smaller pieces. Can help with memory issues/speed of processing.
    '''
    if not crs:
        crs = ee.Projection('EPSG:4326')
    for i in range(minx,maxx,step):
        for j in range(miny,maxy,step):
            imax = i + step
            jmax = j + step
            if jmax > maxy:
                jmax = maxy
            if imax > maxx:
                imax = maxx
            name = namebase + '_%i_%i_%i_%i_WGS84' % (i,imax,j,jmax)
    
            roi = ee.Geometry.Polygon([[i, j], [i, jmax], [imax, jmax], [imax, j]])
            if gdrive:
                if os.path.exists(gdrive + folder + '/' + name + '.tif'):
                    print(gdrive + folder + '/' + name + '.tif exists in gdrive!')
                    pass
                else:
                    print('Starting', name)
                    run_export(image, crs=crs, filename=name, region=roi, scale=scale, folder=folder, **kwaargs)
            else:
                print('Starting', name)
                run_export(image, crs=crs, filename=name, region=roi, scale=scale, folder=folder, **kwaargs)
        
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
    
    #Turn them into a single ImageCollection
    output = ee.ImageCollection.fromImages(merged_images)
    return output

def gee_geometry_from_shapely(geom, crs='epsg:4326'):
    """ 
    Simple helper function to take a shapely geometry and a coordinate system and convert them to a 
    Google Earth Engine Geometry.
    """
    from shapely.geometry import mapping
    ty = geom.geom_type
    if ty == 'Polygon':
        return ee.Geometry.Polygon(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)
    elif ty == 'Point':
        return ee.Geometry.Point(ee.List(mapping(geom)['coordinates']), proj=crs)    
    elif ty == 'MultiPolygon':
        return ee.Geometry.MultiPolygon(ee.List(mapping(geom)['coordinates']), proj=crs, evenOdd=False)
    elif ty == 'MultiPoint':
        return ee.Geometry.MultiPoint(ee.List(mapping(geom)['coordinates']), proj=crs)
    
def geopandas_to_earthengine(gdf, crs='epsg:4326'):
    '''
    Helper to quickly create a set of earth engine geometries from a geodataframe
    '''
    feats = []
    for _, i in gdf.iterrows():
        ee_geom = gee_geometry_from_shapely(i.geometry, crs=crs)
        atts = i.to_dict()
        del atts['geometry'] #Make sure geometry isn't duplicated
        feats.append(ee.Feature(ee_geom, atts))
    
    return ee.FeatureCollection(feats)

def create_ee_fc(fid, buffer=0):
    ''' Secondary (faster) way to create an earth engine feature collection '''
    import json
    import geopandas as gpd
    gdf = gpd.read_file(fid)
    geo_json = gdf.to_json()
    ee_samp = ee.FeatureCollection(json.loads(geo_json))
    if buffer:
        def buff(f):
            return f.buffer(buffer)
        ee_samp = ee_samp.map(buff)
    return ee_samp

def create_ee_fc_split(fid, n_chunks=20):
    ''' Create a list of earth engine feature collections based on chunking '''
    import json
    import geopandas as gpd
    gdf = gpd.read_file(fid)
    dflist = np.array_split(gdf, n_chunks)
    out_fc_list = []
    for df in dflist:
        geo_json = df.to_json()
        ee_samp = ee.FeatureCollection(json.loads(geo_json))
        out_fc_list.append(ee_samp)
    return out_fc_list

def export_features(fc, filename, folder=None, output_format='CSV'):
    '''
    Runs an export function on GEE servers for feature (tabular) data. Can also export GeoJSON.
    '''
    task_config = {'fileNamePrefix': filename, 'fileFormat': output_format}
    if folder:
        task_config['folder'] = folder
    task = ee.batch.Export.table.toDrive(fc, filename, **task_config)
    task.start()
    
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

def reduce_res(ic, out_crs, scale, agg_fx):
    ''' Change the scale of an image based on a given agg_fx (mean, mode, etc)'''
    def red_(image):
        #input_crs = image.projection()
        reduce_image = image.reduceResolution(reducer=agg_fx, maxPixels=65535)
        return reduce_image.reproject(crs=out_crs, scale=scale).set('system:time_start', image.get('system:time_start'))
    return ic.map(red_)

def reduce_imageres(image, out_crs, scale, agg_fx):
    ''' Change the scale of an image based on a given agg_fx (mean, mode, etc)'''
    reduce_image = image.reduceResolution(reducer=agg_fx, maxPixels=65535)
    return reduce_image.reproject(crs=out_crs, scale=scale).set('system:time_start', image.get('system:time_start'))

def rescale(ic, scale, add=0):
    '''
    Rescale all images in an image collection by a given scale factor and addative factor. 
    ASSUMES ALL BANDS SCALED EQUALLY
    '''
    def r(image):
        return image.multiply(scale).add(add).set('system:time_start', image.get('system:time_start'))
    return ic.map(r)

def rename(ic, name):
    '''
    Rename all images in a collection to a specified new name.
    '''
    def rn(image):
        return ee.Image(image).rename(name)
    return ic.map(rn)

def join_fc(f1, f2, on):
    '''
    Quick and dirty join of feature collections on an arbitrary field. 
    '''
    Filter = ee.Filter.equals(leftField=on, rightField=on)
    simpleJoin = ee.Join.simple()
    simpleJoined = simpleJoin.apply(f1, f2, Filter)
    return simpleJoined

def join_c(c1, c2, on='system:index'):
    '''
    Quick and dirty join of Image Collections on an arbitrary field. 
    '''
    filt = ee.Filter.equals(leftField=on, rightField=on) 
    innerJoin = ee.Join.inner() #initialize the join
    innerJoined = innerJoin.apply(c1, c2, filt) #This is a FEATURE COLLECTION
    def combine_joined(feature):
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
    
    joined_collect = ee.ImageCollection(innerJoined.map(combine_joined))
    return joined_collect

def simple_difference(ic, b1, b2, outname='dif'):
    '''
    Simple difference between two bands in an image collection
    '''
    def sub(image):
        i1 = image.select(b1)
        i2 = image.select(b2)
        return i1.subtract(i2).rename(outname).set('system:time_start', image.get('system:time_start')).set('system:footprint', image.get('system:footprint'))
    return ic.map(sub)

def rolling_apply(ic, windowSize, fx, var=None, ts='year'):
    '''
    Apply a function over a moving window (windowsize is in days)
    '''
    if var:
        ic = ic.select(var)

    def apply_roller(dates, fx):
        def roller(t):
            t = ee.Date(t) #Get the date
            #Use the date to create a search window of windowSize in both directions
            window = ic.filterDate(t.advance(-windowSize,ts),t.advance(windowSize,ts))
            
            #Apply the function
            result = fx(window)
            return result
        
        vals = dates.map(roller)
        return vals
    
    dates = ee.List(ic.aggregate_array('system:time_start'))
    vals = apply_roller(dates, fx)
    return ee.ImageCollection.fromImages(vals)

def get_local_utm(geometry):
    '''
    Take a given geometry (shapely geometry object) and get the local UTM code as a EE Projection.
    NOTE: Works on/off-line in tandem, so cannot be mapped over an image collection!
    '''
    def choose_utm_zone(geometry):
        ''' Generate an EPSG code for UTM projection for a given lat/lon'''
        import math
        lon = geometry.centroid.x
        lat = geometry.centroid.y

        #https://gis.stackexchange.com/questions/269518/auto-select-suitable-utm-zone-based-on-grid-intersection
        utm_band = str(int(math.floor((lon + 180) / 6 ) % 60) + 1)
        if len(utm_band) == 1:
            utm_band = '0' + utm_band
        if lat >= 0:
            epsg_code = '326' + utm_band
        else:
            epsg_code = '327' + utm_band

        return int(float(epsg_code))
    
    epsg_code = choose_utm_zone(geometry)
    return ee.Projection('epsg:' + str(epsg_code)) #Turn it into an earth-engine projection

def get_utm_from_image(image):
    '''
    Use an image's geometry to get the proper UTM zone. Useful for processing large batches
    of Landsat/Sentinel data that all have different projections over a large area.
    Returns ee projection item.
    NOTE: This is all done server-side, so can be used within .map() functions without an issue!
    '''
    fp = ee.Geometry(image.get('system:footprint'))
    ct = fp.centroid().coordinates()
    lon, lat = ee.Number(ct.get(0)), ee.Number(ct.get(1)) #x,y
    
    #utm_band = str(int(math.floor((lon + 180) / 6 ) % 60) + 1)
    
    step1 = lon.add(180).divide(6).floor()
    step2 = step1.mod(60).add(1)
    utm_band = ee.Number(step2).int().format('%s')
    
    #if len(utm_band) == 1:
    #    utm_band = ee.String('0').cat(utm_band)
    #if lat >= 0:
    #    epsg_code = ee.String('326').cat(utm_band)
    #else:
    #    epsg_code = ee.String('327').cat(utm_band)
    
    #NOTE: THIS DOESNT SUPPORT POLAR COORDINATE SYSTEMS RIGHT NOW
    epsg_code = ee.Algorithms.If(lat.gte(0), ee.String('326').cat(utm_band), ee.String('327').cat(utm_band))
    
    proj = ee.Projection(ee.String('epsg:').cat(ee.String(epsg_code)))
    return proj

def createConstantBand(image):
    return ee.Image(1).addBands(image).set('system:time_start', image.get('system:time_start'))

#%%
def disentangle_to_daily(ic, ys, ye):
    '''Convert 8-day data to daily data'''
    
    def yearly_data(year):
        return ic.filter(ee.Filter.calendarRange(year, year, 'year'))
    
    #years = ee.List(list(range(ys,ye)))
    #yearly_ic = years.map(yearly_data)
    
    def create_image(date):
        img = ee.Image(0).set('system:time_start', date)
        return img
    
    ly = list(range(1960,2028,4))
    
    ys, ye = int(ys), int(ye)
    daily_data = []
    for yr in range(ys,ye+1):
        yearly_ic = yearly_data(yr)
        yearly_ic_list = ee.ImageCollection.toList(yearly_ic, yearly_ic.size())
        
        if yr in ly:
            nr_days = 366
        else:
            nr_days = 365
        
        #Create empty images for each day of the year
        t0 = ee.Date(str(yr) + '-01-01')
        days = ee.List([t0.advance(x,'day') for x in range(0,nr_days)])
        images = ee.ImageCollection(days.map(create_image))
        
        def add_to_images(images, ds, de, avg):
            subset = images.filterDate(ds, de)
            def add(image):
                return image.add(avg).set('system:time_start', image.get('system:time_start'))
            return subset.map(add)
        
        #Get the image for each 8-day period
        nr_images = ee.Number(yearly_ic.size())
        nr_images_local = nr_images.getInfo()
        for i in range(nr_images_local):
            data = yearly_ic_list.get(i)
            
            if not i == nr_images_local-1:
                avg = ee.Image(data).divide(8)
            else:
                if nr_days == 366:
                    avg = ee.Image(data).divide(6)
                else:
                    avg = ee.Image(data).divide(5)
            
            ds = t0.advance(i*8, 'day')
            de = t0.advance((i+1)*8, 'day')
            updated_images = add_to_images(images, ds, de, avg)
            daily_data.append(updated_images)
        
    output = daily_data[0]
    for i in daily_data[1:]:
        output = output.merge(i)
    return output

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

def piecewise_detrend(ic, start_year, end_year, fit_period=3):
    #Get the length of the sides from the fitting period
    len_sides = int((fit_period - 1) / 2)
    
    def create_yearlist(yr, len_sides):
        '''
        Create a list of years based on years before/after a given center year
        '''
        yearlist = [yr]
        for i in range(1,len_sides+1):
            yearlist.append(yr - i)
            yearlist.append(yr + i)
        yearlist.sort()
        return yearlist
    
    def multi_year_data(yearlist):
        '''
        For a list of years, return only that subset of the data
        '''
        flist = []
        for year in yearlist:
            flist.append(ee.Filter.calendarRange(year, year, 'year'))
        if len(flist) == 1:
            filt = flist[0]
        else:
            filt = ee.Filter.Or(flist)
        return filt
    
    output = []
    for year in range(start_year, end_year + 1):
        #Get this list of years to use as a filter
        yearlist = create_yearlist(year, len_sides)
        filt = multi_year_data(yearlist)
        
        #Subset the data and fit a harmonic of given order
        subset = ic.filter(filt)
        
        #Detrend that subset
        subset_dt = detrend(subset)
        
        #Retain only the given year
        filt = ee.Filter.calendarRange(year, year, 'year')
        dt_year = subset_dt.filter(filt)
        output.append(dt_year)
    
    #Merge the years of data back together
    base = output[0]
    for i in output[1:]:
        base = base.merge(i)
    return base

def detrend_nonlin(collection, order=3):    
    def addVariables(image):
        date = ee.Date(image.get('system:time_start'))
        years = date.difference(ee.Date('1970-01-01'), 'year')
        t = ee.Image(years).rename('t')
        t2 = ee.Image(years).pow(ee.Image(2)).rename('t2').toFloat()
        t3 = ee.Image(years).pow(ee.Image(3)).rename('t3').toFloat()
        t4 = ee.Image(years).pow(ee.Image(4)).rename('t4').toFloat()
        t5 = ee.Image(years).pow(ee.Image(5)).rename('t5').toFloat()
        t6 = ee.Image(years).pow(ee.Image(6)).rename('t6').toFloat()
        t7 = ee.Image(years).pow(ee.Image(7)).rename('t7').toFloat()
        t8 = ee.Image(years).pow(ee.Image(8)).rename('t8').toFloat()
        return image.addBands(ee.Image.constant(1)).addBands(t).addBands(t2)\
                .addBands(t3).addBands(t4).addBands(t5).addBands(t6).addBands(t7).addBands(t8).float()

    idps = ['constant', 't', 't2', 't3', 't4', 't5', 't6', 't7', 't8']
    idps = idps[:order+1]
    independents = ee.List(idps)
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
    return detrended.select('detrend')

def rolling_detrend(ic, fit_period=5, window_unit='year', var=None):
    '''
    Apply a function over a moving window
    '''
    if var:
        ic = ic.select(var)
        
    if window_unit == 'year':
        fit_period = fit_period * 365
        window_unit = 'day'

    len_sides = int((fit_period - 1) / 2)

    def apply_roller(dates):
        def roller(t):
            t = ee.Date(t) #Get the date
            #Use the date to create a search window of windowSize in both directions
            window = ic.filterDate(t.advance(-len_sides,window_unit),t.advance(len_sides,window_unit))
            
            #Detrend over that window
            dt = detrend(window)
            
            #Select only the center date
            img = ee.Image(dt.filterDate(t).first())
            
            return img.set('system:time_start', t).rename('dt').float()
        
        vals = dates.map(roller)
        return vals
    
    dates = ee.List(ic.aggregate_array('system:time_start'))
    vals = apply_roller(dates)
    
    return ee.ImageCollection.fromImages(vals)

def deseason_data(ic, bn=None, detrend_type='moving', deseason=True, detrend_fit_period=5, harmonic='single', start_year=2000, end_year=2022, \
                  harmonic_fit_period=3, harmonic_order=3):
    '''
    Top-level for de-seasoning data. Takes multiple options for detrending and deseasoning. 
    '''
    if not bn:
        img = ic.first()
        bn = img.bandNames().get(0)
     
    #Get rid of linear trends
    if detrend_type == 'single':
        detrend_ic = detrend(ic) #This is a simple non-moving detrender
    elif detrend_type == 'moving':
        detrend_ic = piecewise_detrend(ic, start_year, end_year, detrend_fit_period) #Fit 5-year ramps to detrend, three years for harmonics
    else:
        detrend_ic = ic
        
    if deseason:
        #Find the seasonal component
        if harmonic == 'single':
            harmonic, _, _, _ = fit_multi_harmonic(detrend_ic, harmonics=harmonic_order)
        elif harmonic == 'moving':
            harmonic = moving_harmonic(detrend_ic, start_year, end_year, fit_period=harmonic_fit_period, harmonic_order=harmonic_order)

        #Rename the detrended data as dt
        dt = rename(detrend_ic, 'dt')

        #Get seasonality via harmonics
        seasonality = rename(harmonic, 'seas')

        #Subtract the seasonal signal from the detrended signal to produce the final residuals
        joint_vals = join_c(dt, seasonality, on='system:time_start')
        deseasoned = simple_difference(joint_vals, 'dt', 'seas')

        deseasoned = rename(deseasoned, 'ds')
    else:
        deseasoned = rename(detrend_ic, 'ds')
    
    return deseasoned
    
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
    trend = coll.select(independents.add(dependent)).reduce(ee.Reducer.linearRegression(independents.length(), 1), 4)
    
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
    
def annual_trend(ic, mask=None, first_year=1999, last_year=2023, robust=True, var='median', reducer=ee.Reducer.median()): #Takes in deseasoned/detrended data!
    save = []
    for yr in range(first_year, last_year):
        sy, ey = yr, yr + 1 
        sd, ed = ee.Date(str(sy) + '-10-01'), ee.Date(str(ey) + '-10-01') #Clip out water-years (oct-oct)
        subset = ic.filterDate(sd, ed)
        out_date = ee.Date(str(ey) + '-04-01')
        
        #Get the median for that year
        mn = subset.reduce(reducer).rename(var).set('system:time_start', out_date).float()
        if mask:
            mn = mn.updateMask(mask)
        
        #Create a constant to make regressions with later
        const = ee.Image.constant(1).rename('constant').set('system:time_start', out_date).float()
        base_yr = ee.Image(yr).rename('year').set('system:time_start', out_date).float() #Scale to make regression more stable      
        
        #Create one multiband image
        output = const.addBands(base_yr).addBands(mn)
        if mask:
            output = output.updateMask(mask)
        save.append(output.set('system:time_start', out_date))
    
    #Turn everything into a single collection
    save = ee.ImageCollection.fromImages(save)
    
    #Now do some regressions
    if robust:
        fit = save.select(['constant', 'year', var]).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1), 4)
    else:
        fit = save.select(['constant', 'year', var]).reduce(ee.Reducer.linearRegression(numX=2, numY=1), 4)
    slope = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['constant', 'trend']]).select('trend')
    if mask:
        slope = slope.updateMask(mask)
    #rsq = get_rsq(save, fit, 'count', ee.List(['constant', 'year'])).rename('rsq')
    
    return slope

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
    bn = img.bandNames().get(0)
    
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
def fit_harmonic(collection, bn=None):
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
    if not bn:
        img = collection.first()
        #bn = img.bandNames().getInfo()[0]
        bn = img.bandNames().get(0)
        
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

def fit_multi_harmonic(collection, harmonics=3, bn=None):
    """
    Function to fit a complex harmonic model to an ImageCollection. Uses any number of desired harmonics. 
    Original function adapted from:
    https://code.earthengine.google.com/2669122497313113fc4bb81bc8352828
    
    Returns the harmonic-smoothed dataset, the phase of the harmonic, the amplitude of the harmonic
    and a simple R squared value for the fit. 

    """
    import numpy as np

    #Get the name of the image band
    if not bn:
        img = collection.first()
        #bn = img.bandNames().getInfo()[0]
        bn = img.bandNames().get(0)

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
    harmonicTrend = harmonic_coll.select(independents.add(bn)).reduce(ee.Reducer.linearRegression(independents.length(), 1), 4)
    
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

    return fittedHarmonic.select('fitted'), multiphase, multiamp, rsq

def moving_harmonic(ic, start_year, end_year, bn=None, fit_period=3, harmonic_order=3):
    '''
    Fit a harmonic in multiple parts -- can be done year by year, or by fitting the harmonics to a 
    given multi-year time window. Fit period must be odd (or zero)!
    '''
    
    #Get the length of the sides from the fitting period
    len_sides = int((fit_period - 1) / 2)
    
    def create_yearlist(yr, len_sides):
        '''
        Create a list of years based on years before/after a given center year
        '''
        yearlist = [yr]
        for i in range(1,len_sides+1):
            yearlist.append(yr - i)
            yearlist.append(yr + i)
        yearlist.sort()
        return yearlist
    
    def multi_year_data(yearlist):
        '''
        For a list of years, return only that subset of the data
        '''
        flist = []
        for year in yearlist:
            flist.append(ee.Filter.calendarRange(year, year, 'year'))
        if len(flist) == 1:
            filt = flist[0]
        else:
            filt = ee.Filter.Or(flist)
        return filt
    
    output = []
    for year in range(start_year, end_year + 1):
        yearlist = create_yearlist(year, len_sides)
        filt = multi_year_data(yearlist)
        
        #Subset the data and fit a harmonic of given order
        subset = ee.ImageCollection(ic.filter(filt))
        
        #Detrend that subset
        #subset_dt = detrend(subset)
        harmonic, phase, amplitude, rsq = fit_multi_harmonic(subset, harmonics=harmonic_order, bn=bn)
        
        #Retain only the given year
        filt = ee.Filter.calendarRange(year, year, 'year')
        harmonic_year = harmonic.filter(filt)
        output.append(harmonic_year)
    
    #Merge the years of data back together
    base = output[0]
    for i in output[1:]:
        base = base.merge(i)
    return base

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

def lt_monthly_mean(ic):
    '''
    Long-term average for each month of a given image collection.
    '''
    def apply_filt(m):
        filt = ic.filter(ee.Filter.calendarRange(m, m, 'month'))
        return filt.reduce(ee.Reducer.mean()).set('month', m)
    l = ee.List([1,2,3,4,5,6,7,8,9,10,11,12])
    monthly_means = l.map(apply_filt)
    return monthly_means

def match_time_res(in_ic, match_ic, nr_ts, timeframe, agg_fx):
    '''
    Create a new collection with time spacing based on a reference. Pools the previous 'nr_ts' 'timeframe's (e.g., 16 days)
    with a given aggregator (e.g., ee.Reducer.mean()). Useful for summing up previous weeks of rain, for example.
    '''
    def create_new_composite(image):
        #Get date from input imagecollection
        ed = ee.Date(image.get('system:time_start'))
        #Get previous date timestamp
        sd = ed.advance(-1*nr_ts, timeframe)
        #Filter the match image collection to those dates
        subset_match_ic = match_ic.filterDate(sd, ed)
        #Aggregate
        agg = subset_match_ic.reduce(agg_fx)
        return agg.set('system:time_start', image.get('system:time_start'))
    
    fix_ic = in_ic.map(create_new_composite)
    return fix_ic

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
    bn = ee.Image(collection.reduce(ee.Reducer.mean())).bandNames().get(0)#.getInfo()
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

def reduceFit(collection, time_step='month'):
    t = collection.get('system:time_start')
    bn = collection.first().bandNames().get(0)#.getInfo()[0]
    def createTimeBand(image):
        date = ee.Date(image.get('system:time_start'))
        years = date.difference(ee.Date('1970-01-01'), time_step)
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

    bn = collection.first().bandNames().get(0)#.getInfo()[0]
    
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
    if season == 'MAM':
        for m in [3, 4, 5]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
    if season == 'JJA':
        for m in [6, 7, 8]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
    if season == 'SON':
        for m in [9, 10, 11]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
    if season == 'DJFMAM': 
        for m in [12, 1, 2, 3, 4, 5]:
            allfilts.append(ee.Filter.calendarRange(m, m, 'month'))
    if season == 'JJASON': 
        for m in [6, 7, 8, 9, 10, 11]:
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
        
def linearFit_time(collection, unit='year'):
    bn = collection.first().bandNames()
    
    def createTimeBand(image):
        date = ee.Date(image.get('system:time_start'))
        td = date.difference(ee.Date('1970-01-01'), unit)
        return image.addBands(ee.Image(td).rename('t')).set('system:time_start', image.get('system:time_start'))
    
    def recast(image):
        return image.float().set('system:time_start', image.get('system:time_start'))
    
    prepped = collection.map(createConstantBand).map(recast)
    prepped2 = prepped.map(createTimeBand).map(recast)
    
    var = ee.List(['constant', 't']).cat(bn)
        
    fit = prepped2.select(var).reduce(ee.Reducer.robustLinearRegression(numX=2, numY=1))
    lrImage = fit.select(['coefficients']).arrayProject([0]).arrayFlatten([['constant', 'trend']])

    return lrImage.select('trend'), fit.select(['residuals']).arrayFlatten([['residuals']])
        
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

def otsu(histogram, fixed=False):
    #This function takes a histogram and provides the proper threshold for splitting it into two components.
    #Implementation modified from: https://medium.com/google-earth/otsus-method-for-image-segmentation-f5c48f405e
    if not fixed:
        #E.g., reducer=ee.Reducer.histogram(255,2)
        counts = ee.Array(ee.Dictionary(histogram).get('histogram'))
        means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
    else:
        #E.g., reducer=ee.Reducer.fixedHistogram(min_bin, max_bin, steps)
        both = ee.Array(histogram).transpose()
        means = ee.Array(both.toList().get(0))
        counts = ee.Array(both.toList().get(1))
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sum_ = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean = sum_.divide(total);

    indices = ee.List.sequence(1, size)

    #Compute between sum of squares, where each mean partitions the data.
    def iterate(i):
        aCounts = counts.slice(0, 0, i)
        aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
        aMeans = means.slice(0, 0, i)
        aMean = aMeans.multiply(aCounts)\
            .reduce(ee.Reducer.sum(), [0]).get([0])\
            .divide(aCount)
        bCount = total.subtract(aCount)
        bMean = sum_.subtract(aCount.multiply(aMean)).divide(bCount)
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(bCount.multiply(bMean.subtract(mean).pow(2)))
    bss = indices.map(iterate)

    #Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1])

#%% Data Import and Cleaning Functions
def kNDVI(image):
    ''' Compute kNDVI on a given image -- MUST BE THE ONLY BAND!'''
    ndvi2 = image.pow(2)
    kndvi = ndvi2.tanh()
    return kndvi.set('system:time_start', image.get('system:time_start'))

### GPM
def mask_GPM(image):
    mask = image.gt(0.1)
    precip_filt = image.updateMask(mask).set('system:time_start', image.get('system:time_start'))
    return precip_filt

### Landsat
def applyScaleFactors(image):
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).float()
    thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0).float()
    return image.addBands(opticalBands, overwrite=True).addBands(thermalBands, overwrite=True).set('system:time_start', image.get('system:time_start'))

def NDSI_L7(image):
    ndsi = image.normalizedDifference(['B2', 'B5']).rename('NDSI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndsi)

def NDSI_L8(image):
    ndsi = image.normalizedDifference(['B3', 'B6']).rename('NDSI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndsi)

def NDSI_L8_C2(image):
    ndsi = image.normalizedDifference(['SR_B3', 'SR_B6']).rename('NDSI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndsi)

def NDVI_L7(image):
    ndvi = image.normalizedDifference(['B4', 'B3']).rename('NDVI').set('system:time_start', image.get('system:time_start'))
    return image.addBands(ndvi)

def MSAVI_L8(image):
    msavi2 = image.select('B5').multiply(2).add(1).subtract(image.select('B5').multiply(2).add(1).pow(2)\
                .subtract(image.select('B5').subtract(image.select('B4')).multiply(8)).sqrt())\
                .divide(2)\
                .set('system:time_start', image.get('system:time_start'))\
                .rename('MSAVI2')
    return image.addBands(msavi2)

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

def mask_landsat_c2(image):
    #Bits 3 and 5 are cloud shadow and cloud, respectively.
    cloudShadowBitMask = (1 << 3)
    cloudsBitMask = (1 << 4)
    cirrusBitMask = (1 << 2)
    #Get the pixel QA band.
    qa = image.select('QA_PIXEL')
    #Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0)).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
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

def maskL57toa(image):
    #via: https://gis.stackexchange.com/questions/274612/apply-a-cloud-mask-to-a-landsat8-collection-in-google-earth-engine-time-series
    #via: https://gis.stackexchange.com/questions/292835/using-cloud-confidence-to-create-cloud-mask-from-landsat-8-bqa
    #via: https://gis.stackexchange.com/questions/271483/how-to-apply-a-cloud-mask-in-google-earth-engine-landsat-5-tm-8-day-ndvi-compo
    def extractQABits(qaBand, bitStart, bitEnd):
        numBits = bitEnd - bitStart + 1
        qaBits = qaBand.rightShift(bitStart).mod(2**numBits)
        return qaBits
    
    qa_band = image.select('BQA')
    cloud_shadow = extractQABits(qa_band, 7, 8)
    cloud = extractQABits(qa_band, 4, 4)
    mask_cloud = cloud.neq(1)
    
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
    bn = collection.first().bandNames().get(0)#.getInfo()
    
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
    bn = collection.first().bandNames().get(0)#.getInfo()
    
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

def maskAngGT40(image):
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
    
    NOTE: Final output may need to be flattened, depending on usage. See here: https://developers.google.cn/earth-engine/guides/arrays_array_images
    '''
    #Set up 3x3 kernels 
    weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
    kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

    mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
    variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

    #Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

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
    bn = collection.first().bandNames().get(0)#.getInfo()
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

#%% Quick and dirty filtering
def mask_images(ic, maskband, maskbit):
    '''
    Quick and dirty masking for MODIS data (or other simple contexts where the first QA bit is general)
    '''
    def mask_(image):
        qc = (1 << maskbit) #Bit 0 is the general quality flag
        qa = image.select(maskband) #Get the pixel QA band.
        mask = qa.bitwiseAnd(qc).eq(0)
        return image.updateMask(mask)
    return ic.map(mask_)

def std_filter(ic, threshold=3):
    '''
    Filter out very high/low values based on std from the mean
    '''
    mean = ic.reduce(ee.Reducer.mean())
    std = ic.reduce(ee.Reducer.stdDev())
    mx = mean.add(std.multiply(threshold))
    mn = mean.subtract(std.multiply(threshold))
    
    def mask_(image):
        mask = image.gt(mx)
        mask2 = image.lt(mn)
        combined_mask = (mask.Or(mask2)).Not()
        return image.updateMask(combined_mask).set('system:time_start', image.get('system:time_start'))
    
    return ic.map(mask_)

def add_constant(ic, constant=5):
    '''
    Add a constant to an ImageCollection
    '''
    def add_(image):
        return image.add(constant).set('system:time_start', image.get('system:time_start'))
    return ic.map(add_)

def check_filter_size(orig_ic, filt_ic):
    '''
    Check how much of the data is filtered by a given original and filtered collection.
    NOTE: This also works to check for MASKING -- not just removal of images  
        (e.g., masked pixels aren't included in the count)
    Returns a ratio normalized by the size of the input collection
    '''
    ct1 = ee.Image(orig_ic.count())
    ct2 = ee.Image(filt_ic.count())
    dif = ct1.subtract(ct2) #Dif between raw and filtered
    rat = ee.Image(1).subtract(dif.divide(ct1)) #Normalize by number of measurements
    return rat

def assess_fit(orig_ic, fitted_ic, bn=None):
    '''
    Check how well a model fits to the original data. Returns an RSQ value.
    '''
    #Get the input band name
    if not bn:
        img = orig_ic.first()
        bn = img.bandNames().get(0)
    
    #Get Model fit
    def distance_to_orig(image):
        resids = image.select('fitted').subtract(image.select(bn))
        ss_res = resids.pow(ee.Number(2)).rename('SumSQ')
        dist_mean = image.select(bn).subtract(mn)
        ss_tot = dist_mean.pow(ee.Number(2)).rename('DistMean')
        #rsq = 1 - (ss_res / ss_tot)
        
        return image.addBands(ss_res).addBands(ss_tot)
    
    mn = orig_ic.reduce(ee.Reducer.mean())
    distfit = fitted_ic.map(distance_to_orig)
    
    #Sum of Resids
    sum_resids = distfit.select('SumSQ').reduce(ee.Reducer.sum())
    sum_sqtot = distfit.select('DistMean').reduce(ee.Reducer.sum())
    
    #Divide the summed images and subtract from 1 for a classic RSQ value
    rsq = ee.Image(1).subtract(sum_resids.divide(sum_sqtot)).rename('RSQ')
    return rsq

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
def export_to_pandas(collection, clipper, aggregation_scale, med='median', save_std=True, mask=None):
    '''
    Takes an ImageCollection, an Earth Engine Geometry, and an aggregation scale (e.g., 30m for Landsat, 250m for MODIS, etc)
    
    Returns a pandas time series for the mean/median and standard deviation values over the 
    aggregation area. 
    
    Optionally saves those time series to a CSV file    
    
    '''
    import pandas as pd, numpy as np
    
    def createTS(image):
        date = image.get('system:time_start')
        if mask:
            image = image.updateMask(mask)
        if med == 'median':
            value = image.reduceRegion(ee.Reducer.median(), clipper, aggregation_scale)
        elif med == 'mean':
            value = image.reduceRegion(ee.Reducer.mean(), clipper, aggregation_scale)
        else:
            value = image.reduceRegion(med, clipper, aggregation_scale)
        if save_std:
            std = image.reduceRegion(ee.Reducer.stdDev(), clipper, aggregation_scale)
            ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'Mn': value, 'STD': std})
        else:
            ft = ee.Feature(None, {'system:time_start': date, 'date': ee.Date(date).format('Y/M/d'), 'Mn': value})
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
        out_vals[i] = val
        
        if save_std:
            try:
                std = list(props['STD'].values())[0]
            except:
                std = np.nan
            out_std[i] = std
        out_dates.append(pd.Timestamp(date))
    
    ser = pd.Series(out_vals, index=out_dates)
    ser = ser.sort_index()
    if save_std:
        serstd = pd.Series(out_std, index=out_dates)
        serstd = serstd.sort_index()
        return ser, serstd
    else:
        return ser
    
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
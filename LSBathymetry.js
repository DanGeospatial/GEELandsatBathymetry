//_______________________________
//          PARAMETERS
//_______________________________
var roi = geometry;
//Use lake boundaries for taking depth points for export?
//var boundaries = lakes

//The landsat 8 bands to keep
var Landsat8Bands = ['B2','B3','B4','B5','B6','B7'];
//New names for the landsat 8 bands to make them follow LS5&7 bands
var Landsat8BandRename = ['B1','B2','B3','B4','B5','B6'];

//Set start and end years
var startyear = 1991;
var finalyear = 2020;

//Set the median reducer
var medianreducer = ee.Reducer.median();

//set the chla value for depth processing
//chla value could be adjusted to specific site
var chla = 0.5;

//_______________________________
//          FUNCTIONS - LS8
//_______________________________

//this function is used to build clean water mosaic in the Google Earth Engine
//the threshold value could be revised, the current value is suggested for a common clean coral reefs waters
var cloudMaskL8 = function mask(image){
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var maskimg = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(maskimg);
};

//_______________________________
//          FUNCTIONS - LS7
//_______________________________

var cloudMaskL7 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

//_______________________________
//          FUNCTIONS - LS5
//_______________________________

var cloudMaskL5 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

//_______________________________
//          MAIN SCRIPT
//_______________________________

//Filter all Landsat collections depending on data availability between 1991 to 2020
//Run the cloud mask functions
// Load Landsat 5 Collection for oldest years
var collectionLS5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
    .filterDate('1991-01-01', '1998-12-31').filterBounds(roi).map(cloudMaskL5);
//This is the correct Landsat 7 Collection starting on 1999-01-01
var collectionLS7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                  .filterDate('1999-01-01', '2013-12-31').filterBounds(roi).map(cloudMaskL7);
//filter landsat 8 data by region and collection
var collectionLS8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
    .filterDate('2014-01-01', '2020-12-31').filterBounds(roi).map(cloudMaskL8);
    
//Filter out the ultra blue band in LS8 as it will distrupt the following calculations
var collectionLS8Fixed = collectionLS8.select(Landsat8Bands, Landsat8BandRename);

//Merge each masked collection into a single collection
var LandsatImages = collectionLS5.merge(collectionLS7).merge(collectionLS8Fixed);

//Create a list of years to map over
var years = ee.List.sequence(startyear, finalyear); 

//This function has several tasks
//
var YearDepthComp = years.map(function(year){
  
  var LandsatYear = LandsatImages.filter(ee.Filter.calendarRange(year, year, 'year'));
  
  var LandsatComp = LandsatYear.reduce(medianreducer);
  
  //Add land area filter that will remove pixels with negative NDWI values
  var mndwi = LandsatComp.normalizedDifference(['B2_median', 'B5_median']).rename ('mndwi');
  
  //Remove the negative MNDWI values from imagery
  LandsatComp = LandsatComp.where(mndwi.lt(0), ee.Number(0));
  
  //Remove images with high water turbulence by setting entire scene to 0
  //LandsatComp = LandsatComp.where();
  
  //calculate the big Rrs, rrs,and rrs*1000
  var bigrrs = LandsatComp.divide(ee.Number(31415.926));
  var rrsvec = bigrrs.divide((bigrrs.multiply(ee.Number(1.7))).add(ee.Number(0.52)));
  var rrsvec1k = rrsvec.multiply(ee.Number(1000));
  
  //calculate m0 and m1 parameters
  var m0 = ee.Number(52.073 * Math.exp(0.957*chla));
  var m1 = ee.Number(50.156 * Math.exp(0.957*chla));

  //calculate rrs vec
  var lnrrsvec = rrsvec1k.log();
  
  //calculate depth here
  var depth = ((lnrrsvec.select([1]).divide(lnrrsvec.select([2]))).multiply(m0)).subtract(m1);

  //set boundary, remove nagative value or large value in the result
  var depthA = depth.where(depth.lt(0), ee.Number(0));
  var depth_output = depthA.where(depthA.gt(20), ee.Number(20));
  
  var LandsatList = LandsatYear.aggregate_array('LANDSAT_ID');
  var Bands = LandsatComp.bandNames().size();
  
  //Return the depth map to the image collection along with several image properties
  return depth_output.set({
    'year': year,
    'image_list': LandsatList,
    'n_bands': Bands
  });
});

//Using the above function add the yearly depth maps to an image collection
var YearDepthCol = ee.ImageCollection.fromImages(YearDepthComp);

YearDepthCol = YearDepthCol.filter(ee.Filter.gt('n_bands', 0));

//Get depth values at each image pixel for all points within the lakes 
//For viewing a single year results in the EE map window
var year2015 = YearDepthCol.filter(ee.Filter.eq('year', 2019)).first();

//_______________________________
//      RESULTS AND OUTPUT
//_______________________________

var batch = require('users/fitoprincipe/geetools:batch');

batch.Download.ImageCollection.toDrive(YearDepthCol, 'Depths', {
  scale: 30,
  region: geometry //.getInfo()["coordinates"] // or geometry.getInfo()
});

//_______________________________
//      IMAGE EXPORT
//_______________________________

//depth_output is the final bathymetry
print(year2015);
//plot depth_output
Map.addLayer(year2015, {min: 0, max: 15, palette: ['00FFFF', '0000FF']});

var visParams = {
  bands: ['B3_median', 'B2_median', 'B1_median'],
  min: 0,
  max: 3000,
  gamma: 1.4,
};

//Map.addLayer(year2015, visParams);

Export.image.toDrive({
  image: year2015,
  description: 'WaterBathymetryLS8',
  scale: 30,
  region: geometry
});

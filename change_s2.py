from rasterio import open as rio_open
from xml.etree import ElementTree as ET    

def write_to_NetCDF(self, nc_outpath, compression_level, chunk_size=(1, 91, 99)):
    """ Method writing output NetCDF product.

    Keyword arguments:
    nc_outpath -- output path where NetCDF file should be stored
    compression_level -- compression level on output NetCDF file (1-9)
    chunk_size -- chunk_size
    """

    logger.info("------------START CONVERSION FROM SAFE TO NETCDF-------------")

    # Status
    utils.memory_use(self.t0)

    # Deciding a reference band
    # todo dterreng warning coming from here?
    # yes -> self.src.GetSubDatasets() ok but the gdal.Open does not work
    # add break? how to remove warning from dterr?

    def find_10m_reference_band(xml_path, safe_dir):
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for file_loc in root.findall('.//fileLocation'):
            href = file_loc.attrib['href']
            if '10m' in href:
                return safe_dir / href.split('/')[-1]
        raise ValueError("No 10m band found.")

    reference_band_path = find_10m_reference_band(self.mainXML, self.SAFE_dir)
    with rio_open(reference_band_path) as src:
        self.reference_band = src
        nx, ny = src.width, src.height


    # frequency bands
    nx = self.reference_band.RasterXSize  # number of pixels for 10m spatial resolution
    ny = self.reference_band.RasterYSize  # number of pixels for 10m spatial resolution

    # sun and view angles raster resolution
    nxa, nya = self.sunAndViewAngles[list(self.sunAndViewAngles)[0]].shape

    # output filename
    out_netcdf = (nc_outpath / self.product_id).with_suffix('.nc')

    # Sun and view angle resolutions
    nxa, nya = self.sunAndViewAngles[list(self.sunAndViewAngles)[0]].shape

    # Output filename
    ncout = (nc_outpath / self.product_id).with_suffix('.nc')


    with (netCDF4.Dataset(out_netcdf, 'w', format='NETCDF4')) as ncout:
        ncout.createDimension('time', 0)
        ncout.createDimension('x', nx)
        ncout.createDimension('y', ny)
        ncout.createDimension('raster_band_id', len(cst.s2_bands_order.keys()))
        ncout.createDimension('xa', nxa)
        ncout.createDimension('ya', nya)

        utils.create_time(ncout, self.globalAttribs["PRODUCT_START_TIME"])

        # Add projection coordinates
        ##########################################################
        # Status
        logger.info('Adding projection coordinates')
        utils.memory_use(self.t0)

        xnp, ynp = self.genLatLon(nx, ny, latlon=False)  # Assume gcps are on a regular grid

        ncx = ncout.createVariable('x', 'i4', 'x', zlib=True, complevel=compression_level)
        ncx.units = 'm'
        ncx.standard_name = 'projection_x_coordinate'
        ncx[:] = xnp

        ncy = ncout.createVariable('y', 'i4', 'y', zlib=True, complevel=compression_level)
        ncy.units = 'm'
        ncy.standard_name = 'projection_y_coordinate'
        ncy[:] = ynp

        # Add projection raster band id variable
        ##########################################################
        nc_rasterband_id = ncout.createVariable('band_id', 'i4', 'raster_band_id', zlib=True, complevel=compression_level)
        nc_rasterband_id.long_name = 'raster band id'
        nc_rasterband_id[:] = np.array(list(cst.s2_bands_order.keys()))
        nc_rasterband_id.flag_values = np.array(list(
                        cst.s2_bands_order.keys()),
                                                    dtype=np.int8)
        nc_rasterband_id.flag_meanings = ' '.join(
                        [value for value in list(cst.s2_bands_order.values())])


        # -Document
        ##########################################################
        # Status
        logger.info('Adding frequency bands layers')
        utils.memory_use(self.t0)

        if self.dterrengdata:
            # For DTERR data, gdal fails to properly do the src.GetSubDatasets()
            # so manually read the list of images created beforehand
            images = [[str(i), i.stem] for i in self.image_list_dterreng]
        else:
            images = self.src.GetSubDatasets()
        for k, v in images:
            #subdataset = gdal.Open(k)
            subdataset = xr.open_rasterio(k)
            subdataset_geotransform = subdataset.GetGeoTransform()
            # True color image (8 bit true color image)
            if ("True color image" in v) or ('TCI' in v):
                continue
            # Reflectance data for each band
            else:
                for i in range(1, subdataset.RasterCount + 1):
                    current_band = subdataset.GetRasterBand(i)
                    if self.dterrengdata:
                        band_metadata = None
                        varName = cst.s2_bands_aliases[v[-3::]]
                    else:
                        band_metadata = current_band.GetMetadata()
                        varName = band_metadata['BANDNAME']
                    if varName.startswith('B'):
                        varout = ncout.createVariable(varName, np.uint16, ('time', 'y', 'x'), fill_value=0,
                                                        zlib=True, complevel=compression_level)
                        varout.units = "1"
                        varout.grid_mapping = "UTM_projection"
                        if self.processing_level == 'Level-2A':
                            varout.standard_name = 'surface_bidirectional_reflectance'
                        else:
                            varout.standard_name = 'toa_bidirectional_reflectance'
                        varout.long_name = 'Reflectance in band %s' % varName
                        if band_metadata:
                            varout.bandwidth = band_metadata['BANDWIDTH']
                            varout.bandwidth_unit = band_metadata['BANDWIDTH_UNIT']
                            varout.wavelength = band_metadata['WAVELENGTH']
                            varout.wavelength_unit = band_metadata['WAVELENGTH_UNIT']
                            varout.solar_irradiance = band_metadata['SOLAR_IRRADIANCE']
                            varout.solar_irradiance_unit = band_metadata['SOLAR_IRRADIANCE_UNIT']
                        varout._Unsigned = "true"
                        # from DN to reflectance
                        logger.debug((varName, subdataset_geotransform))
                        if subdataset_geotransform[1] != 10:
                            current_size = current_band.XSize
                            band_measurement = scipy.ndimage.zoom(
                                input=current_band.GetVirtualMemArray(), zoom=nx / current_size,
                                order=0)
                        else:
                            band_measurement = current_band.GetVirtualMemArray()
                        varout[0, :, :] = band_measurement

        # set grid mapping
        ##########################################################
        source_crs = osr.SpatialReference()
        source_crs.ImportFromWkt(self.reference_band.GetProjection())
        nc_crs = ncout.createVariable('UTM_projection', np.int32)
        nc_crs.latitude_of_projection_origin = source_crs.GetProjParm('latitude_of_origin')
        nc_crs.proj4_string = source_crs.ExportToProj4()
        nc_crs.crs_wkt = source_crs.ExportToWkt()
        nc_crs.semi_major_axis = source_crs.GetSemiMajor()
        nc_crs.scale_factor_at_central_meridian = source_crs.GetProjParm('scale_factor')
        nc_crs.longitude_of_central_meridian = source_crs.GetProjParm('central_meridian')
        nc_crs.grid_mapping_name = source_crs.GetAttrValue('PROJECTION').lower()
        nc_crs.semi_minor_axis = source_crs.GetSemiMinor()
        nc_crs.false_easting = source_crs.GetProjParm('false_easting')
        nc_crs.false_northing = source_crs.GetProjParm('false_northing')
        nc_crs.epsg_code = source_crs.GetAttrValue('AUTHORITY', 1)
        nc_crs.crs_wkt = self.reference_band.GetProjection()

        # Add vector layers
        ##########################################################
        # Status
        logger.info('Adding masks layers')
        utils.memory_use(self.t0)
        gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

        # Add specific Level-1C or Level-2A layers
        ##########################################################
        # Status
        specific_layers_kv = {}
        utils.memory_use(self.t0)
        gdal_nc_data_types = {'Byte': 'u1', 'UInt16': 'u2'}

        if self.processing_level == 'Level-1C':
            logger.info('Adding Level-1C specific layers')
            for layer in cst.s2_l1c_layers:
                for k, v in self.imageFiles.items():
                    if layer in k:
                        specific_layers_kv[k] = cst.s2_l1c_layers[k]
                    elif layer in str(v):
                        logger.debug((layer, str(v), k))
                        specific_layers_kv[k] = cst.s2_l1c_layers[layer]


        elif self.processing_level == 'Level-2A':
            logger.info('Adding Level-2A specific layers')
            for layer in cst.s2_l2a_layers:
                for k, v in self.imageFiles.items():
                    if layer in k:
                        specific_layers_kv[k] = cst.s2_l2a_layers[k]
                    elif layer in str(v):
                        logger.debug((layer, str(v), k))
                        specific_layers_kv[k] = cst.s2_l2a_layers[layer]

        for k, v in specific_layers_kv.items():
            logger.debug((k, v))
            varName, longName = v.split(',')
            SourceDS = gdal.Open(str(self.imageFiles[k]), gdal.GA_ReadOnly)
            nb_rasterBands =  SourceDS.RasterCount

            if SourceDS.RasterCount > 1:
                logger.info("Raster data contains more than one layer")

            for i in range(1,nb_rasterBands+1):
                if nb_rasterBands>1:
                    varName =  v.split(',')[0].split()[i-1]
                    longName =  v.split(',')[1].split('-')[i-1]
                NDV = SourceDS.GetRasterBand(i).GetNoDataValue()
                xsize = SourceDS.RasterXSize
                ysize = SourceDS.RasterYSize
                GeoT = SourceDS.GetGeoTransform()
                logger.info("{}".format(GeoT))
                DataType = gdal_nc_data_types[
                    gdal.GetDataTypeName(SourceDS.GetRasterBand(i).DataType)]
                varout = ncout.createVariable(varName, DataType, ('time', 'y', 'x'), fill_value=0,
                                                zlib=True, complevel=compression_level, chunksizes=chunk_size)
                varout.grid_mapping = "UTM_projection"
                varout.long_name = longName
                if varName == "SCL":
                    varout.flag_values = np.array(list(
                        cst.s2_scene_classification_flags.values()),
                                                    dtype=np.int8)
                    varout.flag_meanings = ' '.join(
                        [key for key in list(cst.s2_scene_classification_flags.keys())])

                if GeoT[1] != 10:
                    raster_data = scipy.ndimage.zoom(input=SourceDS.GetRasterBand(i).GetVirtualMemArray(),
                                                        zoom=nx / xsize, order=0)
                else:
                    raster_data = SourceDS.GetRasterBand(i).GetVirtualMemArray()
                varout[0, :] = raster_data

        # Add sun and view angles
        ##########################################################
        # Status
        logger.info('Adding sun and view angles in native resolution')
        utils.memory_use(self.t0)

        varout_view_azimuth = ncout.createVariable('view_azimuth', np.float32, ('time','raster_band_id', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                zlib=True, complevel=compression_level)
        varout_view_azimuth.units = 'degree'
        varout_view_azimuth.long_name = 'Viewing incidence azimuth angle'
        varout_view_azimuth.comment = 'Original 22x22 pixel resolution'

        varout_view_zenith = ncout.createVariable('view_zenith', np.float32, ('time','raster_band_id', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                zlib=True, complevel=compression_level)
        varout_view_zenith.units = 'degree'
        varout_view_zenith.long_name = 'Viewing incidence zenith angle'
        varout_view_zenith.comment = 'Original 22x22 pixel resolution'

        counter = 1
        for k, v in list(self.sunAndViewAngles.items()):
            logger.debug(("Handeling %i of %i" % (counter, len(self.sunAndViewAngles))))

            if 'sun' in k:
                varout = ncout.createVariable(k, np.float32, ('time', 'ya', 'xa'), fill_value=netCDF4.default_fillvals['f4'],
                                                zlib=True, complevel=compression_level)
                varout.units = 'degree'
                varout.long_name = 'Solar %s angle' % k.split('_')[-1]
                varout.comment = 'Original 22x22 pixel resolution'
                varout[0, :, :] = v
            elif 'zenith' in k :
                band_id = k.split('_')[-1]
                varout_view_zenith[0,utils.get_key(cst.s2_bands_order,band_id), :, :] = v
            elif 'azimuth' in k :
                band_id = k.split('_')[-1]
                varout_view_azimuth[0,utils.get_key(cst.s2_bands_order,band_id), :, :] = v

            counter += 1

        # Add orbit specific data
        ##########################################################
        # Status
        logger.info('Adding satellite orbit specific data')
        utils.memory_use(self.t0)

        root = utils.xml_read(self.mainXML)
        if not self.dterrengdata:
            self.globalAttribs['orbitNumber'] = root.find('.//safe:orbitNumber',
                                                            namespaces=root.nsmap).text
        else:
            self.globalAttribs['orbitNumber'] = root.find('.//SENSING_ORBIT_NUMBER').text

        ncout.createDimension('orbit_dim', 3)
        nc_orb = ncout.createVariable('orbit_data', np.int32, ('time', 'orbit_dim'))
        rel_orb_nb = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_NUMBER']
        orb_nb = self.globalAttribs['orbitNumber']
        orb_dir = self.globalAttribs['DATATAKE_1_SENSING_ORBIT_DIRECTION']
        platform = self.globalAttribs['DATATAKE_1_SPACECRAFT_NAME']

        nc_orb.relativeOrbitNumber = rel_orb_nb
        nc_orb.orbitNumber = orb_nb
        nc_orb.orbitDirection = orb_dir
        nc_orb.platform = platform
        nc_orb.description = "Values structured as [relative orbit number, orbit number, " \
                                "platform]. platform corresponds to 0:Sentinel-2A, 1:Sentinel-2B.."
        nc_orb[0, :] = [int(rel_orb_nb), int(orb_nb), cst.platform_id[platform]]

        # Add global attributes

        logger.info('Adding global attributes')
        utils.memory_use(self.t0)

        utils.get_global_attributes(self)
        ncout.setncatts(self.globalAttribs)
        ncout.sync()

        ### Status
        logger.info('Finished.')
        utils.memory_use(self.t0)

        #write_to_geotiff(out_netcdf.is_file())

    return out_netcdf.is_file()
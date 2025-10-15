# script to convert geotiff geoid separation files to NTv2 gsb/asc or WINTER files
import datetime
import pathlib
import geodepy.geodesy as gg
import geodepy.convert as gc
from osgeo import gdal
import math
import struct


# ----------------------------------------------------------------------
# user input and options
# ----------------------------------------------------------------------
# set output type: 'asc', 'gsb', or 'winter'
output_type = 'gsb'

# set directories for geotiff and output files
tif_dir = pathlib.Path(r'')
output_dir = tif_dir / output_type

# set input geotiff file for geoid separations
tif_nval_file = tif_dir / ''

# optional - set deflection of vertical files
tif_dov_pm_file = tif_dir / ''
tif_dov_pv_file = tif_dir / ''

# prints a geotiff metadata summary then exits
print_tif_metadata = False

# set geoid grid metadata
num_orec = 11               # Number of header records in the overview block.
num_srec = 11               # Number of header records in each sub grid block.
num_file = 1                # Number of sub grids contained in the geoid grid file
gs_type = 'SECONDS '        # The units of the grid nodes
version = '1.0.0.0 '        # The geod grid file version
system_f = 'GDA2020 '       # The “from” reference ellipsoid (or geodetic datum with well known ellipsoid)
system_t = 'AHD     '       # The “to” height system or vertical datum.
major_f = gg.grs80.semimaj  # The ellipsoid semi–major axis of the “from” system.
major_t = gg.grs80.semimaj  # The ellipsoid semi–major axis of the “to” system.
minor_f = gg.grs80.semimin  # The ellipsoid semi–minor axis of the “from” system.
minor_t = gg.grs80.semimin  # The ellipsoid semi–minor axis of the “to” system.
sub_name = 'AUSGEOID'       # The name of this particular sub grid.
parent = 'NONE    '         # The parent sub grid name. 'NONE' if this is the parent grid.
created = '18062024'        # Grid file creation date
updated = '18062024'        # Grid file modification date
null_node_val = -999.0      # float value assigned to null data in geotiff


# ----------------------------------------------------------------------
# intial processing
# ----------------------------------------------------------------------

# enable gdal exceptions
gdal.UseExceptions()

# use dov files if they are available
read_dov = True if tif_dov_pv_file.is_file() and tif_dov_pv_file.is_file() else False

# intialise the output file variable and check the directory exists
if output_type == 'asc':
    out_file = output_dir / f'{tif_nval_file.stem}.asc'
elif output_type == 'gsb':
    out_file = output_dir / f'{tif_nval_file.stem}.gsb'
elif output_type == 'winter':
    out_file = output_dir / f'{tif_nval_file.stem}_win.dat'

# check directory exists
if not output_dir.exists():
    output_dir.mkdir()

# check input file exists and has been provided
if tif_nval_file.exists() is False:
    print(f'input file [{tif_nval_file.name}] not found')
    exit()
if tif_nval_file.is_dir():
    print('input file has not been provided')
    exit()

# ----------------------------------------------------------------------
# functions
# ----------------------------------------------------------------------

def check_rounding(dms_angle, precision):
    """
    Function to check if GeodePy DMS angle is being rounded to 60 arc seconds at a given precision and adjust the
    values within accordingly.
    :param dms_angle: angle of latitude/longitude
    :type dms_angle: GeodePy.angles.DMSAngle
    :param precision: decimal places of relevance
    :type precision: int
    :return:
    :rtype:
    """
    if abs(round(dms_angle.second, precision) - 60.0) < 10**-precision:
        dms_angle.second = 0.000
        dms_angle.minute += 1
        if dms_angle.minute >= 60:
            dms_angle.minute -= 60
            dms_angle.degree += 1
            if dms_angle.degree >=360:
                dms_angle.degree -= 360

def datetime_now_str():
    """
    to return the current date and time as string showing 1 decimal place
    :return: current date and time
    :rtype: str
    """

    now = datetime.datetime.now()
    ds = now.second + now.microsecond / 1000000
    dt_str = now.strftime('%Y-%d-%m %H:%M:') + f'{ds:04.1f}'

    return dt_str


def check_nan(in_val, nan_val):
    """
    function to check if a value is NaN and return the null value if true else return input value
    :param in_val: value to test
    :type in_val: float
    :param nan_val: value to return if input if NaN
    :type nan_val: float
    :return:
    :rtype: float
    """
    if math.isnan(in_val):
        return nan_val
    else:
        return in_val


# ----------------------------------------------------------------------
# begin the script
# ----------------------------------------------------------------------

print(f'{datetime_now_str()} opening geotiff files')

# open files
gt_n = gdal.Open(str(tif_nval_file))
if read_dov:
    gt_pm = gdal.Open(str(tif_dov_pm_file))
    gt_pv = gdal.Open(str(tif_dov_pv_file))

# ----------------------------------------------------------------------
# validate geotiff grids of similar size/interval
# ----------------------------------------------------------------------
print(f'{datetime_now_str()} validating grid composition')
gt_n_extents = gt_n.GetGeoTransform()
if read_dov:
    gt_pm_extents = gt_pm.GetGeoTransform()
    gt_pv_extents = gt_pv.GetGeoTransform()

    match_precision = 0.00000000001
    grid_mismatch = []
    for i in gt_n_extents:
        if abs(i - gt_pm_extents[gt_n_extents.index(i)]) > match_precision:
            mismatch_str = f'{tif_nval_file.name}/{tif_dov_pm_file.name} index: {gt_n_extents.index(i)}: ' \
            f'{i}/{gt_pm_extents[gt_n_extents.index(i)]}'
            grid_mismatch.append(mismatch_str)
        if abs(i - gt_pv_extents[gt_n_extents.index(i)]) > match_precision:
            mismatch_str = f'{tif_nval_file.name}/{tif_dov_pv_file.name} index: {gt_n_extents.index(i)}: ' \
            f'{i}/{gt_pv_extents[gt_n_extents.index(i)]}'
            grid_mismatch.append(mismatch_str)

    if grid_mismatch:
        print(' *** mismatch in geotiff grid composition')
        for e in grid_mismatch:
            print(f'     - {e}\n')
        exit()

# ----------------------------------------------------------------------
# compute grid extents
# ----------------------------------------------------------------------

x_ul, x_int, x_rot, y_ul, y_rot, y_int = gt_n_extents

gs_count = gt_n.RasterXSize * gt_n.RasterYSize
lat_inc = round(abs(y_int) * 3600, 0)
lon_inc = round(x_int * 3600, 0)
n_lat = (round(y_ul, 0) * 3600)  # convert to arc seconds
s_lat = (y_ul * 3600 - ((gt_n.RasterYSize - 1) * lat_inc))
e_lon = (-x_ul * 3600 - ((gt_n.RasterXSize - 1) * lon_inc))
w_lon = (-x_ul * 3600)

if print_tif_metadata:
    print(f'file: {tif_nval_file.name}')
    print(f'Projection info: {gt_n.GetProjection()}')
    print(f'Raster X size:   {gt_n.RasterXSize}')
    print(f'Raster Y size:   {gt_n.RasterYSize}')
    print(f'grid node count: {gt_n.RasterXSize * gt_n.RasterYSize}')
    print(f'Upper left X:    {x_ul}')
    print(f'Upper left Y:    {y_ul}')
    print(f'X interval:      {x_int}')
    print(f'Y interval:      {y_int}')
    print(f'Raster count:    {gt_n.RasterCount}')
    print(f'Metadata:        {gt_n.GetMetadata()}')
    exit()

# read in the geotiff values
band_n = gt_n.GetRasterBand(1)
values_n = band_n.ReadAsArray()
if read_dov:
    band_pm = gt_pm.GetRasterBand(1)
    values_pm = band_pm.ReadAsArray()
    band_pv = gt_pv.GetRasterBand(1)
    values_pv = band_pv.ReadAsArray()


# ----------------------------------------------------------------------
# write to file
# ----------------------------------------------------------------------

if output_type == 'asc':
    print(f'{datetime_now_str()} writing asc ntv2 file')
    with open(out_file, 'w') as f:
        asc_header = (
            f'NUM_OREC{num_orec:8d}\n'
            f'NUM_SREC{num_srec:8d}\n'
            f'NUM_FILE{num_file:8d}\n'
            f'GS_TYPE {gs_type:>8s}\n'
            f'VERSION {version:>8s}\n'
            f'SYSTEM_F{system_f:>8s}\n'
            f'SYSTEM_T{system_t:>8s}\n'
            f'MAJOR_F {major_f:16.3f}\n'
            f'MINOR_F {minor_f:16.3f}\n'
            f'MAJOR_T {major_t:16.3f}\n'
            f'MINOR_T {minor_t:16.3f}\n'
        )
        f.write(asc_header)

        sub_header = (
            f'SUB_NAME{sub_name:8s}\n'
            f'PARENT  {parent:8s}\n'
            f'CREATED {created:8s}\n'
            f'UPDATED {updated:8s}\n'
            f'S_LAT   {s_lat:16.6f}\n'
            f'N_LAT   {n_lat:16.6f}\n'
            f'E_LONG  {e_lon:16.6f}\n'
            f'W_LONG  {w_lon:16.6f}\n'
            f'LAT_INC {lat_inc:16.6f}\n'
            f'LONG_INC{lon_inc:16.6f}\n'
            f'GS_COUNT{gs_count:8d}\n'
        )
        f.write(sub_header)

        for y in range(gt_n.RasterYSize - 1, -1, -1):
            for x in range(gt_n.RasterXSize - 1, -1, -1):
                n_val = check_nan(float(values_n[y][x]), null_node_val)
                if read_dov:
                    pm_val = check_nan(float(values_pm[y][x]), null_node_val)
                    pv_val = check_nan(float(values_pv[y][x]), null_node_val)
                    f.write(f'{n_val:10.6f}{pm_val:10.6f}{pv_val:10.6f}{null_node_val:10.6f}\n')
                else:
                    f.write(f'{n_val:10.6f}{null_node_val:10.6f}{null_node_val:10.6f}{null_node_val:10.6f}\n')

        f.write(f'{"END":10s}')

elif output_type == 'gsb':
    print(f'{datetime_now_str()} writing gsb ntv2 file')
    with open(out_file, 'wb') as f:
        pad_4 = struct.pack('4x')

        header_bytes = struct.pack('8s', b'NUM_OREC') + struct.pack('i', num_orec) + pad_4
        header_bytes += struct.pack('8s', b'NUM_SREC') + struct.pack('i', num_srec) + pad_4
        header_bytes += struct.pack('8s', b'NUM_FILE') + struct.pack('i', num_file) + pad_4
        header_bytes += struct.pack('8s', b'GS_TYPE ') + struct.pack('8s', bytes(f'{gs_type[:8]:8s}', 'utf-8'))
        header_bytes += struct.pack('8s', b'VERSION ') + struct.pack('8s', bytes(f'{version[:8]:8s}', 'utf-8'))
        header_bytes += struct.pack('8s', b'SYSTEM_F') + struct.pack('8s', bytes(f'{system_f[:8]:8s}', 'utf-8'))
        header_bytes += struct.pack('8s', b'SYSTEM_T') + struct.pack('8s', bytes(f'{system_t[:8]:8s}', 'utf-8'))
        header_bytes += struct.pack('8s', b'MAJOR_F ') + struct.pack('d', major_f)
        header_bytes += struct.pack('8s', b'MINOR_F ') + struct.pack('d', minor_f)
        header_bytes += struct.pack('8s', b'MAJOR_T ') + struct.pack('d', major_t)
        header_bytes += struct.pack('8s', b'MINOR_T ') + struct.pack('d', minor_t)

        f.write(header_bytes)

        # subgrid
        sub_header_bytes = struct.pack('8s', b'SUB_NAME') + struct.pack('8s', bytes(f'{sub_name[:8]:8s}', 'utf-8'))
        sub_header_bytes += struct.pack('8s', b'PARENT  ') + struct.pack('8s', bytes(f'{parent[:8]:8s}', 'utf-8'))
        sub_header_bytes += struct.pack('8s', b'CREATED ') + struct.pack('8s', bytes(f'{created[:8]:8s}', 'utf-8'))
        sub_header_bytes += struct.pack('8s', b'UPDATED ') + struct.pack('8s', bytes(f'{updated[:8]:8s}', 'utf-8'))
        sub_header_bytes += struct.pack('8s', b'S_LAT   ') + struct.pack('d', s_lat)
        sub_header_bytes += struct.pack('8s', b'N_LAT   ') + struct.pack('d', n_lat)
        sub_header_bytes += struct.pack('8s', b'E_LONG  ') + struct.pack('d', e_lon)
        sub_header_bytes += struct.pack('8s', b'W_LONG  ') + struct.pack('d', w_lon)
        sub_header_bytes += struct.pack('8s', b'LAT_INC ') + struct.pack('d', lat_inc)
        sub_header_bytes += struct.pack('8s', b'LONG_INC') + struct.pack('d', lon_inc)
        sub_header_bytes += struct.pack('8s', b'GS_COUNT') + struct.pack('i', gs_count) + pad_4

        f.write(sub_header_bytes)

        # actual data
        for y in range(gt_n.RasterYSize - 1, -1, -1):
            for x in range(gt_n.RasterXSize - 1, -1, -1):
                n_val = check_nan(float(values_n[y][x]), null_node_val)
                if math.isnan(n_val):
                    n_val = null_node_val
                data_bytes = struct.pack('f', n_val)
                if read_dov:
                    pm_val = check_nan(float(values_pm[y][x]), null_node_val)
                    pv_val = check_nan(float(values_pv[y][x]), null_node_val)

                    data_bytes += struct.pack('f', pm_val)
                    data_bytes += struct.pack('f', pv_val)
                else:
                    data_bytes += struct.pack('f', null_node_val)
                    data_bytes += struct.pack('f', null_node_val)
                data_bytes += struct.pack('f', null_node_val)

                f.write(data_bytes)

        # end of file record
        f.write(struct.pack('8s', b'END     '))
        f.write(struct.pack('d', 3.33e+32))

elif output_type == 'winter':
    print(f'{datetime_now_str()} writing winter format .dat file')
    with open(out_file, 'w') as f:
        f.write(f'{out_file.name[:20]}{"":20s}{"www.ga.gov.au":>21s}\n')

        # write grid nodes
        for y in range(0, gt_n.RasterYSize, +1):
            lat = gc.dec2dms(y_ul + (y * y_int))
            check_rounding(dms_angle=lat, precision=3)
            lat_str = f'{"N" if lat.positive else "S"}{lat.degree:2d} {lat.minute:2d} {lat.second:6.3f}'

            for x in range(0, gt_n.RasterXSize, +1):
                lon = gc.dec2dms(x_ul + (x * x_int))
                check_rounding(dms_angle=lon, precision=3)
                lon_str = f'{"E" if lon.positive else "W"}{lon.degree:3d} {lon.minute:2d} {lon.second:6.3f}'

                n_val = check_nan(float(values_n[y][x]), null_node_val)

                if read_dov:
                    pm_val = check_nan(float(values_pm[y][x]), null_node_val)
                    pv_val = check_nan(float(values_pv[y][x]), null_node_val)

                    f.write(f'GEO {n_val:8.3f} {lat_str} {lon_str} {pm_val:9.2f} {pv_val:9.2f}\n')
                else:
                    f.write(f'GEO {n_val:8.3f} {lat_str} {lon_str} {null_node_val:9.2f} {null_node_val:9.2f}\n')


print(f'{datetime_now_str()} complete')

import arcpy
import os
from arcpy import env
from arcpy import sa
import numpy as np

def FiniteDifference(cr, cc, vm, cellsize):
	# FiniteDifference the third-order finite difference method (also called Horn's method)
	# cr is the row number of the central cell.
	# cc is the column number of the central cell.
	# vm is the vector matrix that consist of (X, Y)
	# cellsize is the resolution of the raster.

    tPara = 8*cellsize
    
    # 	#see Formulas (4) and (5).
    ACV_we = ((vm[:, cr - 1, cc + 1] - vm[:, cr - 1, cc - 1]) \
          + 2*(vm[:, cr, cc + 1]- vm[:, cr, cc - 1]) + \
              (vm[:, cr + 1, cc + 1] - vm[:, cr + 1, cc - 1]))/tPara

    ACV_ns = ((vm[:, cr - 1, cc + 1] - vm[:, cr + 1, cc + 1]) \
          + 2*(vm[:, cr - 1, cc]- vm[:, cr + 1, cc]) + \
              (vm[:, cr - 1, cc - 1] - vm[:, cr + 1, cc - 1]))/tPara

    #see Formula (6)
    out = ACV_we[0] + ACV_ns[1]

    return out


input_path = arcpy.GetParameterAsText(0)#Select the input DEM data file path
output_path = arcpy.GetParameterAsText(1)#Select the output PACV data file path
input_path_List = os.path.split(input_path)
output_path_List = os.path.split(output_path)

env.overwriteOutput = True
inputfilePath = os.path.join(input_path_List[0], input_path_List[1])
outputfilePath = os.path.join(output_path_List[0], output_path_List[1].capitalize())


dem = arcpy.Raster(inputfilePath)
no_data_value = dem.noDataValue
cellsize = dem.meanCellHeight
env.outputCoordinateSystem = dem.spatialReference
lower_left_corner_point = arcpy.Point(dem.extent.XMin, dem.extent.YMin)

A = sa.Aspect(dem)#calculate the aspect by the aspect tool in ArcPro
Theta = np.float64(arcpy.RasterToNumPyArray(A, nodata_to_value = no_data_value))
row_num, col_num = Theta.shape
Alpha = np.zeros((row_num, col_num), dtype = np.float64)

# see the formula (1) in paper.
for i in range(row_num):
    for j in range(col_num):
        if (Theta[i, j] >= 0.0) and (Theta[i, j] <= 90.0):
            Alpha[i, j] = 90.0 - Theta[i, j]
        elif Theta[i, j] < 0.0:
            Alpha[i, j] = -1
        else:
            Alpha[i, j] = 450.0 - Theta[i, j]
 
# see the formula (2) in paper.
Vector_x = cellsize*np.cos(Alpha*(np.pi/180))
Vector_y = cellsize*np.sin(Alpha*(np.pi/180))
Vector_asp = np.array([Vector_x, Vector_y])


PACVmatrix = np.zeros((row_num, col_num), dtype = np.float64)

# slide window
for i in range(1,row_num-1):
    for j in range(1, col_num-1):
        PACVmatrix[i, j] = FiniteDifference(i, j , Vector_asp, cellsize)


new_raster = arcpy.NumPyArrayToRaster(PACVmatrix, lower_left_corner = lower_left_corner_point,\
x_cell_size = cellsize, y_cell_size = cellsize, value_to_nodata = no_data_value)

new_raster.save(outputfilePath)




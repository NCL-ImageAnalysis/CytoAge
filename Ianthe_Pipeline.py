#@ CommandService command

#@ File(label="Image Folder:", style="directory") InputFolder
#@ File(label="ROI Folder:", style="directory") ROIFolder
#@ File(label="Output Folder:", style="directory") OutputFolder
#@ String (choices={"Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="listBox") threshold_method
#@ Integer (label="Minimum lower bounds for thresholding (0 for automatic):", value=0) MinLowerBounds
#@ Float (label="Line Width:", style="format:#####.#####") LineWidth
#@ Integer (label="Minimum length in pixels for ridge detection:", value=0) MinLength
#@ Integer (label="Maximum number of detected ridges (0 for no limit):", value=0) MaxRidges
#@ Integer (label="Angle determination separator:") Separator
#@ Float (label="Maximum angle for extension:") Max_Angle

import os, re, traceback

from java.lang import Math
import java.lang.Exception as JavaException

from ij import IJ, ImagePlus
from ij.io import RoiDecoder
from ij.measure import ResultsTable, Measurements
from ij.plugin import ChannelSplitter
from ij.process import AutoThresholder
from ij.plugin.filter import Analyzer

from sc.fiji.analyzeSkeleton import AnalyzeSkeleton_, Point, Edge


def brightestSlice(Image):
	maxintensity = 0
	highslice = None
	for ImageSlice in range(1, Image.getDimensions()[3]+1):
		Image.setZ(ImageSlice)
		tempmax = Image.getStatistics().max
		if tempmax > maxintensity:
			maxintensity = tempmax
			highslice = ImageSlice
	return highslice

def thresholdRanges(Image, ThresholdMethod):
	processor = Image.getProcessor()
	histogram = processor.getHistogram()
	LowerBounds = AutoThresholder().getThreshold(ThresholdMethod, histogram)
	CappedHistogram = histogram[LowerBounds:]
	CappedHistogram = [0]*(256-len(CappedHistogram)) + list(CappedHistogram)
	UpperBounds = AutoThresholder().getThreshold(ThresholdMethod, CappedHistogram)
	return LowerBounds, UpperBounds

def runRidgeDetection(Image, 
					  line_width, 
					  LowerBounds, 
					  UpperBounds, 
					  min_line_length=0,
					  max_line_length=0, 
					  darkline=False, 
					  correct_position=False, 
					  estimate_width=False, 
					  extend_line=False, 
					  slope_resolution=False, 
					  get_binary=False, 
					  get_results=False
					  ):
	
	settings = str("line_width=" + str(line_width) +
		" high_contrast=" + str(UpperBounds) + 
		" low_contrast=" + str(LowerBounds) + 
		" minimum_line_length=" + str(min_line_length) +
		" maximum=" + str(max_line_length))
	if darkline:
		settings += " darkline"
	if correct_position:
		settings += " correct_position"
	if estimate_width:
		settings += " estimate_width"
	if extend_line:
		settings += " extend_line"
	if slope_resolution:
		settings += " method_for_overlap_resolution=SLOPE"
	else:
		settings += " method_for_overlap_resolution=NONE"
	if get_binary:
		settings += " make_binary"
	if get_results:
		settings += " displayresults"
	if Image.getImageStackSize() > 1:
		settings += " stack"

	IJ.run(Image, "Ridge Detection", settings)
	OutputDict = {}
	if get_results:
		# Need to clone the results tables so we can close the windows
		OutputDict["Summary"] = ResultsTable.getResultsTable("Summary").clone()
		OutputDict["Junctions"] = ResultsTable.getResultsTable("Junctions").clone()
		OutputDict["Results"] = ResultsTable.getResultsTable("Results").clone()
		# Closing the windows
		IJ.selectWindow("Summary")
		IJ.run("Close")
		IJ.selectWindow("Junctions")
		IJ.run("Close")
		IJ.selectWindow("Results")
		IJ.run("Close")
	if get_binary:
		Title = Image.getTitle() + " Detected segments"
		IJ.selectWindow(Title)
		BinaryImage = IJ.getImage()
		# Hiding the window from the user
		BinaryImage.hide()
		OutputDict["Binary"] = BinaryImage
	return OutputDict

def areaPerSlice(Image):
	arealist = []
	IJ.setThreshold(Image, 255, 255)
	for ImageSlice in range(1, Image.getDimensions()[3]+1):
		Image.setSlice(ImageSlice)
		IJ.run(Image, "Create Selection", "")
		arealist.append(Image.getStatistics(Measurements.AREA).area)
	Image.resetRoi()
	return arealist

def get_centroid(points):
	x_sum = sum(p.x for p in points)
	y_sum = sum(p.y for p in points)
	n = float(len(points))
	return Point(Math.round(x_sum / n), Math.round(y_sum / n), 0)

def distanceBetweenPoints(a, b, is_3D=False, scaling=(1,1,1)):
	"""Calculates the distance between two coordinates

	Args:
		a (Point): Start point with x and y attributes
		b (Point): End point with x and y attributes

	Returns:
		float: Distance between the two points
	"""	

	xdiff = (a.x - b.x) * scaling[0]
	ydiff = (a.y - b.y) * scaling[1]

	if is_3D == True:
		zdiff = (a.z - b.z) * scaling[2]
		Distance = Math.sqrt((xdiff*xdiff) + (ydiff*ydiff) + (zdiff*zdiff))
	else:
		Distance = Math.sqrt((xdiff*xdiff) + (ydiff*ydiff))
	return Distance

def reverse_list_if_needed(point_list, reference_point, inverse=False):
	try:
		if distanceBetweenPoints(reference_point, point_list[0]) > distanceBetweenPoints(reference_point, point_list[-1]):
			if not inverse:
				point_list.reverse()
		elif inverse:
			point_list.reverse()
	except IndexError:
		pass
	return point_list

def getAngleBetweenPoints(Point1, Point2):
	"""Gets the angle between two points in radians

	Args:
		Point1 (tuple): X and Y coordinates of first point
		Point2 (tuple): X and Y coordinates of second point

	Returns:
		float: Angle between two points in radians
	"""	
	Angle = Math.toDegrees(Math.atan2(Point2[1] - Point1[1], Point2[0] - Point1[0]))
	return Angle

def getAngleFromLine(point_list, separator, center=None, from_end=False):
	# Needed to avoid modifying the original list
	point_list = point_list[:]
	total_points = len(point_list)
	if center:
		total_points += 1
	if total_points < 2:
		return None
	if from_end:
		point_list = point_list[-separator:]
	if center is None:
		center = point_list[0]
	try:
		Angle = getAngleBetweenPoints((center.x, center.y), (point_list[separator - 1].x, point_list[separator - 1].y))
	except IndexError:
		Angle = getAngleBetweenPoints((center.x, center.y), (point_list[-1].x, point_list[-1].y))
	return Angle

def compareAngles(a, b):
	delta = ((b - a + 180.0) % 360.0) - 180.0
	return abs(delta)

def closestAngle(base_angle, angle_list):
	best_difference = float("infinity")
	best_index = None
	for i in range(0, len(angle_list)):
		try:
			angle_difference = compareAngles(angle_list[i], base_angle)
			if angle_difference < best_difference:
				best_difference = angle_difference
				best_index = i
		except TypeError:
			continue
	return best_index, best_difference

def recursiveExtendEdge(edge, 
						start_vertex, 
						point_list=None, 
						separator=5, 
						max_angle=float("infinity"), 
						vertex_list=None, 
						vertex_locations=None
						):
	# Will initialise the points for the edge and vertex lists if not provided
	if point_list is None:
		point_list = []
	if vertex_list is None:
		vertex_list = []
	if vertex_locations is None:
		vertex_locations = []
	
	# Get all branches and eliminate the current edge
	branches = start_vertex.getBranches()
	try:
		branches.remove(edge)
	except ValueError:
		pass
	# Get the points that from the first edge
	base_slabs = edge.getSlabs()
	# Get the centroid of the starting vertex
	centroid = get_centroid(start_vertex.getPoints())

	# Needed in case it encounters a loop
	if centroid in vertex_list:
		return point_list, vertex_list, vertex_locations

	vertex_list.append(centroid)
	vertex_locations.append(len(point_list))

	# Append the centroid to the point list
	point_list.append(centroid)

	# Reverse the points from the original edge if needed
	base_slabs = reverse_list_if_needed(base_slabs, centroid, inverse=True)
	temp_slabs = base_slabs + [centroid]
	# Get the angle of the original edge
	base_angle = getAngleFromLine(temp_slabs, separator, from_end=True)
	angle_list = []
	for branch in branches:
		# Skip any branches that have already been marked as tree edges
		if branch.getType() == Edge.TREE:
			continue
		# Gets the angle of each branch relative to the vertex centroid
		edge_points = branch.getSlabs()
		edge_points = reverse_list_if_needed(edge_points, centroid)
		angle_list.append(getAngleFromLine(edge_points, separator, center=centroid))
	# If there are no valid branches, return the current point list
	if len(angle_list) < 1:
		return point_list, vertex_list, vertex_locations
	# Find the closest angle to the first edge angle
	closest_index, closest_difference = closestAngle(base_angle, angle_list)
	if closest_index is None:
		return point_list, vertex_list, vertex_locations

	# Need this to avoid modifying the original list
	NoClosest = angle_list[:]
	# Creates a new list without the closest angle for comparison
	NoClosest.pop(closest_index)
	# Gets the closest angle of the edge selected with the closest angle
	other_index, other_difference = closestAngle(angle_list[closest_index], NoClosest)
	
	# If the closest angle is less than the other angle and within the max angle, extend the edge
	if other_difference > closest_difference < max_angle:
		# Sets the next edge for the next layer of the recursion
		next_edge = branches[closest_index]
		# Marks the edge as a tree edge to avoid revisiting
		next_edge.setType(Edge.TREE)
		# Gets the new vertex to continue extending from
		new_vertex = next_edge.getOppositeVertex(start_vertex)
		# Adds the points from the next edge to the point list
		slab_points = next_edge.getSlabs()
		slab_points = reverse_list_if_needed(slab_points, centroid)
		point_list += slab_points
		# Recalls the function with the new edge and vertex
		return recursiveExtendEdge(next_edge, 
							 	   new_vertex, 
								   point_list=point_list, 
								   separator=separator, 
								   max_angle=max_angle,
								   vertex_list=vertex_list,
								   vertex_locations=vertex_locations)
	# If no valid edge is found, return the current point list
	else:
		return point_list, vertex_list, vertex_locations

def extendLineFromEdges(edge, separator=5, max_angle=float("infinity")):
	if edge.getType() == Edge.TREE:
		return
	edge.setType(Edge.TREE)
	initial_slabs = edge.getSlabs()
	v1 = edge.getV1()
	v2 = edge.getV2()
	v1_points, v1_vertices, v1_vertex_indecies = recursiveExtendEdge(edge,
																  	 v1, 
																	 separator=separator, 
																	 max_angle=max_angle 
																	)
	initial_slabs = reverse_list_if_needed(initial_slabs, v1_points[-1])
	offset = len(v1_points) + len(initial_slabs)
	v2_points, v2_vertices, v2_vertex_indecies = recursiveExtendEdge(edge,
																  	 v2, 
																	 separator=separator, 
																	 max_angle=max_angle
																	)
	both_vertices = v1_vertices + v2_vertices
	v2_vertex_indecies = [index + offset for index in v2_vertex_indecies]
	both_vertex_indecies = v1_vertex_indecies + v2_vertex_indecies
	point_list = v1_points + initial_slabs + v2_points
	return point_list, both_vertices, both_vertex_indecies

def getLengthOfPointList(point_list, is_3D=False, scaling=(1,1,1)):
	length = 0
	for ii in range(len(point_list) - 1):
		length += distanceBetweenPoints(point_list[ii], point_list[ii + 1], is_3D=is_3D, scaling=scaling)
	return length

def getCellArea(Image, scaling=(1, 1, 1)):
	CellImage = Image.duplicate()
	IJ.run(CellImage, "Gaussian Blur 3D...", "x=10 y=10 z=" + str(10 * (scaling[0]/scaling[2])))
	BrightSlice = brightestSlice(CellImage)
	CellImage.setSlice(BrightSlice)
	CellImage.resetRoi()
	CellImage.setAutoThreshold("Percentile dark")
	IJ.run(CellImage, "Convert to Mask", "method=Percentile background=Dark black")
	IJ.run(CellImage, "Fill Holes", "stack")
	cell_areas = areaPerSlice(CellImage)
	CellImage.close()
	return cell_areas

def total_3D_area(area_list, scaling=(1,1,1)):
	total_area = 0
	for area in area_list:
		total_area += area * scaling[2]
	return total_area

def dict_2_table(dictionary, key_order=None):
	rt = ResultsTable()
	n = 0
	if key_order is None:
		key_order = list(dictionary.keys())
		key_order.sort()
	for col in key_order:
		if len(dictionary[col]) > n:
			n = len(dictionary[col])
	for row in range(n):
		for col in key_order:
			try:
				rt.addValue(col, dictionary[col][row])
			except IndexError:
				rt.addValue(col, "NaN")
		if row < n-1:
			rt.incrementCounter()
	return rt

def table_2_dict(rt):
	dictionary = {}
	i = 0
	while rt.columnExists(i):
		col = rt.getColumnHeading(i)
		dictionary[col] = list(rt.getColumnAsStrings(col))
		i += 1
	return dictionary

def getRoiMeasurements(SampleRoi, Image, Measurement_Options):
	"""Gets the given measurements of the provided Roi for the given image

	Args:
		SampleRoi (ij.gui.Roi): Roi to be analysed
		Image (ij.ImagePlus): Image to be analysed
		Measurement_Options ([str]) or ([ij.measure.Measurements]): Measurements to be taken in the form of either strings of the column headings or ij.measure.Measurements integers

	Returns:
		[float]: Dictionary of measurements with column headings as titles
	"""	

	# This dictionary converts Measurement_Options to the corresponding column names in the results table
	Measurement_Dict = {
		Measurements.ADD_TO_OVERLAY: [None],
		Measurements.AREA: ['Area'],
		Measurements.AREA_FRACTION: ['%Area'],
		Measurements.CENTER_OF_MASS: ['XM', 'YM'],
		Measurements.CENTROID: ['X', 'Y'],
		Measurements.CIRCULARITY: ['Circ.', 'AR', 'Round', 'Solidity'],
		Measurements.ELLIPSE: ['Major', 'Minor', 'Angle'],
		Measurements.FERET: ['Feret', 'FeretX', 'FeretY', 'FeretAngle', 'MinFeret'],
		Measurements.INTEGRATED_DENSITY: ['IntDen'],
		Measurements.INVERT_Y: [None],
		Measurements.KURTOSIS: ['Kurt'],
		Measurements.LABELS: ['Label'],
		Measurements.LIMIT: [None],
		Measurements.MAX_STANDARDS: [None],
		Measurements.MEAN: ['Mean'],
		Measurements.MEDIAN: ['Median'],
		Measurements.MIN_MAX: ['Min', 'Max'],
		Measurements.MODE: ['Mode'],
		Measurements.NaN_EMPTY_CELLS: [None],
		Measurements.PERIMETER: ['Perim.'],
		Measurements.RECT: ['ROI_X', 'ROI_Y', 'ROI_Width', 'ROI_Height'],
		Measurements.SCIENTIFIC_NOTATION: [None],
		Measurements.SHAPE_DESCRIPTORS: ['Circ.', 'AR', 'Round', 'Solidity'],
		Measurements.SKEWNESS: ['Skew'],
		Measurements.SLICE: [None],
		Measurements.STACK_POSITION: [None],
		Measurements.STD_DEV: ['StdDev']
	}

	# Initialises a new empty results table
	RTable = ResultsTable()
	# Initialises an Analyzer object using 
	# the image and the empty results table
	try:
		# If input list is of ij.measure.Measurements will use those measurements for the analyzer
		Measurement_int = sum(Measurement_Options)
		An = Analyzer(Image, Measurement_int, RTable)
	except TypeError:
		# Otherwise will just use global measurement options
		Measurement_int = None
		An = Analyzer(Image, RTable)
	# Selects the roi on the image
	Image.setRoi(SampleRoi)
	# Takes the measurements
	An.measure()
	# If the measurements were not specified
	# will use input column headings
	if Measurement_int == None:
		Output_List = Measurement_Options
	# Otherwise will get measurement options from dictionary
	else:
		Output_List = []
		for Option in Measurement_Options:
			Output_List += Measurement_Dict[Option]
	# Takes the desired results from the results table and adds to the dictionary
	OutputDict = {}
	for Option in Output_List:
		if Option != None:
			OutputDict[Option] = RTable.getValue(Option, 0)
	# Clears the results table
	RTable.reset()
	# Clears the roi from the image
	Image.resetRoi()
	return OutputDict

def nan_zero_divide(a, b, float_result=True):
	if float_result:
		a = float(a)
		b = float(b)
	try:
		result = a / b
	except ZeroDivisionError:
		result = "NaN"
	return result

def main(InputPath, 
		 OutputPath,
		 roi_path,
		 thresh_method,
		 min_lower_bounds,
		 line_width,
		 min_length,
		 max_ridges,
		 separator,
		 max_angle):
	regexpattern = re.compile(r"\.tif{1,2}$|\.czi$", re.IGNORECASE)
	all_cells_keys = ["Filename",
				   	  "Channel",
					  "Cell Volume (um^3)",
					  "Cell Orientation",
					  "Skeleton Structure Size (um^3)",
					  "Skeleton Density",
					  "Filament Count",
					  "Branch Count",
					  "Branch Ratio",
					  "Branch Point Density (/um^3)",
					  "Average Filament Length 2D (um)",
					  "Average Filament Length 3D (um)",
					  "Average Filament Width (um)",
					  "Thresholding Method",
					  "Lower Threshold",
					  "Upper Threshold"
					  ]
	
	if os.path.exists(os.path.join(OutputPath, "Per_Cell_Measurements.csv")):
		IJ.log("Loading existing Per Cell Measurements CSV...")
		all_cells_rt = ResultsTable.open(os.path.join(OutputPath, "Per_Cell_Measurements.csv"))
		all_cells_dict = table_2_dict(all_cells_rt)
	else:
		all_cells_dict = {key: [] for key in all_cells_keys}

	coords_keys = ["Filename",
				   "Channel",
				   "Filament ID",
				   "Coordinates (um)",
				   "Verticies",
				   "Vertex Locations"
				   ]
	
	if os.path.exists(os.path.join(OutputPath, "Per_Filament_Coordinates.csv")):
		IJ.log("Loading existing Per Filament Coordinates CSV...")
		coords_rt = ResultsTable.open(os.path.join(OutputPath, "Per_Filament_Coordinates.csv"))
		coords_dict = table_2_dict(coords_rt)
	else:
		coords_dict = {key: [] for key in coords_keys}

	if os.path.exists(os.path.join(OutputPath, "Per_Filament_Widths.csv")):
		IJ.log("Loading existing Per Filament Widths CSV...")
		width_rt = ResultsTable.open(os.path.join(OutputPath, "Per_Filament_Widths.csv"))
		width_dict = table_2_dict(width_rt)
		width_keys = sorted(list(width_dict.keys()))
	else:
		width_dict = {}
		width_keys = []

	for root, dirs, files in os.walk(InputPath):
		files = [f for f in files if regexpattern.search(f)]
		for f in files:
			try:
				# Data and metadata loading--------------------------------------------------------v
				# Load image
				filepath = os.path.join(root, f)
				if filepath in all_cells_dict["Filename"]:
					IJ.log("File already processed, skipping: " + filepath)
					continue

				RawImage = ImagePlus(filepath)
				# Get image scaling
				calibration = RawImage.getCalibration()
				scaling = (calibration.pixelWidth, calibration.pixelHeight, calibration.pixelDepth)
				# ROI loading
				roi_filepath = ".".join(filepath.split(".")[:-1]) + ".roi"
				roi_filepath = roi_filepath.replace(InputPath, roi_path)
				if not os.path.exists(roi_filepath):
					IJ.log("ROI file does not exist: " + roi_filepath)
					continue
				roi_decoder = RoiDecoder(roi_filepath)
				cell_roi = roi_decoder.getRoi()
				RawImage.setRoi(cell_roi)
				IJ.run(RawImage, "Clear Outside", "stack")
				#----------------------------------------------------------------------------------^
				channel_images = ChannelSplitter.split(RawImage)
				RawImage.close()
				for channel, Image in enumerate(channel_images):
					try:
						IJ.log("Processing File: " + filepath + " Channel: " + str(channel))
						# Getting cell measurements
						cell_areas = getCellArea(Image, scaling=scaling)
						cell_feret = getRoiMeasurements(cell_roi, Image, [Measurements.FERET])
						cell_3D_area = total_3D_area(cell_areas, scaling=scaling)

						# Ridge detection--------------------------------------------------------------v
						# Preprocessing
						# Sets the contrast to have 0.01% saturation of the pixels (across whole stack)
						Image.setRoi(cell_roi)
						IJ.run(Image, "Enhance Contrast...", "saturated=0.01 use")
						# Converts the image to 8-bit
						IJ.run(Image, "8-bit", "")
						BrightSlice = brightestSlice(Image)
						Image.setSlice(BrightSlice)
						# Determines threshold ranges for ridge detection
						LowerBounds, UpperBounds = thresholdRanges(Image, thresh_method)
						if LowerBounds < min_lower_bounds:
							IJ.log("Lower bound for thresholding is " + str(LowerBounds) + ", likely indicating poor image quality or inappropriate threshold option. Skipping this channel.")
							continue
						IJ.log("Performing ridge detection with lower bound: " + str(LowerBounds) + " and upper bound: " + str(UpperBounds))
						# Runs the ridge detection
						RidgeOutput = runRidgeDetection(Image, 
											line_width,
											LowerBounds,
											UpperBounds,
											correct_position=True,
											estimate_width=True,
											#    slope_resolution=True,
											get_binary=True,
											get_results=True,
											min_line_length=min_length
											)
						# Getting ridge detection measurements
						ridge_area_list = areaPerSlice(RidgeOutput["Binary"])
						ridge_3D_area = total_3D_area(ridge_area_list, scaling=scaling)
						skeleton_density = nan_zero_divide(ridge_3D_area, cell_3D_area, float_result=True)
						ridge_widths = RidgeOutput["Summary"].getColumn("Mean line width")
						# Skip if too many ridges detected
						if len(ridge_widths) > max_ridges > 0:
							IJ.log(str(len(ridge_widths)) + " ridges detected, which is more than the maximum of " + str(max_ridges) + ". Skipping this channel.")
							continue
						avg_ridge_width = nan_zero_divide(sum(ridge_widths), len(ridge_widths), float_result=True)
						#------------------------------------------------------------------------------^
						IJ.log("Running analyze skeleton...")
						IJ.run(RidgeOutput["Binary"], "Skeletonize (2D/3D)", "")
						AS = AnalyzeSkeleton_()
						AS.setup("", RidgeOutput["Binary"])
						skellyres = AS.run(AnalyzeSkeleton_.SHORTEST_BRANCH, # Prune cycles index
										False, # Prune ends
										False, # Calculate shortest path
										RidgeOutput["Binary"],  # ImagePlus skeleton
										True, # Silent mode
										False # Verbose mode
						)
						IJ.log("Analyzing skeleton output...")
						graph = skellyres.getGraph()
						edge_list = []
						for g in graph:
							for e in g.getEdges():
								e.setType(Edge.UNDEFINED)
								edge_list.append(e)
						length_list_2D = []
						length_list_3D = []
						coords_list = []
						vertex_list = []
						v_indeces_list = []
						for ed in edge_list:
							output_tuple = extendLineFromEdges(ed, separator=separator, max_angle=max_angle)
							if output_tuple is None:
								continue
							point_list, vertices, vertex_indecies = output_tuple
							length_list_2D.append(getLengthOfPointList(point_list, is_3D=False, scaling=scaling))
							length_list_3D.append(getLengthOfPointList(point_list, is_3D=True, scaling=scaling))
							coords_list.append(str([(p.x * scaling[0], p.y * scaling[1], p.z * scaling[2]) for p in point_list]))
							vertex_list.append(str([(v.x * scaling[0], v.y * scaling[1], v.z * scaling[2]) for v in vertices]))
							v_indeces_list.append(str(vertex_indecies))
						filament_count = len(length_list_3D)
						avg_filament_length_2D = nan_zero_divide(sum(length_list_2D), filament_count, float_result=True)
						avg_filament_length_3D = nan_zero_divide(sum(length_list_3D), filament_count, float_result=True)
						num_branches = sum(skellyres.getBranches())
						branch_ratio = nan_zero_divide(num_branches, filament_count, float_result=True)
						branch_density = nan_zero_divide(num_branches, cell_3D_area, float_result=True)

						# Data Outputs-------------------------------------------------------------------v
						# All Cells Dictionary----------------------------------------------------------v
						all_cells_dict["Filename"].append(filepath)
						all_cells_dict["Channel"].append(channel)
						all_cells_dict["Cell Volume (um^3)"].append(cell_3D_area)
						all_cells_dict["Cell Orientation"].append(cell_feret["FeretAngle"])
						all_cells_dict["Skeleton Structure Size (um^3)"].append(ridge_3D_area)
						all_cells_dict["Skeleton Density"].append(skeleton_density)
						all_cells_dict["Filament Count"].append(filament_count)
						all_cells_dict["Average Filament Length 2D (um)"].append(avg_filament_length_2D)
						all_cells_dict["Average Filament Length 3D (um)"].append(avg_filament_length_3D)
						all_cells_dict["Average Filament Width (um)"].append(avg_ridge_width)
						all_cells_dict["Branch Count"].append(num_branches)
						all_cells_dict["Branch Ratio"].append(branch_ratio)
						all_cells_dict["Branch Point Density (/um^3)"].append(branch_density)
						all_cells_dict["Thresholding Method"].append(thresh_method)
						all_cells_dict["Lower Threshold"].append(LowerBounds)
						all_cells_dict["Upper Threshold"].append(UpperBounds)
						# Convert to a results table
						all_cells_rt = dict_2_table(all_cells_dict, key_order=all_cells_keys)
						# Save as CSV
						all_cells_rt.save(os.path.join(OutputPath, "Per_Cell_Measurements.csv"))
						#-------------------------------------------------------------------------------^
						# Per Filament Dictionary-------------------------------------------------------v
						for filament_id in range(len(coords_list)):
							coords_dict["Filename"].append(filepath)
							coords_dict["Channel"].append(channel)
							coords_dict["Filament ID"].append(filament_id)
							coords_dict["Coordinates (um)"].append(coords_list[filament_id])
							coords_dict["Verticies"].append(vertex_list[filament_id])
							coords_dict["Vertex Locations"].append(v_indeces_list[filament_id])
						# Convert to results tables
						coords_rt = dict_2_table(coords_dict, key_order=coords_keys)
						width_dict[filepath + "-Channel-" + str(channel)] = ridge_widths
						width_keys.append(filepath + "-Channel-" + str(channel))
						width_rt = dict_2_table(width_dict, key_order=width_keys)
						# Save as CSVs
						coords_rt.save(os.path.join(OutputPath, "Per_Filament_Coordinates.csv"))
						width_rt.save(os.path.join(OutputPath, "Per_Filament_Widths.csv"))
					except (Exception, JavaException) as e:
						IJ.log("Error processing file: " + filepath + " Channel: " + str(channel) + "\n" + traceback.format_exc())

			except (Exception, JavaException) as e:
				IJ.log("Error processing file: " + filepath + "\n" + traceback.format_exc())

inpath = InputFolder.getPath()
outpath = OutputFolder.getPath()
roipath = ROIFolder.getPath()


if __name__ == "__main__":
	main(inpath,
		 outpath,
		 roipath,
		 threshold_method,
		 MinLowerBounds,
		 LineWidth,
		 MinLength,
		 MaxRidges,
		 Separator,
		 Max_Angle
	)
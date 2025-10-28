#@ Float (label="Line Width:", style="format:#####.#####") LineWidth
#@ Integer (label="Minimum length in pixels for ridge detection:", value=0) MinLength
#@ String (choices={"Test All", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="listBox") thresh_method


from ij import IJ
from ij.process import AutoThresholder
from ij.plugin import ImagesToStack

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
					  slope_resolution=False
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
	if Image.getImageStackSize() > 1:
		settings += " stack"

	IJ.run(Image, "Ridge Detection", settings)

if thresh_method == "Test All":
	ThreshList = ["Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"]
else:
	ThreshList = [thresh_method]
RawImage = IJ.getImage()
Image = RawImage.duplicate()
IJ.run(Image, "Enhance Contrast...", "saturated=0.01 use")
# Converts the image to 8-bit
IJ.run(Image, "8-bit", "")
BrightSlice = brightestSlice(Image)
Image.setSlice(BrightSlice)
stack_list = []
for threshold_method in ThreshList:
	ThisImage = Image.duplicate()
	ThisImage.setTitle(threshold_method)
	# Determines threshold ranges for ridge detection
	LowerBounds, UpperBounds = thresholdRanges(ThisImage, threshold_method)
	IJ.log("Determined lower bound: " + str(LowerBounds) + " and upper bound: " + str(UpperBounds) + " for thresholding using " + threshold_method + " method.")
	RidgeOutput = runRidgeDetection(ThisImage, 
							LineWidth,
							LowerBounds,
							UpperBounds,
							correct_position=True,
							estimate_width=True,
							min_line_length=MinLength
							)
	stack_list.append(ThisImage)
if len(stack_list) > 1:
	stacked = ImagesToStack.run(stack_list)
	stacked.show()
else:
	stack_list[0].show()
Image.close()
from ij import IJ
from ij.process import AutoThresholder

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

ThreshList = ["Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"]
RawImage = IJ.getImage()
Image = RawImage.duplicate()
IJ.run(Image, "Enhance Contrast...", "saturated=0.01 use")
# Converts the image to 8-bit
IJ.run(Image, "8-bit", "")
BrightSlice = brightestSlice(Image)
Image.setSlice(BrightSlice)
Image.show()
for threshold_method in ThreshList:
	# Determines threshold ranges for ridge detection
	LowerBounds, UpperBounds = thresholdRanges(Image, threshold_method)
	IJ.log("Determined lower bound: " + str(LowerBounds) + " and upper bound: " + str(UpperBounds) + " for thresholding using " + threshold_method + " method.")
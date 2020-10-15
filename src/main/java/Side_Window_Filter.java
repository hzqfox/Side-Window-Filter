
import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.Arrays;

import ij.*;
import ij.process.*;
import ij.gui.*;
//import ij.plugin.*;
import ij.plugin.filter.*;

/** This plugin implements the side window filtering method on top of: 
 * {"Box", "Gaussian", "Median", "Bilateral", "Guided"} filters
 * 
 * The side window filtering is based on Yuanhao's work
 * , and re-written from the original MATLAB code:
 * 		https://github.com/YuanhaoGong/SideWindowFilter
 * 		CVPR 2019 oral, #5176
 * 
 * Author:	Ziqiang Huang {github.com/hzqfox}
 * Date:	2020.01.02
 * Version:	1.0.0
 *
 */

import ij.*;
import ij.process.*;
import ij.gui.GenericDialog;
import ij.util.ThreadUtil;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.gui.*;
import java.util.concurrent.atomic.AtomicInteger;

/*
 * This plugin implements most of the 3D filters in the Process/Filters submenu.
 * @author Thomas Boudier
 */

public class Side_Window_Filter implements PlugIn {

	public final static int MEAN = 10, MEDIAN = 11, MIN = 12, MAX = 13, VAR = 14, MAXLOCAL = 15;
	private static float xradius = 2, yradius = 2, zradius = 2;

	public void run(String arg) {
		String name = null;
		int filter = 0;
		if (arg.equals("mean")) {
			name = "3D Mean";
			filter = MEAN;
		} else if (arg.equals("median")) {
			name = "3D Median";
			filter = MEDIAN;
		} else if (arg.equals("min")) {
			name = "3D Minimum";
			filter = MIN;
		} else if (arg.equals("max")) {
			name = "3D Maximum";
			filter = MAX;
		} else if (arg.equals("var")) {
			name = "3D Variance";
			filter = VAR;
		} else
			return;
		ImagePlus imp = IJ.getImage();
		if (imp.isComposite() && imp.getNChannels() == imp.getStackSize()) {
			IJ.error(name, "Composite color images not supported");
			return;
		}
		if (!showDialog(name))
			return;
		imp.startTiming();
		run(imp, filter, xradius, yradius, zradius);
		IJ.showTime(imp, imp.getStartTime(), "", imp.getStackSize());
	}

	private boolean showDialog(String name) {
		GenericDialog gd = new GenericDialog(name);
		gd.addNumericField("X radius:", xradius, 1);
		gd.addNumericField("Y radius:", yradius, 1);
		gd.addNumericField("Z radius:", zradius, 1);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return false;
		}
		xradius = (float) gd.getNextNumber();
		yradius = (float) gd.getNextNumber();
		zradius = (float) gd.getNextNumber();
		return true;
	}

	private void run(ImagePlus imp, int filter, float radX, float radY, float radZ) {
		if (imp.isHyperStack()) {
			filterHyperstack(imp, filter, radX, radY, radZ);
			return;
		}
		ImageStack res = filter(imp.getStack(), filter, radX, radY, radZ);
		imp.setStack(res);
	}

	public static ImageStack filter(ImageStack stackorig, int filter, float vx, float vy, float vz) {

		if (stackorig.getBitDepth() == 24)
			return filterRGB(stackorig, filter, vx, vy, vz);

		// get stack info
		final ImageStack stack = stackorig;
		final float voisx = vx;
		final float voisy = vy;
		final float voisz = vz;
		final int width = stack.getWidth();
		final int height = stack.getHeight();
		final int depth = stack.size();
		ImageStack res = null;

		if ((filter == MEAN) || (filter == MEDIAN) || (filter == MIN) || (filter == MAX) || (filter == VAR)) {
			if (filter == VAR)
				res = ImageStack.create(width, height, depth, 32);
			else
				res = ImageStack.create(width, height, depth, stackorig.getBitDepth());
			IJ.showStatus("3D filtering...");
			// PARALLEL
			final ImageStack out = res;
			final AtomicInteger ai = new AtomicInteger(0);
			final int n_cpus = Prefs.getThreads();

			final int f = filter;
			final int dec = (int) Math.ceil((double) stack.size() / (double) n_cpus);
			Thread[] threads = ThreadUtil.createThreadArray(n_cpus);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						StackProcessor processor = new StackProcessor(stack);
						for (int k = ai.getAndIncrement(); k < n_cpus; k = ai.getAndIncrement()) {
							processor.filter3D(out, voisx, voisy, voisz, dec * k, dec * (k + 1), f);
						}
					}
				};
			}
			ThreadUtil.startAndJoin(threads);
		}
		return res;
	}

	private static void filterHyperstack(ImagePlus imp, int filter, float vx, float vy, float vz) {
		if (imp.getNDimensions() > 4) {
			IJ.error("5D hyperstacks are currently not supported");
			return;
		}
		if (imp.getNChannels() == 1) {
			ImageStack stack = filter(imp.getStack(), filter, vx, vy, vz);
			imp.setStack(stack);
			return;
		}
		ImagePlus[] channels = ChannelSplitter.split(imp);
		int n = channels.length;
		for (int i = 0; i < n; i++) {
			ImageStack stack = filter(channels[i].getStack(), filter, vx, vy, vz);
			channels[i].setStack(stack);
		}
		ImagePlus imp2 = RGBStackMerge.mergeChannels(channels, false);
		imp.setImage(imp2);
		imp.setC(1);
	}

	private static ImageStack filterRGB(ImageStack rgb_in, int filter, float vx, float vy, float vz) {
		ImageStack[] channels = ChannelSplitter.splitRGB(rgb_in, false);
		ImageStack red = filter(channels[0], filter, vx, vy, vz);
		ImageStack green = filter(channels[1], filter, vx, vy, vz);
		ImageStack blue = filter(channels[2], filter, vx, vy, vz);
		return RGBStackMerge.mergeStacks(red, green, blue, false);
	}

}

/*
 * public class Side_Window_Filter_ implements ExtendedPlugInFilter,
 * DialogListener {
 * 
 * ImagePlus impInput; ImagePlus imgThr; ColorProcessor ipTemp, ipC;
 * 
 * static int method = -1; String[] methods = {"Box", "Gaussian", "Median",
 * "Bilateral", "Guided"}; int[] idMethods = {0, 1, 2, 3, 4};
 * 
 * static int bSize = 15; double maxValue = 255; static double param1 = 0;
 * 
 * int[] dimensions; boolean ok = false; //whether the ok button was pressed in
 * dialog boolean cancel = false; //whether the cancel button was pressed in
 * dialog boolean output,pv; int maxB = 3;
 * 
 * 
 * public static final int BOX=0, GAU=1, MED=2, BIL=3, GUI=4; private static int
 * HIGHEST_FILTER = GUI;
 * 
 * // Filter parameters private double radius; private int filterType;
 * 
 * // Remember filter parameters for the next time private static double[]
 * lastRadius = new double[HIGHEST_FILTER+1]; //separate for each filter type
 * 
 * // Further class variables private int flags = DOES_ALL | KEEP_PREVIEW;
 * //sotre flags in plugin private ImagePlus imp; private int nPasses = 1; //
 * The number of passes (color channels * stack slices) private
 * PlugInFilterRunner pfr; private boolean previewing = false;
 * 
 * 
 *//** need to keep the instance of ImagePlus */
/*
 * private ImagePlus imagePlus;
 * 
 *//** keep the original image, to restore it after the preview */
/*
 * private ImageProcessor baseImage;
 * 
 *//** Keep instance of result image */
/*
 * private ImageProcessor result;
 * 
 * 
 * 
 * private int pass; // Multithreading - related private int numThreads =
 * Prefs.getThreads(); // Current state of processing is in class variables.
 * Thus, stack parallelization must be done // ONLY with one thread for the
 * image (not using these class variables): private int highestYinCache; // the
 * highest line read into the cache so far private boolean threadWaiting; // a
 * thread waits until it may read data private boolean copyingToCache;
 * 
 *//**
	 * Setup of the PlugInFilter. Returns the flags specifying the capabilities and
	 * needs of the filter.
	 *
	 * @param arg Defines type of filter operation
	 * @param imp The ImagePlus to be processed
	 * @return Flags specifying further action of the PlugInFilterRunner
	 */
/*
 * public int setup(String arg, ImagePlus imp) {
 * 
 * // about... if (arg.equals("about")) { showAbout(); return DONE; }
 * this.impInput = imp; dimensions = impInput.getDimensions();
 * 
 * 
 * if (arg.equals("box")) filterType = BOX; else if (arg.equals("mean"))
 * filterType = BOX; else if (arg.equals("Gaussian")) filterType = GAU; else if
 * (arg.equals("median")) filterType = MED; else if (arg.equals("bilateral"))
 * filterType = BIL; else if (arg.equals("guided")) filterType = GUI; else {
 * IJ.error("RankFilters","Argument missing or undefined: "+arg); return DONE; }
 * 
 * 
 * return flags; }
 * 
 * public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
 * {
 * 
 * ipC = (ColorProcessor)impInput.getProcessor().convertToRGB();
 * 
 * 
 * //GenericDialog gd = new GenericDialog("Adaptive Threshold",
 * IJ.getInstance()); GenericDialog gd = new
 * NonBlockingGenericDialog("Side Window Filter");
 * 
 * 
 * method = method==-1 ? 0 : method; gd.addChoice("Filter Type", methods,
 * methods[method]);
 * 
 * radius = lastRadius[filterType]<=0 ? 2 : lastRadius[filterType];
 * gd.addNumericField("Radius", radius, 1, 6, "pixels");
 * 
 * int digits = imp.getType() == ImagePlus.GRAY32 ? 2 : 0;
 * 
 * gd.addPreviewCheckbox(pfr); //passing pfr makes the filter ready for preview
 * gd.addDialogListener(this); //the DialogItemChanged method will be called on
 * user input gd.showDialog(); if (gd.wasCanceled()) return DONE;
 * 
 * 
 * IJ.register(this.getClass()); //protect static class variables (filter
 * parameters) from garbage collection if (Macro.getOptions() == null) {
 * //interactive only: remember parameters entered lastRadius[filterType] =
 * radius; }
 * 
 * this.pfr = pfr; flags = IJ.setupDialog(imp, flags); //ask whether to process
 * all slices of stack (if a stack)
 * 
 * if ((flags&DOES_STACKS)!=0) { int size = imp.getWidth() * imp.getHeight();
 * 
 * }
 * 
 * return flags;
 * 
 * }
 * 
 * public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
 * 
 * method = idMethods[gd.getNextChoiceIndex()]; radius = gd.getNextNumber(); pv
 * = gd.getPreviewCheckbox().getState();
 * 
 * int maxRadius = 55; //HERE, upperbound for computation
 * 
 * if(gd.invalidNumber() || radius<0 || radius>maxRadius){ return false; }
 * 
 * return true; }
 * 
 * 
 * private void parseDialogParameters(GenericDialog gd) { // extract chosen
 * parameters this.method = idMethods[gd.getNextChoiceIndex()]; this.radius =
 * gd.getNextNumber(); this.previewing = gd.getPreviewCheckbox().getState(); }
 * 
 * 
 * 
 * public void run(ImageProcessor ip) {
 * 
 * swf(ip, radius, filterType);
 * 
 * if (IJ.escapePressed()) // interrupted by user? ip.reset();
 * 
 * pass++;
 * 
 * }
 * 
 * 
 *//**
	 * Filter an image using side window filtering method
	 * 
	 * @param ip         The ImageProcessor that would be filtered ()
	 * @param radius     Determins the kernel size (only rectangular)
	 * @param filterType May be BOX, GAUSSIAN, MEDIAN, BILATERAL, or GUIDED.
	 */
/*
 * public void swf(ImageProcessor ip, double radius, int filterType) {
 * 
 * if (Thread.currentThread().isInterrupted()) return;
 * 
 * Rectangle roi = ip.getRoi(); ImageProcessor mask = ip.getMask(); Rectangle
 * roi1 = null;
 * 
 * //int[] lineRadii = makeLineRadii(radius);
 * 
 * boolean isImagePart = (roi.width<ip.getWidth()) ||
 * (roi.height<ip.getHeight());
 * 
 * boolean[] aborted = new boolean[1]; // returns whether interrupted during
 * preview or ESC pressed
 * 
 * int[][] oldArray = ip.getIntArray(); int width = ip.getWidth(); int height =
 * ip.getHeight(); int[][] newArray = new int[width][height];
 * 
 * for (int h=0; h<height; h++) { for (int w=0; w<width; w++) { newArray =
 * Filter_Array_Computations.Filter2DArray(oldArray, radius, filterType); } }
 * 
 * ip.setIntArray(newArray);
 * 
 * 
 * //Object pixels = ip.getPixels();
 * 
 * //result = op.apply(this.baseImage, strel);
 * 
 * 
 * 
 * if (previewing) { // Fill up the values of original image with values of the
 * result for (int i = 0; i < ip.getPixelCount(); i++) { ip.setf(i,
 * result.getf(i)); } ip.resetMinAndMax(); }
 * 
 * 
 * 
 * 
 * for (int ch=0; ch<ip.getNChannels(); ch++) { int filterType1 = filterType;
 * 
 * 
 * doFiltering(ip, lineRadii, filterType1, ch, aborted);
 * 
 * 
 * if (aborted[0]) break;
 * 
 * }
 * 
 * }
 * 
 * 
 * // Filter a grayscale image or one channel of an RGB image with several
 * threads // Implementation: each thread uses the same input buffer (cache),
 * always works on the next unfiltered line // Usually, one thread reads reads
 * several lines into the cache, while the others are processing the data. //
 * 'aborted[0]' is set if the main thread has been interrupted (during preview)
 * or ESC pressed. // 'aborted' must not be a class variable because it signals
 * the other threads to stop; and this may be caused // by an interrupted
 * preview thread after the main calculation has been started. private void
 * doFiltering(final ImageProcessor ip, final int[] lineRadii, final int
 * filterType, final int colorChannel, final boolean[] aborted) { }
 * 
 * 
 *//**
	 * This method is called by ImageJ to set the number of calls to run(ip)
	 * corresponding to 100% of the progress bar
	 *//*
		 * public void setNPasses (int nPasses) { this.nPasses = nPasses; pass = 0; }
		 * 
		 * private void showProgress(double percent, boolean rgb) { int nPasses2 =
		 * rgb?nPasses*3:nPasses; percent = (double)pass/nPasses2 + percent/nPasses2;
		 * IJ.showProgress(percent); }
		 * 
		 * // About... private void showAbout() { final String title =
		 * "Side Window Filter"; final String body = "Side Window Filter,\n" +
		 * "github.com/YuanhaoGong/SideWindowFilter\n" + "\n" + "by YuanHao Gong\n" +
		 * "(github.com/YuanhaoGong)\n" + "by Ziqiang Huang\n" + "(github.com/hzqfox)" +
		 * "\n"; IJ.showMessage(title, body); }
		 * 
		 * 
		 * 
		 * }
		 */
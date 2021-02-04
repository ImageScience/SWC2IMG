package swc2img;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

import imagescience.feature.Smoother;
import imagescience.image.Image;
import imagescience.random.PoissonGenerator;
import imagescience.shape.Bounds;
import imagescience.utility.FMath;
import imagescience.utility.Formatter;
import imagescience.utility.LineReader;
import imagescience.utility.Progressor;
import imagescience.utility.Timer;

import java.awt.Button;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.StringTokenizer;

/** ImageJ plugin for neuron image synthesis. */
public final class SWC2IMG implements PlugIn, ActionListener, WindowListener {
	
	// Parameters:
	private String file = "";
	private boolean voxel = false;
	private double lateral = 1.0;
	private double axial = 1.0;
	private double background = 10;
	private double correlation = 1;
	private double snr = 10;
	private boolean rescale = false;
	private boolean log = false;
	
	// Limits:
	private final static double MIN_LATERAL = 0.001;
	private final static double MAX_LATERAL = 1000;
	private final static double MIN_AXIAL = 0.001;
	private final static double MAX_AXIAL = 1000;
	private final static double MIN_BACKGROUND = 0;
	private final static double MAX_BACKGROUND = 1000;
	private final static double MIN_CORRELATION = 0;
	private final static double MAX_CORRELATION = 1000;
	private final static double MIN_SNR = 0;
	private final static double MAX_SNR = 1000;
	
	// Constants:
	private final static int BORDER_X = 10;
	private final static int BORDER_Y = 10;
	private final static int BORDER_Z = 5;
	private final static double D2SPI = 2*Math.sqrt(Math.PI);
	private final static int SUPER = 5;
	
	/** Default constructor. */
	public SWC2IMG() {
		
		// Retrieve last parameters:
		file = Prefs.get("swc2img.file",file);
		voxel = Prefs.get("swc2img.voxel",voxel);
		lateral = Prefs.get("swc2img.lateral",lateral);
		axial = Prefs.get("swc2img.axial",axial);
		background = Prefs.get("swc2img.background",background);
		correlation = Prefs.get("swc2img.correlation",correlation);
		snr = Prefs.get("swc2img.snr",snr);
		rescale = Prefs.get("swc2img.rescale",rescale);
		log = Prefs.get("swc2img.log",log);
		
		// Initialize formatters:
		dfm.decs(9); mfm.decs(1);
	}
	
	public void run(String arg) {
		
		// Construct dialog:
		final GenericDialog gd = new GenericDialog(name(),IJ.getInstance());
		
		// Add SWC parameters:
		gd.setInsets(0,0,0); gd.addMessage("SWC");
		gd.setInsets(0,2,0); gd.addCheckbox(" Voxel units for coordinates and radii in file",voxel);
		gd.addStringField("File path:",file,25);
		fileField = (TextField)gd.getStringFields().get(0);
		final Panel panel = new Panel(new FlowLayout(FlowLayout.LEFT,0,0));
		browseButton = new Button("    Browse    ");
		browseButton.addActionListener(this);
		panel.add(browseButton);
		gd.addPanel(panel,GridBagConstraints.WEST,new Insets(0,0,0,0));
		final GridBagLayout gbl = (GridBagLayout)gd.getLayout();
		final GridBagConstraints gbc = gbl.getConstraints(panel);
		gbc.gridx = 1; gbl.setConstraints(panel,gbc);
		
		// Add IMAGE parameters:
		gd.setInsets(10,0,0); gd.addMessage("IMAGE");
		gd.addNumericField("Lateral voxel size:",lateral,decimals(lateral),10,"["+dfm.d2s(MIN_LATERAL)+"..."+dfm.d2s(MAX_LATERAL)+"] micron");
		gd.addNumericField("Axial voxel size:",axial,decimals(axial),10,"["+dfm.d2s(MIN_AXIAL)+"..."+dfm.d2s(MAX_AXIAL)+"] micron");
		
		// Add OPTICS parameters:
		gd.setInsets(10,0,0); gd.addMessage("QUALITY");
		gd.addNumericField("Background:",background,decimals(background),10,"["+dfm.d2s(MIN_BACKGROUND)+"..."+dfm.d2s(MAX_BACKGROUND)+"] intensity");
		gd.addNumericField("Correlation:",correlation,decimals(correlation),10,"["+dfm.d2s(MIN_CORRELATION)+"..."+dfm.d2s(MAX_CORRELATION)+"] micron");
		gd.addNumericField("SNR:",snr,decimals(snr),10,"["+dfm.d2s(MIN_SNR)+"..."+dfm.d2s(MAX_SNR)+"] ratio");
		
		// Add other parameters:
		gd.setInsets(10,2,0); gd.addCheckbox(" Rescale to 8-bit image",rescale);
		gd.setInsets(0,2,0); gd.addCheckbox(" Log messaging",log);
		
		// Add dialog buttons:
		gd.addHelp("https://imagescience.org/meijering/software/swc2img/");
		gd.setHelpLabel("   Help   ");
		gd.setCancelLabel(" Cancel ");
		gd.setOKLabel("    Run    ");
		
		// Show dialog:
		if (position.x >= 0 && position.y >= 0) {
			gd.centerDialog(false);
			gd.setLocation(position);
		} else gd.centerDialog(true);
		gd.addWindowListener(this);
		gd.showDialog();
		
		if (gd.wasCanceled()) return;
		
		// Get parameters:
		file = gd.getNextString();
		voxel = gd.getNextBoolean();
		lateral = gd.getNextNumber();
		axial = gd.getNextNumber();
		background = gd.getNextNumber();
		correlation = gd.getNextNumber();
		snr = gd.getNextNumber();
		rescale = gd.getNextBoolean();
		log = gd.getNextBoolean();
		
		// Store parameters:
		Prefs.set("swc2img.file",file);
		Prefs.set("swc2img.voxel",voxel);
		Prefs.set("swc2img.lateral",lateral);
		Prefs.set("swc2img.axial",axial);
		Prefs.set("swc2img.background",background);
		Prefs.set("swc2img.correlation",correlation);
		Prefs.set("swc2img.snr",snr);
		Prefs.set("swc2img.rescale",rescale);
		Prefs.set("swc2img.log",log);
		
		// Initialize:
		if (IJ.getLog() != null) log(" ");
		log(name()+" "+version());
		
		final Timer timer = new Timer();
		timer.messenger.log(log);
		timer.start();
		
		boolean success = false;
		
		try {
			// Check parameters:
			if (Double.isNaN(lateral)) throw new IllegalArgumentException("The specified lateral voxel size is invalid");
			if (lateral < MIN_LATERAL) throw new IllegalArgumentException("The specified lateral voxel size is too small");
			if (lateral > MAX_LATERAL) throw new IllegalArgumentException("The specified lateral voxel size is too large");
			
			if (Double.isNaN(axial)) throw new IllegalArgumentException("The specified axial voxel size is invalid");
			if (axial < MIN_AXIAL) throw new IllegalArgumentException("The specified axial voxel size is too small");
			if (axial > MAX_AXIAL) throw new IllegalArgumentException("The specified axial voxel size is too large");
			
			if (Double.isNaN(background)) throw new IllegalArgumentException("The specified background is invalid");
			if (background < MIN_BACKGROUND) throw new IllegalArgumentException("The specified background is too small");
			if (background > MAX_BACKGROUND) throw new IllegalArgumentException("The specified background is too large");
			
			if (Double.isNaN(correlation)) throw new IllegalArgumentException("The specified correlation is invalid");
			if (correlation < MIN_CORRELATION) throw new IllegalArgumentException("The specified correlation is too small");
			if (correlation > MAX_CORRELATION) throw new IllegalArgumentException("The specified correlation is too large");
			
			if (Double.isNaN(snr)) throw new IllegalArgumentException("The specified SNR is invalid");
			if (snr < MIN_SNR) throw new IllegalArgumentException("The specified SNR is too small");
			if (snr > MAX_SNR) throw new IllegalArgumentException("The specified SNR is too large");
			
			// Load specified file:
			final Neuron neuron = load(file);
			
			// Image the neuron:
			final ImagePlus image = image(neuron);
			
			// Show image:
			if (image != null) {
				image.show();
				success = true;
			}
			
		} catch (final IllegalArgumentException e) {
			error(e.getMessage());
			
		} catch (final OutOfMemoryError e) {
			error("There is not enough memory to run the program");
			
		} catch (final Throwable e) {
			error("An unidentified error occurred while running the program");
			IJ.handleException(e);
		}
		
		final String finish = "Program finished ";
		if (success) log(finish+"successfully");
		else log(finish+"unsuccessfully");
		
		IJ.showProgress(1);
		IJ.showStatus("");
		
		timer.stop();
	}
	
	static String name() {
		
		return "SWC2IMG";
	}
	
	static String version() {
		
		final String version = SWC2IMG.class.getPackage().getImplementationVersion();
		
		return (version == null) ? "DEV" : version;
	}
	
	private void error(final String message) {
		
		log(message);
		IJ.showMessage(name()+": Error",message+".");
	}
	
	private void log(final String message) {
		
		if (log) IJ.log(message);
	}
	
	private Neuron load(final String file) {
		
		LineReader swc = null;
		
		Neuron neuron = null;
		
		try {
			// Load nodes and links from file:
			swc = new LineReader(new FileReader(file));
			log("Loading \""+file+"\"");
			
			final Hashtable<Integer,Node> nodes = new Hashtable<Integer,Node>();
			final Hashtable<Integer,Integer> links = new Hashtable<Integer,Integer>();
			final Hashtable<Integer,Boolean> used = new Hashtable<Integer,Boolean>();
			
			double xs = 1;
			double ys = 1;
			double zs = 1;
			double rs = 1;
			if (voxel) {
				xs = lateral;
				ys = lateral;
				zs = axial;
				rs = lateral;
			}
			
			double xmin = Double.MAX_VALUE, xmax = -Double.MAX_VALUE;
			double ymin = Double.MAX_VALUE, ymax = -Double.MAX_VALUE;
			double zmin = Double.MAX_VALUE, zmax = -Double.MAX_VALUE;
			
			int roots = 0;
			String line = swc.read();
			while (line != null) {
				line = line.trim();
				if (!line.startsWith("#")) {
					final StringTokenizer st = new StringTokenizer(line);
					final int index = Integer.parseInt(st.nextToken());
					st.nextToken(); // Skip node type parameter
					final double x = xs*Double.parseDouble(st.nextToken());
					final double y = ys*Double.parseDouble(st.nextToken());
					final double z = zs*Double.parseDouble(st.nextToken());
					final double r = rs*Double.parseDouble(st.nextToken());
					final int parent = Integer.parseInt(st.nextToken());
					if (x-r < xmin) xmin = x-r; if (x+r > xmax) xmax = x+r;
					if (y-r < ymin) ymin = y-r; if (y+r > ymax) ymax = y+r;
					if (z-r < zmin) zmin = z-r; if (z+r > zmax) zmax = z+r;
					nodes.put(index,new Node(x,y,z,r));
					links.put(index,parent);
					used.put(index,false);
					if (parent < 0) ++roots;
				}
				line = swc.read();
			}
			
			if (nodes.size() == 0) throw new IllegalStateException("The specified file does not contain any nodes");
			
			log("File contains "+nodes.size()+" nodes of which "+roots+" root"+(roots!=1?"s":""));
			log("Range in x dimension = ["+dfm.d2s(xmin)+","+dfm.d2s(xmax)+"] micron");
			log("Range in y dimension = ["+dfm.d2s(ymin)+","+dfm.d2s(ymax)+"] micron");
			log("Range in z dimension = ["+dfm.d2s(zmin)+","+dfm.d2s(zmax)+"] micron");
			
			// Convert nodes and links to neuron:
			final Neuron N = new Neuron();
			final Set<Integer> indices = nodes.keySet();
			for (int index: indices) {
				final int parent = links.get(index);
				if (parent < 0) continue; // Skip root nodes
				final Node s = nodes.get(parent);
				if (s == null) throw new IllegalStateException("Parent node "+parent+" of node "+index+" is missing in the file");
				final Node e = nodes.get(index);
				N.add(new Segment(s,e));
				used.put(parent,true);
				used.put(index,true);
			}
			for (int index: indices) {
				if (used.get(index) == false) {
					final Node n = nodes.get(index);
					N.add(new Segment(n,n));
				}
			}
			
			neuron = N;
			
		} catch (final FileNotFoundException e) {
			error("The specified file does not exist");
			
		} catch (final IOException e) {
			error("The specified file could not be loaded");
			
		} catch (final IllegalStateException e) {
			error(e.getMessage());
			
		} catch (final NoSuchElementException e) {
			error("Line "+swc.lines()+" of the file is missing information");
			
		} catch (final NumberFormatException e) {
			error("A number in line "+swc.lines()+" of the file could not be interpreted");
			
		} catch (final Throwable e) {
			error("Line "+swc.lines()+" of the file could not be interpreted");
		}
		
		try { swc.close(); } catch (final Throwable e) { }
		
		return neuron;
	}
	
	private ImagePlus image(final Neuron neuron) {
		
		if (neuron == null) return null;
		
		// Determine global bounds:
		final Bounds bounds = neuron.bounds();
		int xm = FMath.floor(bounds.lower.x/lateral);
		int ym = FMath.floor(bounds.lower.y/lateral);
		int zm = FMath.floor(bounds.lower.z/axial);
		int xM = FMath.ceil(bounds.upper.x/lateral);
		int yM = FMath.ceil(bounds.upper.y/lateral);
		int zM = FMath.ceil(bounds.upper.z/axial);
		
		final int w = xM - xm;
		final int h = yM - ym;
		final int d = zM - zm;
		final int xr = BORDER_X;
		final int yr = BORDER_Y;
		final int zr = BORDER_Z;
		final int xR = xr + w;
		final int yR = yr + h;
		final int zR = zr + d;
		final int width = w + 2*BORDER_X;
		final int height = h + 2*BORDER_Y;
		final int depth = d + 2*BORDER_Z;
		final int size = neuron.size();
		
		// Allocate arrays:
		final long volume = 5L*width*height*depth;
		log("Allocating "+memory(volume)+" for mapping and imaging");
		final float[][] fmp = new float[depth][width*height];
		final byte[][] bmp = new byte[depth][width*height];
		
		// Determine offset:
		final int xo = BORDER_X - xm;
		final int yo = BORDER_Y - ym;
		final int zo = BORDER_Z - zm;
		log("Offset in x dimension = "+dfm.d2s(xo*lateral)+" micron ("+dfm.d2s(xo)+" voxel)");
		log("Offset in y dimension = "+dfm.d2s(yo*lateral)+" micron ("+dfm.d2s(yo)+" voxel)");
		log("Offset in z dimension = "+dfm.d2s(zo*axial)+" micron ("+dfm.d2s(zo)+" voxel)");
		
		// Start progress indication:
		log("Making image of neuron structure");
		final Progressor pgs = new Progressor();
		pgs.display(true); pgs.updates(50);
		
		// Make voxel map:
		pgs.status("Voxel mapping...");
		pgs.range(0,0.1);
		pgs.steps(size);
		pgs.start();
		for (int i=0; i<size; ++i, pgs.step()) {
			final Segment s = neuron.get(i);
			final Bounds b = s.bounds();
			xm = xo + FMath.floor(b.lower.x/lateral);
			ym = yo + FMath.floor(b.lower.y/lateral);
			zm = zo + FMath.floor(b.lower.z/axial);
			xM = xo + FMath.ceil(b.upper.x/lateral);
			yM = yo + FMath.ceil(b.upper.y/lateral);
			zM = zo + FMath.ceil(b.upper.z/axial);
			for (int z=zm; z<zM; ++z) {
				final float[] fmpz = fmp[z];
				for (int y=ym; y<yM; ++y) {
					final int ywidth = y*width;
					for (int x=xm; x<xM; ++x)
						fmpz[ywidth+x] = 1;
				}
			}
		}
		
		// Make partial-volume image:
		pgs.status("Partial voluming...");
		pgs.range(0.1,0.6);
		pgs.steps((zR-zr)*(yR-yr));
		pgs.start();
		final double full = SUPER*SUPER*SUPER;
		double max = -Double.MAX_VALUE;
		for (int z=zr; z<zR; ++z) {
			final float[] fmpz = fmp[z];
			for (int y=yr; y<yR; ++y, pgs.step()) {
				final int ywidth = y*width;
				for (int x=xr; x<xR; ++x) {
					final int o = ywidth+x;
					if (fmpz[o] > 0) {
						int count = 0;
						for (int k=0; k<SUPER; ++k) {
							final double dz = (z - zo + (k + 0.5)/SUPER)*axial;
							for (int j=0; j<SUPER; ++j) {
								final double dy = (y - yo + (j + 0.5)/SUPER)*lateral;
								for (int i=0; i<SUPER; ++i) {
									final double dx = (x - xo + (i + 0.5)/SUPER)*lateral;
									if (neuron.contains(dx,dy,dz)) ++count;
								}
							}
						}
						if (count > max) max = count;
						fmpz[o] = count;
					}
				}
			}
		}
		
		if (max <= 0) {
			error("No neuron structure was found with the used sampling density");
			return null;
		}
		
		// Make Poisson image:
		pgs.status("Poisson imaging...");
		pgs.range(0.6,0.95);
		pgs.steps(depth*height);
		pgs.start();
		double snrc = snr;
		if (correlation > 0) {
			pgs.range(0.6,0.8);
			// Compensate SNR for noise reduction due to correlation:
			snrc /= D2SPI*correlation/lateral;
			snrc /= Math.sqrt(D2SPI*correlation/axial);
		}
		final PoissonGenerator poisson = new PoissonGenerator();
		// SNR = (S - B)/Sqrt(S) where S is peak and B is background:
		final double temp = snrc + Math.sqrt(snrc*snrc + 4*background);
		final double peak = 0.25*temp*temp;
		final double clip = peak + 3*Math.sqrt(peak);
		final double ffac = (peak - background)/max;
		for (int z=0; z<depth; ++z) {
			final float[] fmpz = fmp[z];
			for (int y=0, o=0; y<height; ++y, pgs.step()) {
				for (int x=0; x<width; ++x, ++o) {
					fmpz[o] = (float)poisson.next(ffac*fmpz[o] + background);
				}
			}
		}
		
		// Make correlated noise:
		if (correlation > 0) {
			pgs.range(0.8,0.95);
			pgs.status("Correlating noise...");
			final Image corimg = Image.wrap(image(fmp,width,height));
			final Smoother smoother = new Smoother();
			smoother.progressor.parent(pgs);
			smoother.gauss(corimg,correlation);
		}
		
		// Make output image:
		ImagePlus image = null;
		pgs.status("Creating image...");
		pgs.range(0.95,rescale?0.975:1);
		pgs.steps(depth*height);
		pgs.start();
		for (int z=0; z<depth; ++z) {
			final float[] fmpz = fmp[z];
			for (int y=0, o=0; y<height; ++y, pgs.step()) {
				for (int x=0; x<width; ++x, ++o) {
					fmpz[o] = (float)FMath.clip(fmpz[o],0,clip);
				}
			}
		}
		if (rescale) {
			pgs.status("Rescaling image...");
			pgs.range(0.975,1);
			pgs.steps(depth*height);
			pgs.start();
			final double bfac = (clip>0)?(255/clip):0;
			for (int z=0; z<depth; ++z) {
				final byte[] bmpz = bmp[z];
				final float[] fmpz = fmp[z];
				for (int y=0, o=0; y<height; ++y, pgs.step()) {
					for (int x=0; x<width; ++x, ++o) {
						bmpz[o] = (byte)FMath.round(bfac*fmpz[o]);
					}
				}
			}
			image = image(bmp,width,height);
			image.setDisplayRange(0,255);
		} else {
			image = image(fmp,width,height);
			image.setDisplayRange(0,clip);
		}
		
		pgs.stop();
		
		// Return image object:
		image.setTitle((new File(file)).getName());
		return image;
	}
	
	private Button browseButton;
	
	private TextField fileField;
	
	private ImagePlus image(final byte[][] array, final int width, final int height) {
		
		final ImageStack stack = new ImageStack(width,height);
		for (int z=0; z<array.length; ++z) stack.addSlice("",array[z]);
		final ImagePlus image = new ImagePlus("",stack);
		calibrate(image);
		return image;
	}
	
	private ImagePlus image(final float[][] array, final int width, final int height) {
		
		final ImageStack stack = new ImageStack(width,height);
		for (int z=0; z<array.length; ++z) stack.addSlice("",array[z]);
		final ImagePlus image = new ImagePlus("",stack);
		calibrate(image);
		return image;
	}
	
	private void calibrate(final ImagePlus image) {
		
		final Calibration c = new Calibration();
		c.pixelWidth = lateral;
		c.pixelHeight = lateral;
		c.pixelDepth = axial;
		c.setUnit("um");
		image.setCalibration(c);
	}
	
	public void actionPerformed(final ActionEvent e) {
		
		if (e.getSource() == browseButton) {
			final FileDialog fdg = new FileDialog(IJ.getInstance(),name()+": File",FileDialog.LOAD);
			fdg.setFile(""); fdg.setVisible(true);
			final String d = fdg.getDirectory();
			final String f = fdg.getFile();
			fdg.dispose();
			if (d != null && f != null) {
				String path = d + f;
				if (File.separator.equals("\\"))
					path = path.replace('\\','/');
				fileField.setText(path);
			}
		}
	}
	
	public void windowActivated(final WindowEvent e) { }
	
	private static Point position = new Point(-1,-1);
	
	public void windowClosed(final WindowEvent e) {
		
		position.x = e.getWindow().getX();
		position.y = e.getWindow().getY();
	}
	
	public void windowClosing(final WindowEvent e) { }
	
	public void windowDeactivated(final WindowEvent e) { }
	
	public void windowDeiconified(final WindowEvent e) { }
	
	public void windowIconified(final WindowEvent e) { }
	
	public void windowOpened(final WindowEvent e) { }
	
	private final Formatter dfm = new Formatter();
	
	private int decimals(final double d) {
		
		final String value = dfm.d2s(d);
		final int index = value.indexOf('.');
		if (index < 0) return 0;
		return value.length() - index - 1;
	}
	
	private final Formatter mfm = new Formatter();
	
	private String memory(final long volume) {
		
		if (volume < 1000) return volume+" B";
		if (volume < 1000000) return mfm.d2s(volume/1000d)+" kB";
		if (volume < 1000000000) return mfm.d2s(volume/1000000d)+" MB";
		return mfm.d2s(volume/1000000000d)+" GB";
	}
	
}

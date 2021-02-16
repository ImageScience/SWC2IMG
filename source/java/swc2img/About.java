package swc2img;

import ij.IJ;
import ij.ImageJ;
import ij.plugin.BrowserLauncher;
import ij.plugin.PlugIn;

/** ImageJ plugin for loading the <a href="https://imagescience.org/meijering/software/swc2img/" target="_blank">SWC2IMG website</a> in the default internet browser. */
public final class About implements PlugIn {
	
	/** Default constructor. */
	public About() {
	}
	
	public void run(final String arg) {
		
		try { BrowserLauncher.openURL("https://imagescience.org/meijering/software/swc2img/"); }
		
		catch (Throwable e) {
			final String message = "Could not open default internet browser";
			if (IJ.getInstance() != null) IJ.showMessage(SWC2IMG.name()+": Error",message+".");
			else System.out.println(SWC2IMG.name()+": "+message);
		}
	}
	
}

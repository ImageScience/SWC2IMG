package swc2img;

import imagescience.shape.Bounds;
import imagescience.utility.FMath;

import java.util.Vector;

final class Neuron {
	
	Neuron() {
		
		bounds.lower.x = Double.MAX_VALUE;
		bounds.lower.y = Double.MAX_VALUE;
		bounds.lower.z = Double.MAX_VALUE;
		bounds.upper.x = -Double.MAX_VALUE;
		bounds.upper.y = -Double.MAX_VALUE;
		bounds.upper.z = -Double.MAX_VALUE;
	}
	
	void add(final Segment s) {
		
		segments.add(s);
		final Bounds b = s.bounds();
		if (b.lower.x < bounds.lower.x) bounds.lower.x = b.lower.x;
		if (b.lower.y < bounds.lower.y) bounds.lower.y = b.lower.y;
		if (b.lower.z < bounds.lower.z) bounds.lower.z = b.lower.z;
		if (b.upper.x > bounds.upper.x) bounds.upper.x = b.upper.x;
		if (b.upper.y > bounds.upper.y) bounds.upper.y = b.upper.y;
		if (b.upper.z > bounds.upper.z) bounds.upper.z = b.upper.z;
		remap = true;
	}
	
	private final Bounds bounds = new Bounds();
	
	Bounds bounds() {
		
		return bounds;
	}
	
	private final Vector<Segment> segments = new Vector<Segment>();
	
	Segment get(final int index) {
		
		return segments.get(index);
	}
	
	int size() {
		
		return segments.size();
	}
	
	@SuppressWarnings("unchecked")
	boolean contains(final double x, final double y, final double z) {
		
		if (remap) remap();
		
		if (dx == 0 || dy == 0)
			return false;
		
		final int xb = FMath.floor((x - bounds.lower.x)*dx);
		final int yb = FMath.floor((y - bounds.lower.y)*dy);
		
		if (xb < 0 || xb >= BINS ||
			yb < 0 || yb >= BINS ||
			map[yb][xb] == null)
			return false;
		
		for (Segment s: (Vector<Segment>)map[yb][xb])
			if (s.contains(x,y,z))
				return true;
		
		return false;
	}
	
	private String path = null;
	
	void path(final String path) { this.path = path; }
	
	String path() { return path; }
	
	private Vector<String> file = null;
	
	void file(final Vector<String> file) { this.file = file; }
	
	Vector<String> file() { return file; }
	
	private final static int BINS = 1000;
	
	private final Object[][] map = new Object[BINS][BINS];
	
	private boolean remap = true;
	
	private double dx, dy;
	
	@SuppressWarnings("unchecked")
	private void remap() {
		
		remap = false;
		
		for (int y=0; y<BINS; ++y)
			for (int x=0; x<BINS; ++x)
				map[y][x] = null;
		
		dx = dy = 0;
		
		if (bounds.upper.x - bounds.lower.x <= 0 ||
			bounds.upper.y - bounds.lower.y <= 0)
			return;
		
		dx = BINS/(bounds.upper.x - bounds.lower.x);
		dy = BINS/(bounds.upper.y - bounds.lower.y);
		
		for (Segment s: segments) {
			final Bounds b = s.bounds();
			final int xb = FMath.floor((b.lower.x - bounds.lower.x)*dx);
			final int yb = FMath.floor((b.lower.y - bounds.lower.y)*dy);
			int xB = FMath.floor((b.upper.x - bounds.lower.x)*dx);
			int yB = FMath.floor((b.upper.y - bounds.lower.y)*dy);
			if (xB == BINS) xB = BINS - 1;
			if (yB == BINS) yB = BINS - 1;
			for (int y=yb; y<=yB; ++y) {
				for (int x=xb; x<=xB; ++x) {
					if (map[y][x] == null)
						map[y][x] = new Vector<Segment>();
					((Vector<Segment>)map[y][x]).add(s);
				}
			}
		}
	}
	
}

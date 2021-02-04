package swc2img;

import imagescience.shape.Bounds;

final class Segment {
	
	private final Node s, e;
	
	Segment(final Node s, final Node e) {
		
		this.s = s;
		this.e = e;
		
		dx = e.x - s.x;
		dy = e.y - s.y;
		dz = e.z - s.z;
		
		len2 = dx*dx + dy*dy + dz*dz;
		len = Math.sqrt(len2);
		
		bounds = new Bounds();
		bounds.lower.x = s.x - s.r; bounds.upper.x = s.x + s.r;
		bounds.lower.y = s.y - s.r; bounds.upper.y = s.y + s.r;
		bounds.lower.z = s.z - s.r; bounds.upper.z = s.z + s.r;
		if (e.x - e.r < bounds.lower.x) bounds.lower.x = e.x - e.r;
		if (e.y - e.r < bounds.lower.y) bounds.lower.y = e.y - e.r;
		if (e.z - e.r < bounds.lower.z) bounds.lower.z = e.z - e.r;
		if (e.x + e.r > bounds.upper.x) bounds.upper.x = e.x + e.r;
		if (e.y + e.r > bounds.upper.y) bounds.upper.y = e.y + e.r;
		if (e.z + e.r > bounds.upper.z) bounds.upper.z = e.z + e.r;
	}
	
	private final Bounds bounds;
	
	Bounds bounds() {
		
		return bounds;
	}
	
	private final double dx, dy, dz, len, len2;
	
	boolean contains(final double x, final double y, final double z) {
		
		final double vx = x - s.x;
		final double vy = y - s.y;
		final double vz = z - s.z;
		final double in = dx*vx + dy*vy + dz*vz;
		
		if (in > 0) {
			if (in < len2) {
				final double prx = in/len, frx = prx/len;
				final double rad = frx*e.r + (1-frx)*s.r;
				if (vx*vx + vy*vy + vz*vz - prx*prx <= rad*rad)
					return true;
			} else if ((x-e.x)*(x-e.x)+(y-e.y)*(y-e.y)+(z-e.z)*(z-e.z) <= e.r*e.r)
				return true;
		} else if (vx*vx + vy*vy + vz*vz <= s.r*s.r)
			return true;
		
		return false;
	}
	
}

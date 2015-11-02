package Tests;

import static org.junit.Assert.*;

import org.junit.Test;

import DNA.Coordinate;

public class testCoordinate {

	@Test
	public void testAdd() {
		Coordinate a= new Coordinate();
		Coordinate b = new Coordinate(1, -.004, -889);
		assertEquals(a, new Coordinate());
		Coordinate c= b.add(a);
		assertEquals(b,c);
		Coordinate d = new Coordinate (5,10,5);
		Coordinate e= d.add(b);
		assertEquals(e, new Coordinate(6, 9.996, -884));
	}

	@Test
	public void testMultiplyBy() {
		Coordinate a=new Coordinate();
		assertEquals(a.multiplyBy(0),a);
		assertEquals(a.multiplyBy(17.29),a);
		assertEquals(a.multiplyBy(-292),a);
		Coordinate b = new Coordinate(-10, 1.5, 0.0);
		assertEquals(b.multiplyBy(0),a);
		assertEquals(b.multiplyBy(-0.1),new Coordinate(1, -.15,0));
	}

	@Test
	public void testSubtract() {
		Coordinate a= new Coordinate();
		Coordinate b = new Coordinate(1, -.004, -889);
		assertEquals(a, new Coordinate());
		Coordinate c= b.subtract(a);
		assertEquals(b,c);
		Coordinate d = new Coordinate (5,10,5);
		Coordinate e= d.subtract(b);
		assertEquals(e, new Coordinate(4, 10.004, 894));
		assertEquals(e, e.add(d).subtract(d));
	}

	@Test
	public void testSquare() {
		Coordinate a= new Coordinate();
		assertEquals(a.square(),0,.0001);
		Coordinate b = new Coordinate(0,0,-889);
		assertEquals(b.square(),889*889,.0001);
		Coordinate c = new Coordinate(5,-5,10);
		assertEquals(c.square(),c.magnitude()*c.magnitude(),.0001);
		assertEquals(c.square(),150,.0001);
	}

	@Test
	public void testMagnitude() {
		Coordinate a= new Coordinate();
		assertEquals(a.magnitude(),0,.0001);
		Coordinate b = new Coordinate(0,0,-889);
		assertEquals(b.magnitude(),889,.0001);
		Coordinate c = new Coordinate(5,-5,10);
		assertEquals(c.square(),c.magnitude()*c.magnitude(),.0001);
		assertEquals(c.magnitude(),Math.sqrt(150),.0001);
	}
	
	@Test
	public void testDistance() {
		Coordinate a= new Coordinate();
		Coordinate b = new Coordinate(1, -.004, -889);
		assertEquals(a.distance(b), a.subtract(b).magnitude(),.001);
		Coordinate c = new Coordinate(0,0,-889.7329);
		assertEquals(c.magnitude(),c.distance(a),.001);
		assertEquals(889.7326,c.distance(a),.001);
	}


}

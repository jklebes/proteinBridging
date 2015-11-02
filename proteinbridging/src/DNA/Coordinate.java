package DNA;

public class Coordinate {
	
double x;
double y;
double z;

public Coordinate(){
	this.x=0;
	this.y=0;
			this.z=0;
}
public Coordinate (double x, double y, double z){
	this.x=x;
	this.y=y;
	this.z=z;
}
public Coordinate add(Coordinate coordinate) {
	Coordinate sum = new Coordinate(this.x+coordinate.x, this.y+coordinate.y,this.z+coordinate.z);
	return sum;
}
public Coordinate multiplyBy(double factor) {
	Coordinate multiple = new Coordinate(this.x*factor, this.y*factor,this.z*factor);
	return multiple;
	
}
public Coordinate subtract(Coordinate coordinate) {
	Coordinate diff = new Coordinate(this.x-coordinate.x, this.y-coordinate.y,this.z-coordinate.z);
	return diff;
}
public double square() {
	double square=this.magnitude()*this.magnitude();
	return square;
}
public double magnitude() {
	double magnitude = Math.sqrt(x*x + z*z + y *y);
	return magnitude;
}

public double distance(Coordinate coordinate) {
	double distance= Math.sqrt((x-coordinate.x)*(x-coordinate.x) 
			+ (z-coordinate.z)*(z-coordinate.z) 
			+ (y-coordinate.y)*(y-coordinate.y));
	return distance;
}

public String toString(){
	String s = this.x + " "+ this.y+ " "+ this.z;
	return s;
}

@Override
public boolean equals(Object object){
	if (object == null) {
        return false;
    }
	else if (getClass() != object.getClass()) {
        return false;
    
} 
	final Coordinate other = (Coordinate) object;
	double t=.0001;
    if ((this.x<=other.x+t||this.x>=other.x-t) 
    		&& (this.y<=other.y+t||this.y>=other.y-t) 
    		&& (this.z<=other.z+t||this.z>=other.z-t)) {
        return true;
    }
    else {
    return false;
}

}

}

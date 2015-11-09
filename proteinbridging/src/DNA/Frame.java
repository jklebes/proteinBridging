package DNA;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class Frame {

	String filename;
	int n;
	int np;
	//all points including proteins atoms
	Coordinate[] rawpointsall;
	Coordinate[] offsetsall;
	Coordinate[] imagepointsall;
	//just polymer bead atoms
	Coordinate[] rawpoints;
	Coordinate[] offsets;
	Coordinate[] imagepoints;
	int boxsize = 100;
	FileReader filereader;
	BufferedReader bufferedreader;
	ArrayList<String> lines= new ArrayList<String>();
	private Coordinate[] rawpointsprotein;
	private Coordinate[] offsetsprotein;
	private Coordinate[] imagepointsprotein;
	private Coordinate[] proteincenters;	
	private Coordinate[] proteinsides;
	private FileReader thermofilereader;
	private BufferedReader thermobufferedreader;
	private ArrayList<String> thermolines =new ArrayList<String>();
	double[] epair;
	private int patches =1;
	private double[] timeattached;
	private double[] totaltimeattached;
	private int[] detachmentscounter;
	
	public Frame(String filename, int proteinnumber) throws IOException{
		this.filename=filename;
		readFile();
		this.n = Integer.parseInt(lines.get(3).trim());
		setns(proteinnumber);
		this.np=proteinnumber;
	}
	private void readFile() throws IOException{
		this.filereader = new FileReader(filename);
		this.bufferedreader= new BufferedReader(filereader);
		String line;
		while((line = bufferedreader.readLine()) != null) {
               lines.add(line);
            }
	}

	public void readThermoFile(String thermofilename) throws IOException{
		this.epair = new double[1200];
		this.thermofilereader = new FileReader(thermofilename);
		this.thermobufferedreader= new BufferedReader(thermofilereader);
		String line;
		while((line = thermobufferedreader.readLine()) != null) {
               thermolines.add(line);
            }
		thermolines.remove (0);
		for (int i=0;i<n;i++){
			String[] lineelements=thermolines.get(i).split(" ");
			epair[i] = Double.parseDouble(lineelements[2]);
		}
	}

	void setns(int np) {
		this.n = Integer.parseInt(lines.get(3).trim());
		this.np = np;
		this.imagepointsall=new Coordinate[n];
		this.offsetsall=new Coordinate[n];
		this.rawpointsall=new Coordinate[n];

		this.rawpoints= new Coordinate[n-np];
		this.offsets= new Coordinate[n-np];
		this.imagepoints= new Coordinate[n-np];
		
		this.rawpointsprotein= new Coordinate[np];
		
		this.proteincenters=new Coordinate[np/(patches+1)];
		this.proteinsides=new Coordinate[patches*np/(patches+1)];
		
		this.timeattached = new double[proteinsides.length];
		this.totaltimeattached = new double[proteinsides.length];
		this.detachmentscounter = new int[proteinsides.length];
	}

	void readNextPoints(int frame){
		//trim first 9 lines, leaving coordinates
		int trimnumber; //= (frame-1) * (n + 9) +8;
		
		//trim a bit each round
		if (frame==1){ trimnumber=8;}
		else {trimnumber=n+8;}
		for (int i=0;i<=trimnumber;i++){
		lines.remove (0);
		}
		for (int i=0;i<n;i++){
			String[] lineelements=lines.get(i).split(" ");
			int index= Integer.parseInt(lineelements[0])-1;
			rawpointsall[index]=new Coordinate(Double.parseDouble(lineelements[2]), Double.parseDouble(lineelements[3]), Double.parseDouble(lineelements[4]));
			offsetsall[index]=new Coordinate(Double.parseDouble(lineelements[5]), Double.parseDouble(lineelements[6]), Double.parseDouble(lineelements[7]));
		imagepointsall[index]=rawpointsall[index].add(offsetsall[index]);
		}
		convertUnits();
	}
	
	private void convertUnits() {
		for (int i =0; i<n; i++){
			rawpointsall[i] = rawpointsall[i].multiplyBy(boxsize);
			offsetsall[i] = offsetsall[i].multiplyBy(boxsize);
			imagepointsall[i] = imagepointsall[i].multiplyBy(boxsize);
		}
	}
	
	void trimProteins(){
	for (int i=0;i<(n-np);i++){
		rawpoints[i]=rawpointsall[i];
		offsets[i]=offsetsall[i];
	imagepoints[i]=imagepointsall[i];
	}
}
	
	public double calculateRGyr(){
	Coordinate mean = calculateMean();
	double sum=0;
	for (int i=0;i<(n-np);i++){
		sum+= imagepoints[i].subtract(mean).square();
	}
	double rgyr=Math.sqrt(sum / (double)(n-np));
	//System.out.println(frame);
	//System.out.println(rgyr);
	return rgyr;
	}
	
	private Coordinate calculateMean() {
		Coordinate mean=new Coordinate();
		for (int i=0;i<n-np;i++){
		mean=mean.add(imagepoints[i]);
		}
		mean=mean.multiplyBy(1/(double)(n-np));
		return mean;
	}
	
	void makeProteinArrays(){
	for (int i=0;i<np;i++){
		rawpointsprotein[i]=rawpointsall[(i+(n-np))];
		
		//i = 299 last DNA bead 
		//i =300 protein center, i=301, 302 protein outsides
		//(300,301,302,303 in dump file numbering)
		if (i%(patches+1)==0){
			proteincenters[i/(patches+1)]=rawpointsprotein[i];	
		}
		else{
			int index=(int) (patches*i/(patches+1));
			proteinsides[index]= rawpointsprotein[i];
		}
		}
	
	}
	
	//checks 1 point against list of DNA beads
	boolean checkAttached(Coordinate protein){
		double radius =1.0/2.0; //diameter is 1
		boolean attached=false;
		for (int i=0;i<(n-np);i++){
			attached = attached || (rawpoints[i].distance(protein)<radius);
		}
		return attached;
	}
	
	//counts attached protein sides
	int countAttached(){
		int counter=0;
		for (int i=0;i<proteinsides.length;i++){
			if (checkAttached(proteinsides[i])){
				counter += 1;
				timeattached[i]+= 1000*.005;
		}
			else{ //if detached add time to avg calculation, reset counter
				totaltimeattached[i]+= timeattached[i];
				timeattached[i]=0;
				detachmentscounter[i]+=1;
				}
		}
			return counter;
	}
	
	public double avgTimeAttached(){
		double total=0;
		int number = 0;
		for (int i=0;i<proteinsides.length;i++){
			total += totaltimeattached[i];
			number += detachmentscounter[i];
			
		}
		double avg = total/(double)number;
			return avg;
	}
	
	
	public double totalTimeAttached(){
		double total=0;
		for (int i=0;i<proteinsides.length;i++){
			total += totaltimeattached[i];
			
		}
		return total;
	}
	public double errorAvgTimeAttached() {
		double total=0;
		int number = 0;
		double timesquaredtotal=0;
		for (int i=0;i<proteinsides.length;i++){
			total += totaltimeattached[i];
			number += detachmentscounter[i];
			timesquaredtotal += totaltimeattached[i]*totaltimeattached[i];
		}
		double avgtimesquared = (total/(double)number)*(total/(double)number);
		double timesquaredavg = timesquaredtotal / (double)number;
		double diff = timesquaredavg - avgtimesquared;
		double error = Math.sqrt(diff)/Math.sqrt(number);
		return error;
	}
	
}

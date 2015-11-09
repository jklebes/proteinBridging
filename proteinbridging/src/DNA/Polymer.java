package DNA;

import java.io.IOException;

public class Polymer {
	Frame frame;
	String filename;
	double[] rgyrs;
	int[] attached;
	int framenumber;
	int proteinnumber;

	public Polymer (String filename, int proteinnumber, int framenumber){
		this.filename=filename;
		this.framenumber=framenumber;
		rgyrs=new double[framenumber];
		attached=new int[framenumber];
		this.proteinnumber=proteinnumber;
	}

	
	public void CalculateAllRgyr() throws IOException{
		frame = new Frame(filename, proteinnumber);
		for (int i=1 ; i<=framenumber;i++){
			System.out.println(i);
			frame.readNextPoints(i);
			frame.trimProteins();
			rgyrs[i-1]=frame.calculateRGyr();
		}
	}

	public double avgRgyrs(int from, int to){
		double sum=0;
		for (int i=(from-1); i<=(to-1);i++){
			sum += rgyrs[i];
		}
		System.out.println(sum);
		System.out.println((to-from));
		double avg = sum / (double)(to-from);
		return avg;
	}
	
	public double errorRgyrs(int from, int to){
		double squaresum=0;
		double mean=avgRgyrs(from, to);
				for (int i=(from-1); i<=(to-1);i++){
					squaresum += (mean-rgyrs[i])*(mean-rgyrs[i]);
				}
		double number=(double)(to-from);
		double stdev=Math.sqrt(squaresum)/number;
		return stdev;
	}

	public double errorNumberAttached(int from, int to){
		double squaresum=0;
		double mean=avgAttached(from, to);
				for (int i=(from-1); i<=(to-1);i++){
					squaresum += (double)(mean-attached[i])*(double)(mean-attached[i]);
				}
		double number=(double)(to-from);
		double stdev=Math.sqrt(squaresum)/number;
		return stdev;
	}

	public void CalculateAllAttached() throws IOException {
		frame = new Frame(filename, proteinnumber);
		for (int i=1 ; i<=framenumber;i++){
			System.out.println(i);
			frame.readNextPoints(i);
			frame.trimProteins();
			frame.makeProteinArrays();
			attached[i-1]=frame.countAttached();
		}
	}


	public double avgAttached(int from, int to) {
		double sum=0;
		for (int i=(from-1); i<=(to-1);i++){
			sum += (double)(attached[i]);
		}
		System.out.println(sum);
		System.out.println((to-from));
		double avg = sum / (double)(to-from);
		return avg;
	}
	
	//go through each time frame
	public double avgTimeAttached() throws IOException {
		frame = new Frame(filename, proteinnumber);
		for (int i=1 ; i<=framenumber;i++){
			System.out.println(i);
			frame.readNextPoints(i);
			frame.trimProteins();
			frame.makeProteinArrays();
			
			//continuously updates time arrays in frame
			//also saves array of number attached each frame here but this is not needed
			attached[i-1]=frame.countAttached();
		}
		//return avg for this run from total time arrays created in frame]
		return frame.avgTimeAttached();
	}

	
	public double errorAvgTimeAttached() throws IOException {
		frame = new Frame(filename, proteinnumber);
		for (int i=1 ; i<=framenumber;i++){
			System.out.println(i);
			frame.readNextPoints(i);
			frame.trimProteins();
			frame.makeProteinArrays();
			
			//continuously updates time arrays in frame
			//also saves array of number attached each frame here but this is not needed
			attached[i-1]=frame.countAttached();
		}
		//return avg for this run from total time arrays created in frame]
		return frame.errorAvgTimeAttached();
	}
	
	//go through each time frame
	public double totalTimeAttached() throws IOException {
		frame = new Frame(filename, proteinnumber);
		for (int i=1 ; i<=framenumber;i++){
			System.out.println(i);
			frame.readNextPoints(i);
			frame.trimProteins();
			frame.makeProteinArrays();
			
			attached[i-1]=frame.countAttached();
		}
		//total time attached this frame, per protein
		return frame.totalTimeAttached() / proteinnumber;
	}



	public double getRgyr(int i) {

		return rgyrs[i];
	}
	
	public double avgEpair(int from, int to, String thermofile) throws IOException{
		Frame frame = new Frame(filename, proteinnumber);
		frame.setns(proteinnumber);
		frame.readThermoFile(thermofile);
		double sum=0;
		for (int i=(from-1); i<=(to-1);i++){
			sum += frame.epair[i];
		}
		System.out.println(sum);
		System.out.println((to-from));
		double avg = sum / (double)(to-from);
		return avg;
	}
	
}

package DNA;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Correlation {

	String xname;
	double[] xvalues;
	String yname;
	double[] yvalues;
	String line;
	FileWriter filewriter ;
	BufferedWriter bufferedwriter;
	double loopsize;
	private int framenumber;
	private int proteinnumber;
	
	public Correlation(int loopsize, int framenumber){
		this.loopsize = loopsize;
		this.framenumber=framenumber;
		this.xvalues=new double[loopsize];	
		this.yvalues=new double[loopsize];
	}
	
	public void startFile(String outfilename) throws IOException{
		filewriter = new FileWriter(outfilename);
		bufferedwriter = new BufferedWriter(filewriter);
	}
	
	public void energyRgyr(int from, double eIncrement) throws IOException{
		double[] errorvalues= new double[framenumber-from];
		for (int i=1 ; i<=loopsize;i++){
			xvalues[i-1]=i*eIncrement;
		String infilename ="/home/jsk/proteinBridging/numberLoopDNA/dump_run_"+i+"_200.DNAdumbells";
		Polymer p = new Polymer(infilename,proteinnumber,framenumber);
		p.CalculateAllRgyr();
		double rgyr=p.avgRgyrs(from, framenumber);
		System.out.println(rgyr);
		yvalues[i-1]=rgyr;
		errorvalues[i-1]=p.errorRgyrs(from, framenumber);
	}
		for (int i=0 ; i<loopsize;i++){
			line = xvalues[i]+ " "+ yvalues[i]+ " "+ errorvalues[i];
			bufferedwriter.write(line);
			bufferedwriter.newLine();
		}
			bufferedwriter.close();
			}
	
		
		public void numberRgyr(int from,double eIncrement) throws IOException{
			double[] errorvalues= new double[framenumber-from];
			for (int i=0 ; i<loopsize;i++){
				xvalues[i]=200- i*eIncrement;
			String infilename ="/home/jsk/proteinBridging/numerLoopChromatin/dump_run_"+i+"_200.DNAdumbells";
			Polymer p = new Polymer(infilename,(int) (600-i*eIncrement*3), framenumber);
			p.CalculateAllRgyr();
			double rgyr=p.avgRgyrs(from, framenumber);
			System.out.println(rgyr);
			yvalues[i]=rgyr;
			errorvalues[i]=p.errorRgyrs(from, framenumber);
			}
			for (int i=0 ; i<loopsize;i++){
				line = xvalues[i]+ " "+ yvalues[i]+ " "+ errorvalues[i];
				bufferedwriter.write(line);
				bufferedwriter.newLine();
			}
				bufferedwriter.close();
				
		}
			public void numberAttached(int from, double eIncrement) throws IOException{

				double[] errorvalues= new double[framenumber-from];
				for (int i=0 ; i<loopsize;i++){
					xvalues[i]=200- i*eIncrement;
				String infilename ="/home/jsk/proteinBridging/numberLoopDNA/dump_run_"+i+"_200.DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (600-i*eIncrement*3), framenumber);
				p.CalculateAllAttached();
				double attached=p.avgAttached(from, framenumber);
				System.out.println(attached);
				yvalues[i]=attached;
				errorvalues[i]=p.errorNumberAttached(from, framenumber);
				}
				for (int i=0 ; i<loopsize;i++){
					line = xvalues[i]+ " "+ yvalues[i]+ " "+ errorvalues[i];
					bufferedwriter.write(line);
					bufferedwriter.newLine();
				}
					bufferedwriter.close();
			}
		
			public void numberThermo(int from,double eIncrement) throws IOException{
				for (int i=0 ; i<loopsize;i++){
					xvalues[i]=200- i*eIncrement;
				String infilename ="/home/jsk/proteinBridging/numberLoopDNA/dump_run_"+i+"_200.DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (600-i*eIncrement*3), framenumber);
				double epair=p.avgEpair(from, framenumber, "/home/jsk/proteinBridging/numberLoopDNA/thermo_run_"+i+".DNAdumbells");
				System.out.println(epair);
				yvalues[i]=epair;
				}
				for (int i=0 ; i<loopsize;i++){
					line = xvalues[i]+ " "+ yvalues[i];
					bufferedwriter.write(line);
					bufferedwriter.newLine();
				}
					bufferedwriter.close();
					
			}
			
			public void forceTimeAttached(int from,double DIncrement) throws IOException{
				for (int i=0 ; i<loopsize;i++){
					xvalues[i]=DIncrement*(i+1);
				String infilename ="/home/jsk/proteinBridging/timeAttached/dump_run__D"+(i+1)+".DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (10), framenumber);
				double avg=p.avgTimeAttached();
				System.out.println(avg);
				yvalues[i]=avg;
				}
			}
			
			public void timeAttachedErrors(int from,int DIncrement) throws IOException{
				for (int i=0 ; i<loopsize;i++){
					xvalues[i]=DIncrement*(i+1);
				String infilename ="/home/s1203908/proteinBridging/timeAttached/dump_run__D"+(i+1)+".DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (10), framenumber);
				double avg=p.errorAvgTimeAttached();
				System.out.println(avg);
				yvalues[i]=avg;
				}
			}
			
			public void forceFractionAttached(int from,double DIncrement) throws IOException{
				double runtime = 6000000 * .005 ;
				for (int i=(int) (loopsize/2) ; i<loopsize;i++){
					xvalues[i]=DIncrement*(i+1);
				String infilename ="/home/s1203908/proteinBridging/timeAttached/dump_run__D"+(i+1)+".DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (10), framenumber);
				double total=p.totalTimeAttached();
				System.out.println(total);
				yvalues[i]=total / runtime;
				}
			}
			
			public void singleRgyr(double timeIncrement, int filenumber) throws IOException{
				this.xvalues=new double[framenumber];	
				this.yvalues=new double[framenumber];
				String infilename ="/home/jsk/proteinBridging/numberLoopChromatin/dump_run_"+filenumber+".DNAdumbells";
				Polymer p = new Polymer(infilename,(int) (300-filenumber*10*3),framenumber);
				p.CalculateAllRgyr();
				for (int i=0 ; i<framenumber;i++){
					xvalues[i]=i*timeIncrement;
				double rgyr=p.getRgyr(i);
				System.out.println(rgyr);
				yvalues[i]=rgyr;
				}
				for (int i=0 ; i<framenumber;i++){
					line = xvalues[i]+ " "+ yvalues[i];
					bufferedwriter.write(line);
					bufferedwriter.newLine();
				}
					bufferedwriter.close();}
			
		
		public void writeFile() throws IOException{
			for (int i=0 ; i<loopsize;i++){
			line = xvalues[i]+ " "+ yvalues[i];
			bufferedwriter.write(line);
			bufferedwriter.newLine();
		}
			bufferedwriter.close();}

		public void setFrameNumber(int f) {
			this.framenumber=f;
		}

		public void setProteinNumber(int np) {
			this.proteinnumber=np;
		}
}

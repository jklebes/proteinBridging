package DNA;

import java.io.IOException;

public class RunCorrelation {

	/**
	 * @param args
	 * @throws IOException 
	 */
	static int framenumber=6000;
	private static int loopnumber=20;
	private static int from = 1000;
	private static int numberincrement=1;
	
	public static void main(String[] args) throws IOException {
		//Correlation c= new Correlation(400, framenumber) ;
		//c.setProteinNumber(proteinnumber);
		//c.startFile("/home/jsk/proteinBridging/numberLoopChromatin/energyrgyr2.txt");
		//c.energyRgyr(10,1);
		//c.writeFile();
		//Correlation nr =new Correlation(loopnumber, framenumber) ;
		//nr.startFile("/home/jsk/proteinBridging/numerLoopChromatin/numberrgyr.txt");
		//nr.numberRgyr(from,numberincrement);
		//nr.writeFile();
		//Correlation na =new Correlation(loopnumber, framenumber) ;
		//na.startFile("/home/jsk/proteinBridging/numberLoopDNA/numberattached.txt");
		//na.numberAttached(from, numberincrement);
		//na.writeFile();
		//Correlation thermo =new Correlation(loopnumber, framenumber) ;
		//thermo.startFile("/home/jsk/proteinBridging/numberLoopDNA/numberthermo.txt");
		//thermo.numberThermo(from, numberincrement);
		//Correlation rt =new Correlation(11, framenumber) ;
		//rt.startFile("/home/jsk/proteinBridging/numberLoopChromatin/rgyrtime22.txt");
		//rt.singleRgyr(1,2);
		//
		//rt.writeFile();
		Correlation ne= new Correlation(loopnumber, framenumber) ;
ne.startFile("/home/s1203908/proteinBridging/D15/rgyrnp.txt");
		ne.numberRgyr(from,numberincrement);
		//dt.writeFile();
	}

}

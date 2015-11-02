package DNA;

import java.io.IOException;

public class FrameTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		for (int i=1 ; i<=51;i++){
		Frame f = new Frame("/home/jsk/proteinBridging/dump_run_7.DNA_only", 100);
	f.readNextPoints(i);
	f.trimProteins();
	f.calculateRGyr();
	//write in file
		}
	}

	
}

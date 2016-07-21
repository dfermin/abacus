/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package mainFunction;

import abacus.abacus;
import abacus.abacusUI.abacusUI;
import abacus.globals;
import java.io.File;
import java.io.IOException;
import org.hsqldb.util.DatabaseManagerSwing;

/**
 *
 * The sole purpose of this class is to determine what type of interface
 * the user wants. With the proper command line options, the program starts
 * without using the gui. Under all other circumstances, the gui is started.
 *
 * @author dfermin
 */
public class mainFunction {

	public static void main(String[] args) throws IOException {
		abacusUI gui = null;
		abacus cmdLine = null;
		String inputFile = null;
		DatabaseManagerSwing dbGUI = null;

		/*
		 * Insufficient command line arguments, start the gui
		 */
		if( (args.length == 1) && (args[0].equals("-dbgui")) ) {
			// execute HyperSQL gui
			dbGUI = new DatabaseManagerSwing();
			dbGUI.main();
		}
                else if( (args.length == 1) && (args[0].equals("-t")) ) {
                    globals.writeTemplate();
                    System.exit(0);
                }
		else if( (args.length == 2) && (args[0].equals("-p")) ) {
			// run command line version using input file
			inputFile = args[1];
			File tf = new File(inputFile);
			if( tf.exists() ) {
				cmdLine = new abacus();
				cmdLine.main(args);
			}
		}
		else { 
			// execute Abacus GUI
			
			// Check to see if you can run the gui. Need to be running a 
			// desktop environment for that
			if( java.awt.GraphicsEnvironment.isHeadless() ) {
				cmdLine = new abacus();
				System.err.print("\n\n" + cmdLine.printHeader() + "\n");
				System.err.print("\nERROR!\n" +
					"I was unable to start the GUI. Perhaps you are using a " +
					"remote terminal connection?\n\n" +
					"Recommended command line usage: java -jar abacus.jar -p <parameter_file.txt>\n" +
					"For database interface type: java -jar abacus.jar -dbgui\n\n"
					//"For database interface type: java -cp \".:./abacus.jar\" org.hsqldb.util.DatabaseManagerSwing\n\n"
				);
				System.exit(-1);
			}
			else {
				gui = new abacusUI();
				gui.main(args);
			}
		}
//		else if(args.length != 2) {
//			if( java.awt.GraphicsEnvironment.isHeadless() ) {
//				cmdLine = new abacus();
//				System.err.print("\n\n" + cmdLine.printHeader() + "\n");
//				System.err.print("\nERROR!\n" +
//					"I was unable to start the GUI. Perhaps you are using a " +
//					"remote terminal connection?\n\n" +
//					"Recommended command line usage: java -jar abacus.jar -p <parameter_file.txt>\n" +
//					"For database interface type: java -cp \".:./abacus.jar\" org.hsqldb.util.DatabaseManagerSwing\n\n"
//				);
//				System.exit(-1);
//			}
//			gui = new abacusUI();
//			gui.main(args);
//		}
//		 else {
//			// just run command line version
//			if(args[0].equals("-p")) {
//				inputFile = args[1];
//				File tf = new File(inputFile);
//				if( tf.exists() ) {
//					cmdLine = new abacus();
//					cmdLine.main(args);
//				}
//			}
//			if(args[0].equals("-db")) {
//				// start GUI interface to HyperSQL DB
//				dbGUI = new DatabaseManagerSwing();
//				dbGUI.main();
//			}
//			else {
//				// provided command line options were incorrect, start gui
//				gui = new abacusUI();
//				gui.main(args);
//			}
//		 }
	}

}

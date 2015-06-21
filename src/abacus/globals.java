package abacus;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.xml.stream.XMLStreamReader;
import javax.swing.GroupLayout.*;
import abacus_textArea.abacus_textArea;
import java.text.DecimalFormat;


/*
 * Global functions and variables go here
 */
public class globals {

	public static String OStype = null;
	public static String srcDir = null;
	public static String paramFile = null;
	public static String DBname = null;
	public static String gene2protFile = null;

	public static String outputFilePath = null;  // full path to output file
	public static String outputPath = null;  // holds the parent path of the output file

	public static String combinedFilePath = null;
	public static String combinedFile = null;

	public static String decoyTag = null;
	public static String pepXMLsuffix = "pep.xml";
	public static String protXMLsuffix = "prot.xml";

	// these variables are used to record peptide modifications the user wants
	// consider or discard in the data set.
	public static String pepRegexText = null;
	public static String[ ] pepMods_plus = null;
	public static String[ ] pepMods_minus = null;


	public static boolean keepDB = false;  // true means the user wants to keep the database file
	public static boolean recalcPepWts = false; // true means that the peptide weights in each expt. should be adjusted
	public static boolean byGene = false;  // true means the user wants gene-centric output
	public static boolean byPeptide = false; // true means the user wants peptide-level results
	public static int outputFormat = -1; // determines what kind of output the user wants.
	public static boolean genesHaveDescriptions = false; // true means the parsed gene2prot file has gene descriptions in it
	public static boolean doNSAF = false; // true means the user wants the NSAF formatted spectral counts

	public static double NSAF_FACTOR = -1; // used to deal with rounding errors in computing NSAF
	public static double maxIniProbTH = -100;
	public static double iniProbTH = -100;
	public static double minCombinedFilePw = -100;
	public static double minPw = -100;
	public static double epiThreshold = -100; // epi = experimental peptide probability inclusion threshold
	
	public static String fastaFile = null;
	public static HashMap<String, Integer> protLen;

	public static boolean makeVerboseOutput = false;

	// maps pepXML files to their protXML files
	// key = pepXML, value = parent protXML file
	public static HashMap<String, String> pepTagHash = new HashMap<String, String>();
	
	// maps protXML files to their tag
	// key = protXML, value = short tag
	public static HashMap<String, String> protTagHash = new HashMap<String, String>();

	public static List<String> pepXmlFiles = new ArrayList<String>();
	public static List<String> protXmlFiles = new ArrayList<String>();
	
	public static String fileSepChar = System.getProperty("file.separator"); // get either '\' or '/'
	
	// constants used for determining what type of output the user wants
	public static final int defaultOutput = 0;
	public static final int protQspecFormat = 1;
	public static final int geneQspecFormat = 4;
	public static final int customOutput = 2;
	public static final int geneOutput = 3;
	public static final int peptideOutput = 5;
	
	
	public static boolean proceedWithQuery = false; //used to deterin if a prepared statement should proceed or not 
	
	// variables used for custom output.
	public static Set<String> printC = new HashSet<String>();
	public static Set<String> printE = new HashSet<String>();
 
	//ERROR CODES
	public static int noDBname = 1;
	public static int noCombinedTag = 2;
	public static int noSrcDir = 3;
	public static int noDecoyTag = 4;
	public static int noMaxIniProbTH = 5;
	public static int noIniProbTH = 6;
	public static int noMinCombinedFilePw = 7;
	public static int noMinPw = 8;
	public static int DirError = 9;
	public static int noFastaFile = 10;
	public static int fastaNotFound = 11;
	public static int paramFileNotFound = 12;
	public static int mapFileNotFound = 13;
	public static int paramFileNull = 14;
	public static int outputPathNotFound = 15;

	
	
	/***************************************************************************
	 * 
	 *                        Functions
	 *
	 **************************************************************************/


	/*
	 *  Function reads in FASTA file and records each protein's length
	 */
	public static boolean parseFasta(abacus_textArea console) {
		File fastaF = null;
		InputStream in = null;
		DataInputStream dis = null;
		BufferedReader br = null;
		
		protLen = new HashMap<String, Integer>();
		String k = null;
		String seq = null;
		String line = null;
		String err = null; //holds stderr text

		int len = 0;
		boolean firstLine = true; // false means the line is NOT the first one
		boolean everythingOK = true;
		boolean printOnce = true; // false means the line has been printed once already

		try {
			fastaF = new File( fastaFile );
			if( !fastaF.exists() ) return true;
			
			in = new FileInputStream(fastaF);
			dis = new DataInputStream(dis);
			br = new BufferedReader(new InputStreamReader(in));
			k = "";
			seq = "";
			
			while( (line = br.readLine()) != null ) {

				if(firstLine) {
					if( !line.startsWith(">") ) { // poorly formatted fasta file
						err = "ERROR! '" + fastaF.getName() + "' is not a properly formatted FASTA File.\n";
						
						if (console == null) { System.err.print(err); System.exit(-1); }
						else console.append(err);
						
						everythingOK = false;
						return true;
						//break; // get out of while loop
					}
					firstLine = false;
				}

				
				if(line.startsWith(">")) {
					
					if(k.isEmpty() || k == null) {
						k = globals.formatProtId(line.substring(1));
					}
					else { // record current protein before continuing
						
						len = seq.trim().length();
						protLen.put( k, len);
						
						k = null;
						k = "";
						k = globals.formatProtId(line.substring(1));
						
						seq = null;
						len = -1;
					}
					seq = "";
				}
				else seq += line.trim();
			}
			br.close();
			
			
			if( (seq != null) && !seq.isEmpty() ) { // record the last protein
				len = seq.trim().length();
				protLen.put(k, len);
			}
			
			if(everythingOK == false) { return true; } // get out of function
		}
		catch (Exception e) {
			//print nothing
		}

        return false;
	}

	
	
	// Function assigns command line arguments to global variables
	public static void parseCommandLineArgs(String[] argv) throws IOException {
		boolean a = false;
		boolean b = false;
		
		if(argv[0].equals("-p")) paramFile = argv[1];
		
		if(paramFile == null || paramFile.isEmpty() ) {
			printError(paramFileNull);
			System.exit(0);
		}
		
		fastaFile = ""; // initalize the variable to be empty
		parseParametersFile();
		
		/*
		 * Check to make sure all of these options were set
		 */
		if(DBname == null) printError(noDBname);
		if(combinedFile == null) printError(noCombinedTag);
		if(srcDir == null) printError(noSrcDir);
		if(decoyTag == null) printError(noDecoyTag);
		//if(fastaFile == null) printError(noFastaFile);
		if(maxIniProbTH == -100) printError(noMaxIniProbTH);
		if(iniProbTH == -100) printError(noIniProbTH);
		if(minCombinedFilePw == -100) printError(noMinCombinedFilePw);
		//if(minPw == -100) printError(noMinPw);
		
		if(gene2protFile == null && byGene) printError(mapFileNotFound);
	}

	
	
	/*
	 * Function assigns the operating system type to OStype
	 */
	public static void getOStype() {
				
		// determine the operating system being used
		String OS = System.getProperty("os.name").toLowerCase();
		
		if( OS.indexOf("win") >= 0) OStype = "windows";
		if( OS.indexOf("mac") >= 0) OStype = "mac";
		else OStype = "nix";
	}
	
	
	/*
	 * Function to parse the parameters file and record its values
	 */
	public static void parseParametersFile() {
		
		File inputFile = new File(paramFile);
		if(!inputFile.exists()) printError(paramFileNotFound);
		Pattern regex = Pattern.compile("\\s+"); // for matching white space in ary[1] element
		
		try {
			BufferedReader input = new BufferedReader( new FileReader(inputFile) );
			String line = null;
			String ary[] = new String[ 2 ]; // key, value
			
			while( (line = input.readLine()) != null ) {
				if(line.startsWith("#")) continue; // ignore comment lines
				if( line.matches("^[^\\w]*") ) continue; // blank line

				// this just handle's error where user didn't provide a field after the '=' sign
				if( line.split("=").length < 2 ) {
					ary[0] = line.split("=")[0];
					ary[1] = "ERROR";
				}
				else ary = line.split("=");

				
				// special handling of the required AA mods text filed
				if(ary[0].trim().equals("reqAAmods")) {
					String tmp = ary[1].trim().replaceAll("\\s", "");
					if(tmp.length() < 3) ary[1] = "ERROR";
					else ary[1] = tmp;
				}
				else if( ary[0].trim().contains("Prob") || ary[0].trim().contains("Pw") ) {
					// different handling for fields that are suppose to be numeric
					ary[1].trim(); // just keep the value as is (could be negative number)
				}
				else {
					// just take the first 'word' in ary[1]
					ary[1].trim();
					String tmp[] = regex.split(ary[1]);
					ary[1] = tmp[0];
				}
				
				if(ary[0].equals("dbName")) DBname = (ary[1].equals("ERROR")) ? "ABACUS" : ary[1];
				if(ary[0].equals("combinedFile")) combinedFilePath = (ary[1].equals("ERROR")) ? "" : ary[1];
				if(ary[0].equals("srcDir")) srcDir = (ary[1].equals("ERROR")) ? "" : ary[1];
				if(ary[0].equals("fasta")) fastaFile = (ary[1].equals("ERROR")) ? "" : ary[1];
				if(ary[0].equals("decoyTag")) decoyTag = (ary[1].equals("ERROR")) ? "" :ary[1];
				
				if(ary[0].equals("maxIniProbTH")) maxIniProbTH = (ary[1].equals("ERROR")) ? 0.99 : Double.parseDouble(ary[1]);
				if(ary[0].equals("iniProbTH")) iniProbTH = (ary[1].equals("ERROR")) ? 0.50 : Double.parseDouble(ary[1]);
				if(ary[0].equals("epiTH")) epiThreshold = (ary[1].equals("ERROR")) ? 0.50 : Double.parseDouble(ary[1]);
				if(ary[0].equals("minCombinedFilePw")) minCombinedFilePw = (ary[1].equals("ERROR")) ? 0.90 : Double.parseDouble(ary[1]);
				
				//if(ary[0].equals("minPw")) minPw = (ary[1].equals("ERROR")) ? 0 : Double.parseDouble(ary[1]);

				if(ary[0].equals("verboseResults")) makeVerboseOutput = (ary[1].equals("true")) ? true : false;

				if(ary[0].equals("outputFile")) outputFilePath = (ary[1].equals("ERROR")) ? "ABACUS_output.tsv" : ary[1];

				if(ary[0].equals("recalcPepWts")) recalcPepWts = (ary[1].equals("true")) ? true : false;
				
				// these two options are legacy code that we keep here just in case we need it again in the future
				// technically they should not be needed. 
				if(ary[0].equals("protXMLsuffix")) protXMLsuffix = (ary[1].equals("ERROR")) ? "prot.xml" : ary[1];
				if(ary[0].equals("pepXMLsuffix")) pepXMLsuffix = (ary[1].equals("ERROR")) ? "pep.xml" : ary[1];
				
				if(ary[0].equals("keepDB")) {
					if(ary[1].equals("true")) keepDB = true;
					else { keepDB = false; }
                }

				if(ary[0].equals("asNSAF")) {
					if(ary[1].equals("true")) doNSAF = true;
					else { doNSAF = false; }
				}


				if(ary[0].equals("reqAAmods")) pepRegexText = (ary[1].equals("ERROR")) ? "" : ary[1];

				/*
				 * Some quick error checking, you can't have iniProbTH be greater
				 * that maxIniProbTH.
				 */
				if(iniProbTH > maxIniProbTH) {
					double tmp1 = maxIniProbTH;
					maxIniProbTH = iniProbTH;
					iniProbTH = tmp1;
				}

				if(fastaFile == null) fastaFile = "";
				
				
				if( fastaFile.isEmpty() && doNSAF) {
					System.err.printf("\nERROR: NSAF output requires a FASTA file.\n\n");
					System.exit(-1000);
				}
				
				
				
				
				//determine output type
				if(globals.outputFormat == -1) { // -1 means 'unset yet'
					if(ary[0].equals("output")) {
						if(ary[1].equals("Custom")) outputFormat = customOutput;
						else if(ary[1].equals("GeneQspec")) outputFormat = geneQspecFormat;
						else if(ary[1].equals("ProtQspec")) outputFormat = protQspecFormat;
						else if(ary[1].equals("Default")) outputFormat = defaultOutput;
						else if(ary[1].equals("Peptide")) {
							byPeptide = true;
							outputFormat = peptideOutput;
						}
						else if(ary[1].equals("Gene")) {
							byGene = true;
							outputFormat = geneOutput;
						}
						
					}
					
					if( 
						((outputFormat == geneQspecFormat) || (outputFormat == protQspecFormat)) &&
						(fastaFile.isEmpty())
					) {
						System.err.println("\nERROR: QSpec output requires a FATA file.\n");
						System.exit(-1000);
					}
				}

				if(outputFormat == customOutput) {
					//This code deals with generating custom outputs
					if(ary[0].equals("printC")) parseCustomOutputOptions(ary);
					if(ary[0].equals("printE")) parseCustomOutputOptions(ary);

					if( (outputFormat == customOutput) &&
							(printE == null || printC == null)
					  ) {
							System.err.println("printC: " + printC.size());
							System.err.println("printE: " + printE.size());

							System.err.printf("\nERROR: You indicated you wanted a custom output");
							System.err.printf(" but I was unable to interpret the options for either ");
							System.err.printf("'printC' or 'printE'.\nPlease check your param_file and try again.\n\n");
							System.exit(-1000);
					}
				}


				if(outputFormat == defaultOutput) byGene = false; // ensures protein-centric output

				
				/*
				 * Path to to text file mapping protids to geneSymbols.
				 * This option implies the user has provided a file that maps
				 * protids to their gene symbols. This file must be tab delimited
				 * with the format gene ID <TAB> protein ID.
				 */
				if(ary[0].equals("gene2prot")) gene2protFile = ary[1];
			}
			input.close();


			if(globals.epiThreshold == -1) globals.epiThreshold = 0;

			// extract the combined file name from combinedFilePath
			if( globals.combinedFilePath != null ) {
				File cf = new File(globals.combinedFilePath);
				combinedFile = cf.getName();
				cf = null;
			}


			/*
			 * If the user provided a gene2prot mapping file AND elected 
			 * to generate GeneQspec data, we assume the user wants gene-centric
			 * Qspec data.
			 */
			if( (gene2protFile != null) && (outputFormat == geneQspecFormat) ) byGene = true;
			

			// Handle case where user provides required peptide modification features
			if( (pepRegexText != null) && (!pepRegexText.isEmpty()) ) formatPepRegex();
			
		} catch (FileNotFoundException e) {
			System.err.print("\n#\n#FileNotFoundException Error at: globals.parseParameterFile()\n#\n#\n\n");
			System.err.print(e.toString());
		} catch (IOException e) {
			System.err.print("\n#\n#IOException Error at: globals.parseParameterFile()\n#\n#\n\n");
			System.err.print(e.toString());
		}
		
		// If the user didn't provide a decoyTag, make up one that won't occur
		// in the data. This allows the code to run and labels all proteins as
		// being forward sequences.
		if(decoyTag == null) decoyTag = UUID.randomUUID().toString().replace('-', 'x');
	}
	
	
	/*
	 * Function to print the parameters that will be used for this run of Abacus
	 */
	public static String printParameters() {
		String ret = "\n\nParameters for this execution:\n"
				   + "\tSource directory: '" + srcDir + "'\n"
				   + "\tDB name:          '" + DBname + "'\n"
				   + "\tOutput file:      '" + outputFilePath + "'\n"
				   + "\tCombined file P:   " + minCombinedFilePw + "\n"
				   + "\tiniProb threshold: " + iniProbTH + "\n"
				   + "\tmaxIniProb:        " + maxIniProbTH + "\n"
				   + "\tKeep DB files:     " + keepDB + "\n"
				   + "\tRecalc Pep Wts:    " + recalcPepWts + "\n";


		// code used print output type
		String outputTxt = null;
		switch (outputFormat) {
			case 1: outputTxt = "Protein Qspec"; break;
			case 2: outputTxt = "Custom"; break;
			case 3: outputTxt = "Gene"; break;
			case 4: outputTxt = "Gene Qspec"; break;
			case 5: outputTxt = "Peptide"; break;
			default: outputTxt = "Default"; break;
		}

		ret += "\tOutput format:     " + outputTxt + "\n";


		if( (null != pepMods_plus) && (pepMods_plus.length > 0) ){
			String x = "";
			for(int i = 0; i < pepMods_plus.length - 1; i++) { x += pepMods_plus[i] + ", "; }
			x += pepMods_plus[ (pepMods_plus.length - 1) ];
			ret += "\tAA mods to keep:   " + x + "\n";
		}

		if( (null != pepMods_minus) && (pepMods_minus.length > 0) ) {
			String x = "";
			for(int i = 0; i < pepMods_minus.length - 1; i++) { x += pepMods_minus[i] + ", "; }
			x += pepMods_minus[ (pepMods_minus.length - 1) ];
			ret += "\tAA mods to avoid:  " + x + "\n";
		}

		ret += "\n";

		return ret;
	}

	/*
	 * Function generates the current time as a string.
	 */
	public static String formatCurrentTime() {
		String ret = null;
		
		String[] monthName = {
				"Jan", "Feb", "Mar", 
				"Apr", "May", "Jun", 
				"Jul", "Aug", "Sep", 
				"Oct", "Nov", "Dec"
		};
		
		Calendar rightNow = Calendar.getInstance();
		String month = monthName[ rightNow.get(Calendar.MONDAY)];
		
		int min = rightNow.get(Calendar.MINUTE);
		String minString = null;
		if(min < 10) minString = "0" + Integer.toString(min);
		else minString = Integer.toString(min);
		
		int hour = rightNow.get(Calendar.HOUR_OF_DAY);
		String hourString = null;
		if(hour < 10) hourString = "0" + Integer.toString(hour);
		else hourString = Integer.toString(hour);
		
		int day = rightNow.get(Calendar.DAY_OF_MONTH);
		String dayString = null;
		if(day < 10) dayString = "0" + Integer.toString(day);
		else dayString = Integer.toString(day);
		
		ret = "" + rightNow.get(Calendar.YEAR)
		    + month
		    + dayString 
		    + "_"
		    + hourString
		    + minString;
		
		return ret;
	}

	
	public static void printError(int err) {
		
		switch (err) {
			case 1: System.err.println("\nError in " + paramFile + ": dbName=?\n"); break;
			case 2:	System.err.println("\nError in " + paramFile + ": combinedTag=?\n"); break;
			case 3:	System.err.println("\nError in " + paramFile + ": srcDir=?\n"); break;
			case 4: System.err.println("\nError in " + paramFile + ": decoyTag=?\n"); break;
			case 5: System.err.println("\nError in " + paramFile + ": maxIniProbTH=?\n"); break;
			case 6: System.err.println("\nError in " + paramFile + ": iniProbTH=?\n"); break;
			case 7: System.err.println("\nError in " + paramFile + ": minCombinedFilePw=?\n"); break;
			case 8: System.err.println("\nError in " + paramFile + ": minPw=?\n"); break;
			case 9: System.err.println("\nError: srcdir='" + srcDir + "' was not found.\n"); break;
			case 10: System.err.println("\nError in " + paramFile + ": fasta=?\n"); break;
			case 11: System.err.println("\nError: fastaFile='" + fastaFile + "' was not found.\n"); break;
			case 12: System.err.println("\nError: paramFile='" + paramFile + "' was not found.\n"); break;
			case 13: System.err.println("\nError: gene2protFile='" + gene2protFile + "'. You didn't specify a gene-to-protein ID mapping file\n"); break;
			case 14: System.err.println("\nError: No parameter file was read in. Did you forget the '-p' option?\n"); break;
			case 15: System.err.println("\nError: The path for the output file does not exist.\n"); break;
			default: System.err.println("Undefined error."); break;
		}
		
		System.exit(-1000);
	}
	
	/**************
	 * 
	 * Function parses the 'printX' options from the param_file and records
	 * what columns the user wants reported in the final output.
	 *
	 */
	private static void parseCustomOutputOptions(String[] ary) {
		
		if(ary[0].equals("printC")) printC = new HashSet<String>();
		if(ary[0].equals("printE")) printE = new HashSet<String>();
		
		
		String opts[] = ary[1].split(",");
		for(int i = 0; i < opts.length; i++) {
			
			// options available to printE
			if( ary[0].equals("printE") ) {
				if(opts[i].equals("id")) printE.add("_id"); 
				if(opts[i].equals("Pw")) printE.add("_Pw");
				if(opts[i].equals("numPepsTot")) printE.add("_numPepsTot");
				if(opts[i].equals("numPepsUniq")) printE.add("_numPepsUniq");
				if(opts[i].equals("numSpecsTot")) printE.add("_numSpecsTot");
				if(opts[i].equals("numSpecsUniq")) printE.add("_numSpecsUniq");
				if(opts[i].equals("numSpecsAdj")) printE.add("_numSpecsAdj");

				if(globals.doNSAF) {
					if(opts[i].equals("numSpecsTot")) printE.add("_totNSAF");
					if(opts[i].equals("numSpecsUniq")) printE.add("_uniqNSAF");
					if(opts[i].equals("numSpecsAdj")) printE.add("_adjNSAF");
				}
			}
			
			// options available to printC
			if(ary[0].equals("printC")) {
				if(opts[i].equals("id")) printC.add("ALL_ID");
				if(opts[i].equals("allPw")) printC.add("ALL_PW");
				if(opts[i].equals("localPw")) printC.add("ALL_LOCALPW");
				if(opts[i].equals("numPepsTot")) printC.add("ALL_NUMPEPSTOT");
				if(opts[i].equals("numPepsUniq")) printC.add("ALL_NUMPEPSUNIQ");
				if(opts[i].equals("numSpecsTot")) printC.add("ALL_NUMSPECSTOT");
				if(opts[i].equals("numSpecsUniq")) printC.add("ALL_NUMSPECSUNIQ");
				
				if(opts[i].equals("maxPw")) printC.add("MAXPW");
				if(opts[i].equals("maxIniProb")) printC.add("MAXINIPROB");
				if(opts[i].equals("wt_maxIniProb")) printC.add("WT_MAXINIPROB");
				if(opts[i].equals("maxIniProbUniq")) printC.add("MAXINIPROBUNIQ");
				
				if(opts[i].equals("protid")) printC.add("PROTID");
				if(opts[i].equals("isFwd")) printC.add("ISFWD");
				if(opts[i].equals("defline")) printC.add("DEFLINE");
				if(opts[i].equals("numXML")) printC.add("NUMXML");
				if(opts[i].equals("protLen")) printC.add("PROTLEN");
				
				if(opts[i].equals("geneid")) {
					// make sure the data will contain a gene id field
					if(gene2protFile != null) printC.add("GENEID");
				}
			}
		}
	}
	
	
	/********************************
	 * 
	 * This just prints a progress report to the terminal to let you know the
	 * computer is doing something.
	 * 
	 */
	public static void cursorStatus(int i, String msg) throws Exception {
		String anim = "|/-\\";
		int r = i % anim.length();

		String data = "\r" + msg + "  [ " + anim.charAt(r) + " " + i + " Working... ]";
		System.err.print(data);
	}
	
	
	
	/************************
	 * 
	 * Function to strip '#' character from strings. We need to do this because
	 * the replaceAll() function of Strings doesn't seem to work
	 * 
	 */
	public static String replaceAll(String src, char badChar, char goodChar) {
		String ret = null;
		StringBuilder sb = new StringBuilder();
		
		for(int i = 0; i < src.length(); i++) {
			if(src.charAt(i) == badChar) sb.append(goodChar);
			else sb.append(src.charAt(i));
		}
		
		ret = sb.toString();
		return ret;
	}


	/*****************
	 * 
	 * Function extracts the pepXML files used to create the given protXML file
	 * Data is stored in pepxmlTagHash and protxmlTagHash
	 * 
	 * @param xmlStreamReader
	 * @param protXMLfile 
	 * 
	 */
	public static boolean parseProtXML_header(XMLStreamReader xmlStreamReader, String protXMLfile, abacus_textArea console) {
		
		String[] ary = null;
		String tag = null;
		String origTag = null;
		String tmp = null;
		boolean status = false;

		
		// the combined file is never parsed for this information since it
		// will always contain all of the pepXML files
		if( !protXMLfile.contains(combinedFile) ) {
			
			// this code extracts the protXML files unique identifier
			Pattern protPattern = Pattern.compile("interact-(.+).prot.xml");
			Matcher protMatcher = protPattern.matcher(protXMLfile);
			if(protMatcher.find()) {
				origTag = protMatcher.group(1);
				
				// hyperSQL database cannot handle column names with hyphens or 
				// first characters that are digits. This could prevents such
				// problems.
				if(origTag.matches("^\\d.*")) {
					tmp = "x" + origTag;
					tag = globals.replaceAll(tmp, '-', '_');
				}
				else { 
					tag = globals.replaceAll(origTag, '-', '_'); 
				}
					
			}
			else {
				// the user's files don't follow the pattern: interact-<TAG>.prot.xml'
				// so we make one up based upon the given file name.
				// take the name of the file up to the first dot character
				int dotIdx = protXMLfile.indexOf('.');
				tag = globals.replaceAll(protXMLfile.substring(0, dotIdx), '-', '_');
				origTag = tag;
			}
			protTagHash.put(protXMLfile, tag);

			// check to see if we have a pepXML file with the same 'tag' as this
			// protXML file. If we do, we don't have to parse the protXML header line.
			String pXML = null;
			if( !search_srcDir_for_pepXML(origTag, tag) ) { // we have to parse the protXML file header line
			
				// this code identifies all the pepXML files that went into making
				// the current protXML file and associates them with the current protXML's tag
				//Pattern regexPattern = Pattern.compile(".*[\\/](.+.pep.xml)$");
				Pattern regexPattern = Pattern.compile(".*[\\/](.+." + pepXMLsuffix +")$");
				
				Matcher matcher = null;
				String err;

				for(int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
					String n = xmlStreamReader.getAttributeLocalName(i);
					String v = xmlStreamReader.getAttributeValue(i);

					if(n.equals("source_files")) {
						ary = v.split("\\s+"); // each element should be 1 pepXML file

						for(int j = 0; j < ary.length; j++) {
							matcher = regexPattern.matcher(ary[j]);
							if(matcher.find()) { //pepXML file extracted
								pXML = matcher.group(1);
								pepTagHash.put(pXML, tag);
							}
						}
					}
				}
			}
			
			// check to make sure that there there is at least one distinct
			// pepXML file for each protXML file
			if(pepTagHash.size() < protTagHash.size()) status = true;
		}
		return status; // false means everything is okay
	}



	/*******************
	 * Function searches the srcDir for the given pepXML file name.
	 * If it can't find it. The function then searches for a pepXML file with the
	 * same name as the protXML file.
	 * If that fails, the function returns false.
	 * 
	 * @param origProtXML_tag protXML_tag
	 * @return
	 */
	private static boolean search_srcDir_for_pepXML(String origProtXML_tag, String protXML_tag) {

		boolean ret = false;
		File dir = new File(globals.srcDir);
		String alt_pepXML = null;

		if(origProtXML_tag.equals(protXML_tag)) {
			alt_pepXML = "interact-" + protXML_tag + "."  + pepXMLsuffix;
		}
		else {
			// the original protXML tag started with a number. We padded the
			// beginning of the tag with an 'x' to may HyperSQL happy
			alt_pepXML = "interact-" + origProtXML_tag + "." + pepXMLsuffix;
		}

		if( globals.pepXmlFiles.contains(alt_pepXML) ) {
			pepTagHash.put(alt_pepXML, protXML_tag);
			ret = true;
		}
		
		return ret;
	}


	
	/****************
	 *  Function rounds a double to the desired number of decimal places
	 * @param d
	 * @param numDecimalPlaces
	 * @return
	 * 
	 */
	public static double roundDbl(double d, int numDecimalPlaces) {
		String pattern = "#.";
		double ret = 0;

		for(int i = 0; i < numDecimalPlaces; i++) {
			pattern += "#";
		}

		DecimalFormat df = new DecimalFormat(pattern);

		// in case the value of df is infinite (pos or neg) return 0
		try {
			ret = Double.valueOf(df.format(d));
		}catch(NumberFormatException e) {
			ret = 0.0;
		}
		
		return ret;
	}

	
	/***************
	 * Function takes the given string and returns a substring of it up to
	 * the first non-alphanumeric character it encounters
	 */
	public static String formatProtId(String line) {
		String ret = null;
		int p;
		char ch;

		// some default regex patterns to use for proteins
		Pattern uniprotPattern = Pattern.compile("^(sp|tr)\\|([^\\|]+).*");
		Matcher uniprotMatcher = uniprotPattern.matcher(line);
		
		Pattern refseqPattern = Pattern.compile("^gi\\|\\d+\\|ref\\|([^\\|]+).*");
		Matcher refseqMatcher = refseqPattern.matcher(line);
		
		Pattern ipiPattern = Pattern.compile("^IPI:([^\\|]+).*");
		Matcher ipiMatcher = ipiPattern.matcher(line);
		
		ret = "";
		
		if(uniprotMatcher.find()) { // this record is a uniprot match
			
			if(line.startsWith(decoyTag)) ret = decoyTag;
			ret += uniprotMatcher.group(2);
			return ret;
		}
		
		if(refseqMatcher.find()) { // this record is a refseq match
			
			if(line.startsWith(decoyTag)) ret = decoyTag;
			ret += refseqMatcher.group(1);
			return ret;
		}
		
		if(ipiMatcher.find()) { // this record is an IPI match
			
			if(line.startsWith(decoyTag)) ret = decoyTag;
			ret += ipiMatcher.group(1);
			return ret;
		}
		
		
		// use this code when 'line' doesn't conform to any of the regex's above.
		p = 0;
		for(int i = 0; i < line.length(); i++) {
			p++;
			ch = line.charAt(i);
			if(ch == ' ') break;
		}
		// if the protein defline is really, long take name up to
		// first non-alphanumeric character
		if(p >= 100) {
			p = 0;
			for(int i = 0; i < line.length(); i++) {
				p++;
				ch = line.charAt(i);
				if( !Character.isLetterOrDigit(ch) ) break;
				if(i >= 100) break;
			}
		}
		String tmpId = line.substring(0, p).trim();
		
		ret = line.substring(0, p).trim();
		return ret;
	}


	// Function returns a string reporting how long a given step in the program
	// took to run
	public static String formatTime(long elapsed_time) {
		String ret = null;
		
		double seconds = Math.floor((elapsed_time / 1000));
		double minutes = Math.floor((elapsed_time / (60*1000) ));
		double hours   = Math.floor((elapsed_time / (60*60*1000) ));

		if(seconds > 60) {
			double x =  Math.floor((seconds / 60));
			minutes += x;
			x = seconds % 60;
			seconds = x;
		}

		int hh = (int) hours;
		int mm = (int) minutes;
		int ss = (int) seconds;

		ret = Integer.toString( hh ) + ":";

		if(mm < 10) ret += "0";
		ret += Integer.toString( mm ) + ":";

		if(ss < 10) ret += "0";
		ret += Integer.toString( mm );
		
		return ret;
	}


	/************
	 *
	 * Function parses the pepRegexText value and stores the modifications
	 * listed in it as array elements in pepMods_plus, pepMods_minus
	 */
	public static void formatPepRegex() {

		String localregex = pepRegexText.replaceAll("\\s", "");
		String[ ] x = localregex.split(";");
		String curTxt = null;
		int p = 0, m = 0;
		
		
		int N = x.length;
		int nP = countChar(localregex, '+');
		int nM = countChar(localregex, '-');
		
		pepMods_minus = new String[ nM ];
		pepMods_plus = new String[ nP ];
		
		for(int i = 0; i < x.length; i++) {
			curTxt = x[i].trim();
			
			if( curTxt.startsWith("-") ) { // we want to avoid this modification
				if(null == pepMods_minus) pepMods_minus = new String[ nM ];

				pepMods_minus[ m ] = curTxt.substring(1).toUpperCase();
				m++;
			}
			else if( curTxt.startsWith("+") ) { // we want to keep this modification
				if(null == pepMods_plus) pepMods_plus = new String[ nP ];

				pepMods_plus[ p ] = curTxt.substring(1).toUpperCase();
				p++;
			}
		}
	}


	// Function returns true if the given modPeptide should be included in the
	// database. The decision for inclusion is done based upon the final score.
	// For every matched modification from pepMods_plus score goes up by 1
	// For every matched modification from pepMods_minus score goes down by 1
	// At the end if the score > 0, we return true
	public static boolean check_modPeptide(String modPep) {
		boolean status = false;
		int score;

		if( (null == pepMods_plus)  && (null == pepMods_minus) ) status = true;
		else {
			score = 0;
			if( (pepMods_plus != null) && (pepMods_plus.length > 0) ) {
				for(int i = 0; i < pepMods_plus.length; i++)
					if( modPep.contains( pepMods_plus[i] ) ) score++;
			}

			if( (pepMods_minus != null) && (pepMods_minus.length > 0) ) {
				for(int i = 0; i < pepMods_minus.length; i++)
					if( modPep.contains( pepMods_minus[i] ) ) score--;
			}

			// the current modPep value contains more 'plus' modifications
			// than negative ones so we return true
			if( score > 0 ) status = true;
			else if( score == 0 ) status = false;
		}
		
		return status;
	}

	
	
	// function returns the frequency of the given character in the given string
	public static int countChar(String haystack, char needle) {
		int ret = 0;
		
		for(int i = 0; i < haystack.length(); i++) {
			if(haystack.charAt(i) == needle) ret++;
		}
		
		return ret;
	}
	
	
	// Function records the pepXML file names and their tags
	public static void recordPepXMLtags() {
		String curFile = null;
		String curTag = null;
		Pattern pat = Pattern.compile("(.+)\\." + pepXMLsuffix + "$");
		Matcher matcher = null;
		
		for(int i = 0; i < pepXmlFiles.size(); i++) {
			curFile = pepXmlFiles.get(i);
			matcher = pat.matcher(curFile);
			if(matcher.find()) {
				curTag = matcher.group(1).replaceAll("interact-", "");
				pepTagHash.put(curFile, curTag);
			}
		}
	}
	
}

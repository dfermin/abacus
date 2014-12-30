package abacus;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.HashMap;
import javax.xml.stream.XMLStreamReader;

/*
 * Class defining the pepXML object
 */
public class pepXML {
	//private String searchEngine;   // stores the name/type of search tool used
	private String srcFile;
	private String specId;
	private int hitRank;
	private double mass;
	private int charge;
	private String peptide;
	private char prevAA;
	private char nextAA;
	private String modPeptide;
	private double iniProb;
	private boolean isProphetData;
	private double wt; // variable used in protXML files
	private double nsp; // variable used in protXML files
	private int ntt; // variable used in protXML files
	private int nspecs; // variable used in protXML files

	private HashMap<Integer, Integer> aaMods; // holds the AA modification
											  // positions
	
	
	// Variables for XTANDEM search results
	private double hyperscore;
	private double nextscore;
	private double xtandem_expect;
	
	// Variables for MASCOT search results
	private double mascot_ionscore;
	private double mascot_identityscore;
	private int mascot_star;
	private double mascot_homologyscore;
	private double mascot_expect;
	
	// Variables for SEQUEST search results
	private double sequest_xcorr;
	private double sequest_deltacn;
	private double sequest_deltacnstar;
	private double sequest_spscore;
	private double sequest_sprank;

	public pepXML() {
	}; // default constructor

	public pepXML(String txt, boolean b) {
		this.srcFile = txt;
		this.isProphetData = b;
		this.hitRank = -1; // once you have read in the top hit this goes to 1
	}

	// public SET functions (need them for parsing protXML files)
	public void setPeptide(String txt) {
		this.peptide = txt;
	}
	
	public void setCharge(String txt) {
		this.charge = Integer.parseInt(txt);
	}

	public void setIniProb(String txt) {
		this.iniProb = Double.parseDouble(txt);
	}

	public void setNSP(String txt) {
		this.nsp = Double.parseDouble(txt);
	}

	public void setWt(String txt) {
		this.wt = Double.parseDouble(txt);
	}

	public void setMass(String txt) {
		this.mass = Double.parseDouble(txt);
	}

	public void setNTT(String txt) {
		this.ntt = Integer.parseInt(txt);
	}

	public void setNspecs(String txt) {
		this.nspecs = Integer.parseInt(txt);
	}

	// public GET functions
	public String getSpecId() {
		return specId;
	}

	public double getMass() {
		return mass;
	}

	public int getCharge() {
		return charge;
	}

	public String getPeptide() {
		return peptide;
	}

	public char getPrevAA() {
		return prevAA;
	}

	public char getNextAA() {
		return nextAA;
	}

	public String getModPeptide() {
		return modPeptide;
	}

	public double getHyperscore() {
		return hyperscore;
	}

	public double getNextscore() {
		return nextscore;
	}

	public double getXtandem_expect() {
		return xtandem_expect;
	}

	public double getIniProb() {
		return iniProb;
	}

	public double getWt() {
		return wt;
	}

	public double getNSP() {
		return nsp;
	}

	public int getNTT() {
		return ntt;
	}

	public int getNspecs() {
		return nspecs;
	}

	// MASCOT variables
	public double getMascot_ionscore() {
		return mascot_ionscore;
	}
	
	public double getMascot_identityscore() {
		return mascot_identityscore;
	}
	
	public int getMascot_star() {
		return mascot_star;
	}
	
	public double getMascot_homologyscore() {
		return mascot_homologyscore;
	}
	
	public double getMascot_expect() {
		return mascot_expect;
	}
	
	
	
	// SEQUEST variables
	public double getSequest_xcorr() {
		return sequest_xcorr;
	}
	
	public double getSequest_deltacn() {
		return sequest_deltacn;
	}
	
	public double getSequest_deltacnstar() {
		return sequest_deltacnstar;
	}
	
	public double getSequest_spscore() {
		return sequest_spscore;
	}
	
	public double getSequest_sprank() {
		return sequest_sprank;
	}
	
	
	
	/*
	 * Function parses the given XML stream and records the relevant information
	 * found in it.
	 */
	public void parse_pepXML_line(XMLStreamReader xmlStreamReader) {
		String attrName = null;
		String attrValue = null;

		if(this.hitRank == 1) return; // this means we have already recorded the best hit for this PSM
		
		for (int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if (attrName.equals("hit_rank")) 
				this.hitRank = Integer.parseInt(attrValue);
			
			if (attrName.equals("spectrum"))
				this.specId = attrValue;
			if (attrName.equals("assumed_charge"))
				this.charge = Integer.parseInt(attrValue);
			if (attrName.equals("precursor_neutral_mass"))
				this.mass = Double.parseDouble(attrValue);

			if (attrName.equals("peptide"))
				this.peptide = attrValue;
			if (attrName.equals("peptide_prev_aa"))
				this.prevAA = attrValue.charAt(0);
			if (attrName.equals("peptide_next_aa"))
				this.nextAA = attrValue.charAt(0);
		}

		if(this.isProphetData) this.charge = 0;
	}

	/*
	 * Function parses amino acid modifications into aaMods variable
	 */
	public void record_AA_mod(XMLStreamReader xmlStreamReader) {
		String attrName = null;
		String attrValue = null;
		int k = -1;
		int v = 0;

		if (this.aaMods == null)
			this.aaMods = new HashMap<Integer, Integer>();

		for (int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if (attrName.equals("mod_nterm_mass")) { // N-terminal modification
				k = -100;
				v = 43;
				this.aaMods.put(k,v);
			}
			else { // not an N-terminal modification
				if (attrName.equals("position"))
					k = Integer.parseInt(attrValue) - 1;
				if (attrName.equals("mass")) {
					v = (int) Math.round(Double.parseDouble(attrValue));

					if (k > -1 && v > 0)
						this.aaMods.put(k, v);
					else {
						System.err.printf("\nERROR: mod_aminoacid_mass line pepXML::record_AA_mod()\n");
						System.err.println(this.specId + "\n");
						System.err.println(xmlStreamReader.toString());
						System.exit(-1);
					}
				}
			}
		}
	}

	/*
	 * Function parses search_score lines
	 */
	public void parse_search_score_line(XMLStreamReader xmlStreamReader) {
		String attrValue = null;

		for (int i = 0, j = 1; i < xmlStreamReader.getAttributeCount(); i++, j++) {
			attrValue = xmlStreamReader.getAttributeValue(i);

			/*
			 *   X!Tandem search scores
			 */
			if (attrValue.equals("hyperscore")) {
				//GLOBALS.searchEngineSet.add("XTANDEM");
				this.hyperscore = Double.parseDouble(xmlStreamReader
						.getAttributeValue(j));
			}
			if (attrValue.equals("nextscore"))
				this.nextscore = Double.parseDouble(xmlStreamReader
						.getAttributeValue(j));
			if (attrValue.equals("expect"))
				this.xtandem_expect = Double.parseDouble(xmlStreamReader
						.getAttributeValue(j));
			
			
			/*
			 *   Mascot search scores
			 */
			if (attrValue.equals("ionscore")) {
				//GLOBALS.searchEngineSet.add("MASCOT");
				this.mascot_ionscore = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			}
			if (attrValue.equals("identityscore"))
				this.mascot_identityscore = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("star"))
				this.mascot_star = Integer.parseInt(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("homologyscore")) 
				this.mascot_homologyscore = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("expect"))
				this.mascot_expect = Double.parseDouble(xmlStreamReader
						.getAttributeValue(j));
			
			
			/*
			 *   Sequest search scores
			 */
			if (attrValue.equals("xcorr")) {
				//GLOBALS.searchEngineSet.add("SEQUEST");
				this.sequest_xcorr = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			}
			if (attrValue.equals("deltacn"))
				this.sequest_deltacn = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("deltacnstar"))
				this.sequest_deltacnstar = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("spscore"))
				this.sequest_spscore = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
			if (attrValue.equals("sprank"))
				this.sequest_sprank = Double.parseDouble(xmlStreamReader.
						getAttributeValue(j));
		}
	}
	
	
	/*
	 * Function parses out peptide probability
	 */
	public void record_iniProb(XMLStreamReader xmlStreamReader) {
		String attrName = null;
		String attrValue = null;

		for (int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if (attrName.equals("probability"))
				this.iniProb = Double.parseDouble(attrValue);
		}
	}

	/*
	 * Function annotates modPeptide
	 */
	public void annotate_modPeptide() {
		if (this.aaMods == null)
			this.modPeptide = this.peptide; // peptide has no modifications
		else {
			modPeptide = "";
			// check for N-terminal modification
			if(this.aaMods.containsKey(-100)) {
				this.modPeptide += "n[43]";
				this.aaMods.remove(-100);
			}

			for (int i = 0; i < this.peptide.length(); i++) {
				this.modPeptide += this.peptide.charAt(i);
				if (this.aaMods.containsKey(i))
					this.modPeptide += "[" + this.aaMods.get(i) + "]";
			}
		}
	}

	public void write_to_db(PreparedStatement prep) {

		if( !globals.check_modPeptide(this.modPeptide) ) ; // do nothing (don't add it into database
		else if(this.iniProb > 0) {
			try {
				prep.setString(1, this.srcFile.toUpperCase());
				prep.setString(2, this.specId);
				prep.setInt(3, this.charge);
				prep.setString(4, this.peptide);
				prep.setString(5, this.modPeptide);
				prep.setDouble(6, this.iniProb);
				prep.addBatch();
				globals.proceedWithQuery = true; // true means you are implying the insertion worked

			} catch (SQLException e) {
				e.printStackTrace();
				System.exit(-1);

			}
			
		}
	}
}

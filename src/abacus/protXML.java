package abacus;

import java.sql.PreparedStatement;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import javax.xml.stream.XMLStreamReader;


/*
 * Class defining the protXML object
*/
public class protXML {
	private String srcFile;
	private int groupid;
	private double Pw;
	private double localPw;
	private String siblingGroup;
	private int isFwdGroup;
	private boolean isProphet;
	private HashMap<String, Integer> protIdClass; //holds protId and isFwd (0/1)
	private HashMap<String, String> protIds; // holds protIds and deflines
	private HashMap<String, Integer> protLen; // holds protIds and protein lengths
	private HashMap<String, pepXML> peptides; //holds modPeps and their pepXML objects

	public protXML() {}; // default constructor

	public protXML(String txt, boolean iProphet_status) {

		this.srcFile = txt;
		this.protIds = new HashMap<String, String>();
		this.peptides = new HashMap<String, pepXML>();
		this.protIdClass = new HashMap<String, Integer>();
		this.protLen = new HashMap<String, Integer>();
		this.isFwdGroup = 0;
		this.isProphet = iProphet_status; // true means the file is an i-Prophet output file
	}


	// public GET functions
	public int getGroupId() { return groupid; }
	public double getPw() { return Pw; }
	public double getLocalPw() { return localPw; }
	public String getSiblingGroup() { return siblingGroup; }
	public int getFwdStatus() { return isFwdGroup; }

	/*
	 *  Function to parse protXML file lines
	 */
	public void parse_protGroup_line(XMLStreamReader xmlStreamReader) {
		String attrName  = null;
		String attrValue = null;

		for(int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName  = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if(attrName.equals("group_number")) this.groupid = Integer.parseInt(attrValue);
			if(attrName.equals("probability")) this.Pw = Double.parseDouble(attrValue);
		}
	}


	/*
	 *  Function to parse protein line in protXML file
	 */
	public String parse_protein_line(XMLStreamReader xmlStreamReader) {
		String attrName  = null;
		String attrValue = null;
		String protid_  = null;

		for(int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName  = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if(attrName.equals("protein_name")) protid_ = globals.formatProtId(attrValue);
			if(attrName.equals("probability")) this.localPw = Double.parseDouble(attrValue);
			if(attrName.equals("group_sibling_id")) this.siblingGroup = attrValue;
		}

		return(protid_);
	}


	/*
	 *  Function records the current protein ID and it's description
	 */
	public void setProtId(String defline, String protid_) {
		// we limit protId lines to 500 characters.
		String x = defline.replace('\'', '_');
		defline = x;
		if(defline.length() == 0) this.protIds.put(protid_.trim(), "No Description");
		else if(defline.length() > 500) this.protIds.put(protid_.trim(), defline.substring(0, 500));
		else this.protIds.put(protid_.trim(), defline);
	}


	/*
	 *  Function records the current peptide into protXML object
	 */
	public String parse_peptide_line(XMLStreamReader xmlStreamReader) {
		pepXML curPep = null;
		curPep = new pepXML();
		String attrName  = null;
		String attrValue = null;
		int pepCtr = 0;
		
		for(int i = 0; i < xmlStreamReader.getAttributeCount(); i++) {
			attrName  = xmlStreamReader.getAttributeLocalName(i);
			attrValue = xmlStreamReader.getAttributeValue(i);

			if(attrName.equals("peptide_sequence")) curPep.setPeptide(attrValue);
			if(attrName.equals("charge")) curPep.setCharge(attrValue);
			if(attrName.equals("initial_probability")) curPep.setIniProb(attrValue);
			if(attrName.equals("nsp_adjusted_probability")) curPep.setNSP(attrValue);
			if(attrName.equals("weight")) curPep.setWt(attrValue);
			if(attrName.equals("n_enzymatic_termini")) curPep.setNTT(attrValue);
			if(attrName.equals("n_instances")) curPep.setNspecs(attrValue);
			if(attrName.equals("calc_neutral_pep_mass")) curPep.setMass(attrValue);
		}

		if(this.isProphet) curPep.setCharge("0"); // special case for i-prophet data

		pepCtr = peptides.size(); // need this in case the same peptide sequence occurs twice with different charge states
		
		String k = null;
		k = "" + curPep.getPeptide() +
		    "-" + curPep.getMass() +
		    "-" + curPep.getCharge() +
		    "-" + curPep.getIniProb() +
			"-#" + Integer.toString(pepCtr);
		
		peptides.put(k, curPep);
		curPep = null;

		return(k);
	}



	/*
	 *  Function clears out variables for next protein block in group
	 */
	public void clear_variables() {
		this.protIds.clear();
		this.peptides.clear();
		this.protLen.clear();
		this.localPw = 0.0;
		this.siblingGroup = null;
	}



	/*
	 *  Function records the current peptide's modifications (if any)
	 */
	public void record_AA_mod_protXML(XMLStreamReader xmlStreamReader, String k) {
		pepXML curPep = peptides.get(k);
		curPep.record_AA_mod(xmlStreamReader);
		peptides.put(k, curPep);
		curPep = null;
	}



	/*
	 *  Function to annotate modPeptide string
	 */
	public void annotate_modPeptide_protXML(String k) {
		pepXML curPep = peptides.get(k);
		
		if(curPep.getModPeptide() == null) {
			curPep.annotate_modPeptide();
			peptides.put(k,curPep);
			curPep = null;
		}
	}


	/*
	 * The protein length is included now by default in protXML files
	 */
	public void record_protLen(String curProtid, String pl) {
		int protLen = 0;
		protLen = Integer.parseInt(pl);
		this.protLen.put(curProtid, protLen);
	}



	/*
	 *  Function determines if the group is a Forward or Decoy group.
	 *  If it is a forward group, all decoy matches in it will be removed.
	 *  (This program takes the optimistic viewpoint that if a group contains
	 *   a mix of forward and decoy proteins, the decoys are a fluke and should
	 *   be deleted).
	 */
	public void classify_group() {

		Set<String> protKeys = this.protIds.keySet();
		Iterator<String> protIter = protKeys.iterator();

		String k = "";
		int isFwd_ctr = 0;


		while(protIter.hasNext()) {
			k = (String) protIter.next();

			if( k.startsWith( globals.decoyTag ) ) this.protIdClass.put(k, 0);
			else {
				this.protIdClass.put(k, 1);
				isFwd_ctr++;
			}
		}


		if(isFwd_ctr > 0) {  // the group (overall) is a forward group
			this.isFwdGroup = 1;

			HashMap<String, String> newProtids = new HashMap<String, String>();
			HashMap<String, Integer> newProtidsClass = new HashMap<String, Integer>();


			protIter = null;
			protIter = protKeys.iterator();


			while(protIter.hasNext()) { // remove decoys from group
				k = (String) protIter.next();

				if( !k.startsWith(globals.decoyTag) ) {
					newProtids.put(k, protIds.get(k));
					newProtidsClass.put(k, 1);
				}
			}

			this.protIdClass.clear();
			this.protIdClass.putAll(newProtidsClass);
			this.protIds.clear();
			this.protIds.putAll(newProtids);

			newProtidsClass.clear(); newProtidsClass = null;
			newProtids.clear(); newProtids = null;
		}
	}



	/*
	 *  Function writes the collected data for the current protein group to
	 *  the hyperSQL database (RAWprotXML table)
	 */
	public void write_to_db(PreparedStatement prep) {

		Set<String> protKeys = this.protIds.keySet();
		Iterator<String> protIter = protKeys.iterator();

		String defline = "";
		String k = "";
		int protClass = 0;

		while(protIter.hasNext()) { // iterate over protein Ids
			k = (String) protIter.next();
			defline = this.protIds.get(k);

			if(this.protIdClass.containsKey(k)) protClass = this.protIdClass.get(k);
			else protClass = 0;

			/*
			 *  This code has to be here to iterate over each peptide for every
			 *  protein in the protIds hashMap
			 */
			Set<String> pepKeys = peptides.keySet();
			Iterator<String> pepIter = pepKeys.iterator();
			String pepId = "";

			pepXML curPep = null;

			while(pepIter.hasNext()) { // iterate over peptides
				pepId = (String) pepIter.next();
				curPep = null;
				curPep = this.peptides.get(pepId);

				if( !globals.check_modPeptide( curPep.getModPeptide() ) ) continue;

				if(curPep.getIniProb() < globals.iniProbTH) continue;
				try {
					prep.setString(1, this.srcFile.toUpperCase());
					prep.setInt(2, this.groupid);
					prep.setString(3, this.siblingGroup);
					prep.setDouble(4, this.Pw);
					prep.setDouble(5, this.localPw);
					prep.setString(6, k);
					prep.setInt(7, protClass);

					//peptide level data
					prep.setString(8, curPep.getPeptide() );
					prep.setString(9, curPep.getModPeptide() );
					prep.setInt(10, curPep.getCharge() );
					prep.setDouble(11, curPep.getIniProb() );
					prep.setDouble(12, curPep.getWt() );

					prep.setString(13, defline);
					prep.addBatch();
					globals.proceedWithQuery = true; // by incrementing this you imply the that insertion worked
				}
				catch (Exception e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
			
			pepIter = null;
		}
		protIter = null;
	}


}
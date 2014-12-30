package abacus;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.DecimalFormat;
import abacus_textArea.abacus_textArea;
import java.util.Iterator;



/******************
 *
 * This class inherits from hyperSQLObject and is used to run gene-centric
 * queries.
 *
 * @author dfermin
 *
 */
public class hyperSQLObject_gene extends hyperSQLObject {


	/*****************
	 * Function creates a combined table from the combined file that is based
	 * on genes not proteins.
	 *
	 * @param conn
	 * @throws Exception
	 *
	 */
	public void makeGeneCombined(Connection conn, abacus_textArea console) throws Exception {

		if(console != null) console.append("Creating gene-centric combined table (this can take a while)...\n");
		else System.err.print("Creating gene-centric combined table (this can take a while)...\n");

		Statement stmt  = conn.createStatement();
		Statement stmt2 = conn.createStatement();

		ResultSet rs  = null;
		ResultSet rs2 = null;
		PreparedStatement prep = null;

		String query = null;

		query = "CREATE CACHED TABLE geneCombined ( "
			  + "  geneid, isFwd, modPeptide " //, charge "
			  + ") AS ( "
			  + "SELECT gn.geneid, c.isFwd, c.modPeptide " //, c.charge "
			  + "FROM combined c, gene2prot gn "
			  + "WHERE c.protid = gn.protid "
			  + "AND c.isFwd = 1 "
			  + "GROUP BY gn.geneid, c.isFwd, c.modPeptide " //, c.charge "
			  + "ORDER BY gn.geneid ASC "
			  + ") WITH DATA ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("ALTER TABLE geneCombined ADD COLUMN max_local_Pw DECIMAL(8,6) BEFORE modPeptide");
		stmt.executeUpdate("ALTER TABLE geneCombined ADD COLUMN maxPw DECIMAL(8,6) BEFORE max_local_Pw");
		stmt.executeUpdate("ALTER TABLE geneCombined ADD COLUMN iniProb DECIMAL(8,6)");

		stmt.executeUpdate("CREATE INDEX gc_idx1 ON geneCombined(geneid)");
		stmt.executeUpdate("CREATE INDEX gc_idx2 ON geneCombined(modPeptide)"); //, charge)");


		/*
		 * MAX(Pw) and MAX(localPw) update
		 */
		if(console != null) console.append("  Updating maxPw\n");
		else System.err.print("  Updating maxPw\n");

		query = "CREATE CACHED TABLE t1_ ( "
			  + "  geneid VARCHAR(100), "
			  + "  maxPw DECIMAL(8,6), "
			  + "  max_localPw DECIMAL(8,6) "
			  + ")";
		stmt.executeUpdate(query);

		query = "INSERT INTO t1_ ( "
			  + "SELECT gn.geneid, MAX(c.Pw), MAX(c.localPw) "
			  + "FROM gene2prot gn, combined c "
			  + "WHERE gn.protid = c.protid "
			  + "AND c.isFwd = 1 "
			  + "GROUP BY gn.geneid "
			  + ")";
		stmt.executeUpdate(query);

		stmt.executeUpdate("CREATE INDEX t1_idx1 ON t1_(geneid)");

		query = "UPDATE geneCombined AS gc "
			  + "  SET (maxPw, max_local_Pw) = ( "
			  + "    SELECT maxPw, max_LocalPw "
			  + "    FROM t1_ "
			  + "    WHERE t1_.geneid = gc.geneid "
			  + ")";
		stmt.executeUpdate(query);


		/*
		 * MAX(iniProb) update
		 */
		if(console != null) console.append("  Updating Peptide Probabilities\n");
		else System.err.print("  Updating Peptide Probabilities\n");


		query = "CREATE CACHED TABLE t2_ ("
			  + "  geneid VARCHAR(200),"
			  + "  modPeptide VARCHAR(250), "
//			  + "  charge INT, "
			  + "  iniProb DECIMAL(8,6) "
			  +")";
		stmt.executeUpdate(query);

		query = "INSERT INTO t2_ ("
//			  + "  SELECT gc.geneid, gc.modPeptide, gc.charge, MAX(px.iniProb) "
			  + "  SELECT gc.geneid, gc.modPeptide, MAX(px.iniProb) "
			  + "  FROM geneCombined gc, pepXML px "
			  + "  WHERE gc.modPeptide = px.modPeptide "
//			  + "  AND gc.charge = px.charge "
			  + "  GROUP BY gc.geneid, gc.modPeptide "//, gc.charge "
			  + "  ORDER BY gc.geneid "
			  + ")";
		stmt.executeUpdate(query);
		stmt.executeUpdate("CREATE INDEX t2_idx1 ON t2_(geneid)");
		stmt.executeUpdate("CREATE INDEX t2_idx2 ON t2_(modPeptide)"); //, charge)");
		stmt.executeUpdate("CREATE INDEX t2_idx3 ON t2_(geneid, modPeptide)"); //, charge)");

		query = "UPDATE geneCombined AS gc "
			  + "  SET iniProb = ( "
			  + "    SELECT x.iniProb "
			  + "    FROM t2_ x "
			  + "    WHERE x.geneid = gc.geneid "
			  + "    AND x.modPeptide = gc.modPeptide "
//			  + "    AND x.charge = gc.charge "
			  + ")";
		stmt.executeUpdate(query);


		/*
		 *  Need to insert decoy proteins as a separate query. By definition,
		 *  decoy proteins have now genes, so we can't map them.
		 */
		if(console != null) console.append("  Accounting for decoy protein matches (if any)\n");
		else System.err.print("  Accounting for decoy protein matches (if any)\n");

		query = "INSERT INTO geneCombined ("
			  + "SELECT CONCAT('decoy-',c.groupid), c.isFwd, "
			  + "  MAX(c.Pw), MAX(c.localPw), c.modPeptide, "
			  //+ "  c.charge, MAX(c.iniProb) "
			  + "  MAX(c.iniProb) "
			  + "FROM combined c "
			  + "WHERE c.isFwd = 0 "
			  + "GROUP BY c.groupid, c.isFwd, c.modPeptide "//, c.charge "
			  + ") ";
		stmt.executeUpdate(query);


		/*
		 * Record how many protein groups matched each geneid
		 */
		stmt.executeUpdate("ALTER TABLE geneCombined ADD COLUMN numGroups INT DEFAULT 1");

		if(console != null) console.append("  Calculating gene id usage\n");
		else System.err.print("  Calculating gene id usage\n");

		query = "CREATE CACHED TABLE t3_("
			  + " geneid VARCHAR(100), "
			  + " freq INT "
			  + ")";
		stmt.executeUpdate(query);

		query = "INSERT INTO t3_ ("
			  + "  SELECT gn.geneid, COUNT(DISTINCT c.groupid) "
			  + "  FROM combined c, gene2prot gn "
			  + "  WHERE c.protid = gn.protid "
			  + "  GROUP BY gn.geneid "
			  + "  ORDER BY gn.geneid "
			  + ")";
		stmt.executeUpdate(query);
		stmt.executeUpdate("CREATE INDEX t3_idx1 ON t3_(geneid)");


		query = "UPDATE geneCombined AS gc "
			  + "  SET numGroups = ( "
			  + "    SELECT x.freq "
			  + "    FROM t3_ x "
			  + "    WHERE x.geneid = gc.geneid "
			  + ")";
		stmt.executeUpdate(query);

		/*
		 * Recompute peptide weights to be gene-centric instead of protein-centric
		 */
		DecimalFormat df = new DecimalFormat("#0.0000");
		double n = 0;

		stmt.executeUpdate("ALTER TABLE geneCombined ADD COLUMN wt DECIMAL(8,6) DEFAULT 0");

		query = "UPDATE geneCombined "
			  + "  SET wt = ? "
			  + "WHERE modPeptide = ? ";
		prep = conn.prepareStatement(query);

		if(console != null) console.append("  Adjusting peptide weights on gene basis\n");
		else System.err.print("  Adjusting peptide weights on gene basis\n");

		rs = stmt.executeQuery("SELECT DISTINCT modPeptide FROM geneCombined");
		while(rs.next()) {
			String modPep = rs.getString(1);
			query = "SELECT COUNT(DISTINCT geneid) "
				  + "FROM geneCombined "
				  + "WHERE modPeptide = '" + modPep + "' ";
			rs2 = stmt2.executeQuery(query);
			rs2.next();

			n = 1.0000 / (double)rs2.getInt(1);
			String newWT = df.format(n);
			n = Double.parseDouble(newWT);

			prep.setDouble(1, n);
			prep.setString(2, modPep);
			prep.addBatch();
		}
		conn.setAutoCommit(false);
		prep.executeBatch();
		conn.setAutoCommit(true);

		//just cleaning house
		stmt.executeUpdate("DROP INDEX t1_idx1");
		stmt.executeUpdate("DROP TABLE t1_");
		stmt.executeUpdate("DROP INDEX t2_idx1");
		stmt.executeUpdate("DROP INDEX t2_idx2");
		stmt.executeUpdate("DROP INDEX t2_idx3");
		stmt.executeUpdate("DROP TABLE t2_");
		stmt.executeUpdate("DROP INDEX t3_idx1");
		stmt.executeUpdate("DROP TABLE t3_");


		prep.close(); prep = null;
		rs.close(); rs = null;
		rs2.close(); rs2 = null;
		stmt.close(); stmt = null;
		stmt2.close(); stmt2 = null;
	}


	/************
	 *
	 * Function creates geneXML table from individual protXML files
	 *
	 */
	public void makeGeneXML(Connection conn, abacus_textArea console) throws Exception {
		if(console != null) console.append("\nCreating geneXML table\n");
		else System.err.print("\nCreating geneXML table\n");


		Statement stmt = conn.createStatement();
		ResultSet rs = null;
		String query = null;

		query = "CREATE CACHED TABLE geneXML ( "
			  + "  tag VARCHAR(250), "
			  + "  geneid VARCHAR(100), "
			  + "  isFwd INT, "
			  + "  maxPw DECIMAL(8,6), "
			  + "  max_localPw DECIMAL(8,6), "
			  + "  modPeptide VARCHAR(250), "
//			  + "  charge INT, "
			  + "  iniProb DECIMAL(8,6) "
			  + ")";
		stmt.executeUpdate(query);

		query = "INSERT INTO geneXML ( "
			  + "SELECT r.tag, gn.geneid, r.isFwd, MAX(r.Pw), MAX(r.localPw), "
//			  + " r.modPeptide, r.charge, MAX(px.iniProb) "
			  + " r.modPeptide, MAX(px.iniProb) "
			  + "FROM protXML r, gene2prot gn, pepXML px "
			  + "WHERE r.protId = gn.protid "
			  + "AND r.isFwd = 1 "
			  + "AND r.tag = px.tag "
			  + "AND r.modPeptide = px.modPeptide "
//			  + "AND r.charge = px.charge "
              + "AND px.iniProb >= " + globals.iniProbTH + " "
              + "AND r.iniProb >= " + globals.iniProbTH + " "
			  + "GROUP BY r.tag, gn.geneid, r.isFwd, r.modPeptide "//, r.charge "
			  + "HAVING max(r.localPw) > " + this.minPw + " "
			  + ") ";
		stmt.executeUpdate(query);

		// decoy matches
		query = "INSERT INTO geneXML ("
			  + "SELECT r.tag , CONCAT('decoy-',r.groupid), 0, MAX(r.Pw), "
//			  + "  MAX(r.localPw), r.modPeptide, r.charge, MAX(r.iniProb) "
			  + "  MAX(r.localPw), r.modPeptide, MAX(r.iniProb) "
			  + "FROM protXML r "
			  + "WHERE r.isFwd = 0 "
              + "AND r.iniProb >= " + globals.iniProbTH + " "
			  + "GROUP BY r.tag, r.groupid, r.modPeptide" //, r.charge "
			  + ") ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("ALTER TABLE geneXML ADD COLUMN wt DECIMAL(8,6) DEFAULT 0 NOT NULL");

		stmt.executeUpdate("CREATE INDEX gx_idx1 ON geneXML(tag, geneid)");
//		stmt.executeUpdate("CREATE INDEX gx_idx2 ON geneXML(modPeptide, charge)");
		stmt.executeUpdate("CREATE INDEX gx_idx3 ON geneXML(modPeptide)");
		stmt.executeUpdate("CREATE INDEX gx_idx4 ON geneXML(tag)");
		stmt.executeUpdate("CREATE INDEX gx_idx5 ON geneXML(geneid)");


		// if the epiThreshold > iniProbTH then we have to remove
		// all entries where the gene does not have at least 1 peptide with a
		// probability >= epiThreshold
		if(globals.epiThreshold > globals.iniProbTH) {
			query = "CREATE MEMORY TABLE x_ ("
				  + " tag, geneid, maxIniProb "
				  + ") AS ( "
				  + "SELECT tag, geneid, MAX(iniProb) "
				  + "FROM geneXML "
				  + "GROUP BY tag, geneid "
				  + ") WITH DATA ";

			stmt.executeUpdate(query);
			stmt.executeUpdate("CREATE INDEX x_1 ON x_(tag, geneid)");
			stmt.executeUpdate("CREATE INDEX x_2 ON x_(maxIniProb)");

			query = "SELECT * FROM x_ WHERE maxIniProb < " + globals.epiThreshold;
			rs = stmt.executeQuery(query);

			Statement stmt2 = conn.createStatement();
			while(rs.next()) {
				String tag = rs.getString(1);
				String gid = rs.getString(2);

				query = "DELETE FROM geneXML "
					  + "WHERE tag = '" + tag + "' "
					  + "AND geneid = '" + gid + "' ";
				stmt2.executeUpdate(query);
			}
			stmt2.close();
			rs.close();
			rs = null;
			stmt2 = null;

			stmt.executeUpdate("DROP INDEX IF EXISTS x_2");
			stmt.executeUpdate("DROP INDEX IF EXISTS x_1");
			stmt.executeUpdate("DROP TABLE IF EXISTS x_");
		}

		stmt.close();
		stmt = null;
		query = null;
	}



	/************
	 *
	 * Function adjusts weights of peptides on a gene-basis
	 *
	 */
	public void adjustGenePeptideWT(Connection conn, abacus_textArea console) throws Exception {
		Statement stmt  = conn.createStatement();
		Statement stmt2 = conn.createStatement();
		Statement stmt3 = conn.createStatement();
		Statement stmt4 = conn.createStatement();
		Statement stmt5 = conn.createStatement();
		PreparedStatement prep = null;
		String query = null;
		ResultSet rs = null;
		ResultSet rs2 = null;
		ResultSet rs3 = null;
		ResultSet rs4 = null;
		

		int N = 0, ctr = 1;

		if(console != null) console.append("  Adjusting peptide weights on gene basis\n");
		else System.err.print("  Adjusting peptide weights on gene basis\n");

		query = "SELECT COUNT(DISTINCT tag) FROM srcFileTags WHERE fileType = 'prot'";
		rs = stmt.executeQuery(query);
		rs.next();
		N = rs.getInt(1);

		String msg = "Computing peptide weights.";
		if(console != null) console.monitorBoxInit(N, msg);

		query = "SELECT DISTINCT tag FROM srcFileTags WHERE fileType = 'prot'";
		rs = stmt.executeQuery(query);

		while(rs.next()) {
			String tag = rs.getString(1);

			query = "CREATE CACHED TABLE gpwt_" + tag + " (" // make this 'MEMORY
				  + " modPeptide, numGenes "
				  + ") AS ( "
				  + "  SELECT modPeptide, COUNT(DISTINCT geneid) "
				  + "  FROM geneXML "
				  + "  WHERE tag = '" + tag + "' "
				  + "  GROUP BY modPeptide "
				  + ") WITH DATA ";
			stmt2.executeUpdate(query);


			query = "CREATE INDEX gpwt_" + tag + "_idx1 ON gpwt_" + tag + "(modPeptide)";
			stmt2.executeUpdate(query);

			query = "ALTER TABLE gpwt_" + tag + " ADD COLUMN wt DECIMAL(8,6) DEFAULT 0 NOT NULL";
			stmt2.executeUpdate(query);

			query = "UPDATE gpwt_" + tag + " "
				  + "  SET gpwt_"+tag+".wt = ROUND( (1.0000 / CAST(numGenes AS DECIMAL(16,6))), 4) ";
			stmt2.executeUpdate(query);

			
			query = "SELECT DISTINCT modPeptide "
				  + "FROM geneXML "
				  + "WHERE tag = '" + tag + "' ";
			rs2 = stmt3.executeQuery(query);
			
			while(rs2.next()) {
				String modPep = rs2.getString(1);
				
				query = "SELECT gpwt_" + tag + ".wt "
					  + "FROM gpwt_" + tag + " " 
					  + "WHERE modPeptide = '" + modPep + "' ";
				rs3 = stmt4.executeQuery(query);
				rs3.next();
				double wt_ = rs3.getDouble(1);
				
				query = "UPDATE geneXML "
					  + "  SET geneXML.wt = " + wt_ + " "
					  + "WHERE modPeptide = '" + modPep + "' ";
				stmt5.executeUpdate(query);
				
			}
					
			
//			query = "UPDATE geneXML "
//				  + "  SET geneXML.wt = ( "
//				  + "    SELECT gpwt_" + tag + ".wt "
//				  + "    FROM gpwt_" + tag + " "
//				  + "    JOIN geneXML ON gpwt_" + tag + ".modPeptide = geneXML.modPeptide "
//				  + ") ";
//			System.err.println("\n"+query+"\n");
//			stmt2.executeUpdate(query);
			//stmt2.executeUpdate("DROP INDEX IF EXISTS gpwt_" + tag + "_idx1");
			//stmt2.executeUpdate("DROP TABLE IF EXISTS gpwt_" + tag + " ");

			//stmt2.executeUpdate("UPDATE geneXML SET wt = 0 WHERE wt IS NULL");

			if(console != null) {
				console.monitorBoxUpdate(ctr);
				String x = "  " + tag + "\n";
				console.append(x);
			}
			ctr++;
		}

		stmt.close();  stmt  = null;
        stmt2.close(); stmt2 = null;
        stmt3.close(); stmt3 = null;
		if(console != null) {
			console.closeMonitorBox();
			console.append("\n");
        }
		else System.err.print("\n");
	}


	/**************
	 *
	 * Creates the geneidSummary table
	 *
	 */
	public void makeGeneidSummary(Connection conn, abacus_textArea console) throws Exception {

		String msg = "Creating geneidSummary table (this can take a while)...";
		if(console != null) console.append(msg + "\n");

		Statement stmt = conn.createStatement();
		Statement stmt2 = conn.createStatement();

		ResultSet rs  = null;
		ResultSet rs2 = null;

		String query = null;

		query = "CREATE CACHED TABLE geneidSummary ( "
			  + "  geneid VARCHAR(100), "
			  + "  isFwd INT, "
			  + "  maxPw DECIMAL(8,6), "
			  + "  max_localPw DECIMAL(8,6), "
			  + "  maxIniProb DECIMAL(8,6), "
			  + "  numGroups INT "
			  + ") ";
		stmt.executeUpdate(query);

		query = "INSERT INTO geneidSummary ( "
			  + "  SELECT geneid, isFwd, maxPw, max_local_Pw, MAX(iniProb), numGroups "
			  + "  FROM geneCombined "
			  + "  GROUP BY geneid, isFwd, maxPw, max_local_Pw, numGroups "
			  + "  ORDER BY geneid "
			  + ") ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("CREATE INDEX geneSum_idx1 ON geneidSummary(geneid)");

		stmt.executeUpdate("ALTER TABLE geneidSummary ADD COLUMN numXML INT");
		stmt.executeUpdate("ALTER TABLE geneidSummary ADD COLUMN numSpecsTot INT DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE geneidSummary ADD COLUMN numSpecsUniq INT DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE geneidSummary ADD COLUMN numPepsTot INT DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE geneidSummary ADD COLUMN numPepsUniq INT DEFAULT 0");


		int iter = 0; // for progress bar

		rs = stmt.executeQuery("SELECT geneid FROM geneidSummary");
		while(rs.next()) {
			String geneid = rs.getString(1);
			query = "SELECT COUNT(DISTINCT tag) FROM geneXML WHERE geneid = '" + geneid + "'";
			rs2 = stmt2.executeQuery(query);
			rs2.next();
			int numXML = rs2.getInt(1);
			rs2.close(); rs2 = null;

			stmt2.executeUpdate("UPDATE geneidSummary SET numXML = " + numXML + " WHERE geneid = '" + geneid + "' ");

			int nspecsTot  = getNumSpecs_GC(geneid, this.combinedFile, conn, 0.0);
			int nspecsUniq = getNumSpecs_GC(geneid, this.combinedFile, conn, this.wtTH);

			query = "UPDATE geneidSummary "
				  + "  SET numSpecsTot =  " + nspecsTot + ", "
				  + "      numSpecsUniq = " + nspecsUniq + " "
				  + "WHERE geneid = '" + geneid + "' ";
			stmt2.executeUpdate(query);

			int npepsTot  = getNumPeps_GC(geneid, this.combinedFile, conn, 0.0);
			int npepsUniq = getNumPeps_GC(geneid, this.combinedFile, conn, this.wtTH);

			query = "UPDATE geneidSummary "
				  + "  SET numPepsTot = " + npepsTot + ", "
				  + "      numPepsUniq = " + npepsUniq + " "
				  + "WHERE geneid = '" + geneid + "' ";
			stmt2.executeUpdate(query);

			if(console == null) {
				globals.cursorStatus(iter, msg);
				iter++;
			}
		}

		if(console != null) console.append("\n");
		else System.err.print("\n");

		stmt.close(); stmt = null;
		stmt2.close(); stmt2 = null;
	}



	/*************
	 *
	 * Function returns the number of peptides assigned to a given geneid
	 *
	 */
	public int getNumPeps_GC(String geneid, String sft, Connection conn, double wt) throws SQLException {
		Statement stmt  = conn.createStatement();

		ResultSet rs = null;
		String query = null;
		int retVal = 0;

		String tag = null;

		if(sft.equals(globals.combinedFile))  tag = "COMBINED";
		else tag = sft;

		query = "SELECT COUNT(*) "
			  + "FROM ("
			  + "  SELECT DISTINCT modPeptide "//, charge "
			  + "  FROM g2pep_ "
			  + "  WHERE tag = '" + tag + "' "
			  + "  AND geneid = '" + geneid + "' "
			  + "  AND wt >= " + wt + " "
			  + ")";
		rs = stmt.executeQuery(query);
		rs.next();
		retVal = rs.getInt(1);

		return retVal;
	}



	/***************
	 *
	 * Function returns the number of spectra assigned to a given geneid
	 *
	 */
	public int getNumSpecs_GC(String geneid, String sft, Connection conn, double wt) throws SQLException {

		Statement stmt  = conn.createStatement();

		ResultSet rs = null;
		String query = null;
		int retVal = 0;

		String tag = null;

		if(sft.equals(globals.combinedFile))  tag = "COMBINED";
		else tag = sft;


		query = "SELECT SUM(nspec) "
			  + "FROM g2pep_ "
			  + "WHERE tag = '" + tag + "' "
			  + "AND geneid = '" + geneid + "' "
			  + "AND wt >= " + wt + " ";
		rs = stmt.executeQuery(query);
		rs.next();
		retVal = rs.getInt(1);

		stmt.close();
		stmt = null;

		return retVal;
	}




	/****************
	 *
	 * Creates a results table that is gene-centric
	 *
	 */
	public void makeGeneResults(Connection conn, abacus_textArea console) throws SQLException {

		if(console != null) console.append("Creating gene-centric results table\n");
		else System.err.print("Creating gene-centric results table\n");

		Statement stmt = conn.createStatement();
		ResultSet rs  = null;
		String query = null;


		query = "CREATE CACHED TABLE geneResults ("
			  + "  geneid, isFwd, numXML, numGroups, maxPw, max_localPw, "
			  + "  maxIniProb, ALL_numSpecsTot, ALL_numSpecsUniq, "
			  + "  ALL_numPepsTot, ALL_numPepsUniq "
			  + ") AS ( "
			  + "SELECT geneid, isFwd, numXML, numGroups, maxPw, max_localPw, "
			  + "  maxIniProb, numSpecsTot, numSpecsUniq, numPepsTot, numPepsUniq "
			  + "FROM geneidSummary "
			  + "WHERE maxIniProb > " + this.maxIniProbTH + " "
			  + "AND numXML > 0 " 
			  + "GROUP BY geneid, isFwd, numXML, numGroups, maxPw, max_localPw, "
			  + "  maxIniProb, numSpecsTot, numSpecsUniq, numPepsTot , numPepsUniq  "
			  + "ORDER BY geneid ASC "
			  + ")WITH DATA ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("CREATE INDEX gr_idx1 ON geneResults(geneid)");
		stmt.executeUpdate("ALTER TABLE geneResults ADD COLUMN numProts INT DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE geneResults ADD COLUMN avgProtLen INT DEFAULT 0");

		query = "UPDATE geneResults "
			  + "  SET numProts = ?, "
			  + "      avgProtLen = ? "
			  + "WHERE geneid = ? ";
		PreparedStatement prep = conn.prepareStatement(query);

		rs = stmt.executeQuery("SELECT geneid FROM geneResults");
		while(rs.next()) {
			String geneid = rs.getString(1);

			int numProts = get_numProts(conn, geneid);
			int avgProtLen = 0;
			
			if( globals.fastaFile == null || globals.fastaFile.isEmpty() ) avgProtLen = 0;
			else {
				avgProtLen = get_avgProtLen(conn, geneid);
			}

			prep.setInt(1, numProts);
			prep.setInt(2, avgProtLen);
			prep.setString(3, geneid);
			prep.addBatch();
		}
		conn.setAutoCommit(false);
		prep.executeBatch();
		conn.setAutoCommit(true);

		rs.close(); rs = null;

		if(console != null) console.append("\n");
		else System.err.print("\n");

	}


	/*
	 * For gene-centric output this returns the number of proteins a given
	 * geneid has associated with it
	 */
	private int get_numProts(Connection conn, String geneid) throws SQLException {
		Statement stmt = conn.createStatement();
		ResultSet rs = null;
		String query = null;
		int retVal = 0;

		/*
		 * Determine if this geneid is a forward or decoy case
		 */
		if(geneid.startsWith("decoy-")) {

			int groupid = Integer.parseInt( geneid.substring(6) );

			query = "SELECT COUNT(DISTINCT protid) "
				  + "FROM combined "
				  + "WHERE groupid = " + groupid + " ";
			rs = stmt.executeQuery(query);
			rs.next();

			retVal = rs.getInt(1);
			rs.close(); rs = null;
		}
		else { // forward case
			query = "SELECT COUNT(DISTINCT protid) "
				  + "FROM gene2prot "
				  + "WHERE geneid = '" + geneid + "' ";
			rs = stmt.executeQuery(query);
			rs.next();

			retVal = rs.getInt(1);
			rs.close(); rs = null;
		}

		return retVal;
	}


	/*
	 * For gene-centric output this returns the average length of all the
	 * proteins derived from the given geneid
	 */
	private int get_avgProtLen(Connection conn, String geneid) throws SQLException {
		int retVal = 0;
		Statement stmt = conn.createStatement();
		ResultSet rs = null;
		String query = null;

		int sumLen = 0; // sum of all protein lengths
		int N = 0; // number of proteins
		float avg = 0; // average protein length

		/*
		 * Determine if the geneid is a forward or decoy case
		 */
		if(geneid.startsWith("decoy-")) { //decoy

			int groupid = Integer.parseInt( geneid.substring(6) );
			query = "SELECT DISTINCT protid "
				  + "FROM combined "
				  + "WHERE groupid = " + groupid + " ";
			rs = stmt.executeQuery(query);

			while(rs.next()) {
				String protid = rs.getString(1);

				if( globals.protLen.containsKey(protid) ) {
					N++;
					sumLen += globals.protLen.get(protid);
				}
			}
			avg = (float)sumLen / (float) N;
			retVal = Math.round(avg);

			rs.close(); rs = null;
		}
		else { // regular geneid

			query = "SELECT DISTINCT protid "
				  + "FROM gene2prot "
				  + "WHERE geneid = '" + geneid + "' ";
			rs = stmt.executeQuery(query);

			while(rs.next()) {
				String protid = rs.getString(1);

				if( globals.protLen.containsKey(protid) ) {
					N++;
					sumLen += globals.protLen.get(protid);
				}
			}
			avg = (float)sumLen / (float) N;
			retVal = Math.round(avg);

			rs.close(); rs = null;
		}

		return retVal;
	}


	/************
	 *
	 * Function creates a peptide usage table that is gene-centric
	 *
	 */
	public void makeGenePepUsageTable(Connection conn, abacus_textArea console) throws Exception {

		String msg = "Creating gene-centric peptide usage table (this could take a while)...";
		if(console != null) console.append(msg + "\n");

		Statement stmt = conn.createStatement();
		Statement stmt2 = conn.createStatement();

		ResultSet rs = null;
		ResultSet rs2 = null;
		PreparedStatement prep = null;
		String query = null;

		query = "CREATE CACHED TABLE genePepUsage_ ( "
			  + "  tag VARCHAR(250), "
			  + "  geneid VARCHAR(100), "
			  + "  modPeptide VARCHAR(250), "
//			  + "  charge INT, "
			  + "  nspecs INT, "
			  + "  wt DECIMAL(8,6) "
			  + ")";
		stmt.executeUpdate(query);

		query = "INSERT INTO genePepUsage_ ("
			  + "  SELECT gr.tag, r.geneid, gr.modPeptide, " //, gr.charge,"
			  + "    COUNT(DISTINCT px.specId), gr.wt "
			  + "  FROM geneXML gr, geneResults r, pepXML px "
			  + "  WHERE gr.tag = px.tag  "
			  + "  AND gr.geneid = r.geneid "
			  + "  AND gr.modPeptide = px.modPeptide "
//			  + "  AND gr.charge = px.charge "
			  + "  GROUP BY gr.tag, r.geneid, gr.modPeptide, gr.wt " //gr.charge, gr.wt "
			  + "  ORDER BY gr.tag, r.geneid "
			  + ") ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("ALTER TABLE genePepUsage_ ADD COLUMN numer INT DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE genePepUsage_ ADD COLUMN denom INT DEFAULT 0 ");
		stmt.executeUpdate("ALTER TABLE genePepUsage_ ADD COLUMN alpha DECIMAL(8,6) DEFAULT 0");
		stmt.executeUpdate("ALTER TABLE genePepUsage_ ADD COLUMN adjSpecs INT DEFAULT 0");

		stmt.executeUpdate("CREATE INDEX gpu_idx1 ON genePepUsage_(geneid)");
		stmt.executeUpdate("CREATE INDEX gpu_idx2 ON genePepUsage_(tag)");
		stmt.executeUpdate("CREATE INDEX gpu_idx3 ON genePepUsage_(tag,geneid)");
//		stmt.executeUpdate("CREATE INDEX gpu_idx4 ON genePepUsage_(modPeptide,charge)");
		stmt.executeUpdate("CREATE INDEX gpu_idx5 ON genePepUsage_(modPeptide)");
		stmt.executeUpdate("CREATE INDEX gpu_idx6 ON genePepUsage_(tag,modPeptide)");
		stmt.executeUpdate("CREATE INDEX gpu_idx7 ON genePepUsage_(tag,geneid,modPeptide)");


		query = "CREATE CACHED TABLE gWts_ ( "
			  + " tag, geneid, nspecsUniq "
			  + ") AS ("
			  + "  SELECT tag, geneid, SUM(nspecs) "
			  + "  FROM genePepUsage_ "
			  + "  WHERE wt > " + this.wtTH + " "
			  + "  GROUP BY tag, geneid "
			  + "  ORDER BY tag, geneid "
			  + ") WITH DATA ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("CREATE INDEX gw_idx1 ON gWts_(tag, geneid)");


		query = "UPDATE genePepUsage_ "
			  + "  SET numer = ( "
			  + "   SELECT gWts_.nspecsUniq "
			  + "   FROM gWts_ "
			  + "   WHERE gWts_.tag = genePepUsage_.tag "
			  + "   AND gWts_.geneid = genePepUsage_.geneid "
			  + " ) ";
		stmt.executeUpdate(query);
		stmt.executeUpdate("UPDATE genePepUsage_ SET numer = 0 WHERE numer IS NULL");


		query = "UPDATE genePepUsage_ "
			  + "  SET denom = ? "
			  + "WHERE tag = ? "
			  + "AND modPeptide = ? ";
		prep = conn.prepareStatement(query);

		query = "SELECT tag, modPeptide "
			  + "FROM genePepUsage_ "
			  + "GROUP BY tag, modPeptide ";
		rs = stmt.executeQuery(query);

		int iter = 0;
		while(rs.next()) {
			String tag = rs.getString(1);
			String modPep = rs.getString(2);

			query = "SELECT DISTINCT geneid, numer "
				  + "FROM genePepUsage_ "
				  + "WHERE tag = '" + tag + "' "
				  + "AND modPeptide = '" + modPep + "' ";
			rs2 = stmt2.executeQuery(query);

			int nspec = 0;
			while(rs2.next()) {
				nspec += rs2.getInt(2);
			}

			prep.setInt(1, nspec);
			prep.setString(2, tag);
			prep.setString(3, modPep);
			prep.addBatch();

			if(console == null) globals.cursorStatus(iter, msg);
			iter++;
		}
		conn.setAutoCommit(false);
		prep.executeBatch();
		conn.setAutoCommit(true);

		rs.close(); rs = null;

		// this prevents division by zero
		stmt.executeUpdate("UPDATE genePepUsage_ SET denom = 1 WHERE denom is NULL");
		stmt.executeUpdate("UPDATE genePepUsage_ SET denom = 1 WHERE denom = 0");



		query = "UPDATE genePepUsage_ "
			  + "  SET alpha = (CAST(numer AS DECIMAL(16,6)) / CAST(denom AS DECIMAL(16,6)))";
		stmt.executeUpdate(query);

		query = "UPDATE genePepUsage_ "
			  + "  SET adjSpecs = ( CAST("
			  + "    ROUND( (CAST(nspecs AS DECIMAL(16,6)) * alpha), 0) "
			  + "    AS INT) "
			  + ") ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("DROP INDEX IF EXISTS gw_idx1");
		stmt.executeUpdate("DROP TABLE IF EXISTS gWts_");

		stmt.close(); stmt = null;
		stmt2.close(); stmt2 = null;
		prep.close(); prep = null;

		if(console != null) console.append("\n\n");
		else System.err.print("\n\n");
	}



	/*************
	 *
	 * Appends independent experiment data to gene-centric results table
	 *
	 */
	public void appendIndividualExpts_GC(Connection conn, abacus_textArea console) throws Exception {

		if(console != null) console.append("Appending individual experiment results\n");
		else System.err.print("Appending individual experiment results\n");

		Statement stmt  = conn.createStatement();
		Statement stmt2 = conn.createStatement();
		Statement stmt3 = conn.createStatement();
		Statement stmt4 = conn.createStatement();
		ResultSet rs  = null;
		ResultSet rs2 = null;
		String query = null;


		/*
		 * Iterate over protXML files
		 */
		rs = stmt.executeQuery("SELECT DISTINCT tag FROM srcFileTags WHERE fileType = 'prot' ORDER BY tag ASC");

		while(rs.next()) {
			String tag = rs.getString(1);
			String msg = "  Adding columns for " + tag;
			if(console != null) console.append(msg);

			int iter = 0;

			/*
			 * Add necessary columns to geneResults
			 */
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_maxPw DECIMAL(8,6) DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_max_localPw DECIMAL(8,6) DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_maxIniProb DECIMAL(8,6) DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_numSpecsTot INT DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_numSpecsUniq INT DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_numSpecsAdj INT DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_numPepsTot INT DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE geneResults ADD COLUMN " + tag + "_numPepsUniq INT DEFAULT 0");


			query = "SELECT geneid, maxPw, max_localPw, MAX(iniProb) "
				  + "FROM geneXML "
				  + "WHERE tag = '" + tag + "' "
				  + "GROUP BY geneid, maxPw, max_localPw ";
			rs2 = stmt3.executeQuery(query);

			while(rs2.next()) {
				String geneid = rs2.getString(1);

				query = "UPDATE geneResults "
					  + "  SET (" + tag + "_maxPw, "
					  + "       " + tag + "_max_localPw, "
					  + "       " + tag + "_maxIniProb "
					  + "      ) = ( "
					  + "  " + rs2.getDouble(2) // maxPw
					  + " ," + rs2.getDouble(3) // maxLocalPw
					  + " ," + rs2.getDouble(4) // maxIniProb
					  + ") "
					  + "WHERE geneid = '" + geneid + "' ";
				stmt4.executeUpdate(query);

				int nsT = getNumSpecs_GC(geneid, tag, conn, 0); // total num specs
				int nsU = getNumSpecs_GC(geneid, tag, conn, this.wtTH); // unique num specs

				int npT = getNumPeps_GC(geneid, tag, conn, 0); // total num peps
				int npU = getNumPeps_GC(geneid, tag, conn, this.wtTH); // unique num peps


				query = "UPDATE geneResults "
					  + "  SET (" + tag + "_numSpecsTot, "
					  + "       " + tag + "_numSpecsUniq,"
					  + "       " + tag + "_numPepsTot,"
					  + "       " + tag + "_numPepsUniq "
					  + "      ) = ( "
					  + "    " + nsT + ", " + nsU + ", " + npT + ", " + npU + " "
					  + ") "
					  + "WHERE geneid = '" + geneid + "' ";
				stmt4.executeUpdate(query);

				if(console == null) {
					globals.cursorStatus(iter, msg);
					iter++;
				}
			}

			if(console != null) console.append("\n");
			else System.err.print("\n");
		}


		/*
		 * Now fill in the adjusted spectral count fields in geneResults
		 */
		query = "SELECT tag, geneid, SUM(adjSpecs) "
			  + "FROM genePepUsage_ "
			  + "GROUP BY tag, geneid "
			  + "ORDER BY tag, geneid ";
		rs = stmt.executeQuery(query);

		while(rs.next()) {
			query = "UPDATE geneResults "
				  + "  SET " + rs.getString(1) + "_numSpecsAdj = " + rs.getInt(3)
				  + "WHERE geneid = '" + rs.getString(2) + "' ";
			stmt2.executeUpdate(query);
		}
		rs.close(); rs = null;

		stmt.close(); stmt = null;
	}


	/***************
	 *
	 * Function creates a temporary table called by 'appendIndividualExpts'
	 * This same table is dropped within that function.
	 *
	 * @param conn
	 * @throws SQLException
	 */
	public void makeTempGene2pepTable(Connection conn) throws SQLException {
		Statement stmt = conn.createStatement();
		String query = null;

		/*
		 * Create temporary table for peptide and gene counting
		 */
		query = "CREATE CACHED TABLE g2pep_ ("
			  + "  tag VARCHAR(250), "
			  + "  geneid VARCHAR(100), "
			  + "  modPeptide VARCHAR(250), "
//			  + "  charge INT, "
			  + "  wt DECIMAL(8,6), "
			  + "  nspec INT DEFAULT 0 "
			  + ")";
		stmt.executeUpdate(query);

		// geneCombined peptides
		query = "INSERT INTO g2pep_ ( "
			  + "SELECT 'COMBINED', c.geneid, c.modPeptide, "
//			  + "  c.charge, c.wt, COUNT(DISTINCT px.specId) "
			  + "  c.wt, COUNT(DISTINCT px.specId) "
			  + "FROM geneCombined c, pepXML px "
			  + "WHERE c.modPeptide = px.modPeptide "
//			  + "AND c.charge = px.charge "
			  + "AND px.iniProb >= " + globals.iniProbTH + " "
			  + "GROUP BY c.geneid, c.modPeptide, c.wt " //c.charge, c.wt "
			  + "ORDER BY c.geneid, c.modPeptide "
			  + ")";
		stmt.executeUpdate(query);

		// geneXML peptides
		query = "INSERT INTO g2pep_ ( "
			  + "SELECT gx.tag, gx.geneid, gx.modPeptide, "
//			  + "  gx.charge, gx.wt, COUNT(DISTINCT px.specId) "
			  + "  gx.wt, COUNT(DISTINCT px.specId) "
			  + "FROM geneXML gx, pepXML px "
			  + "WHERE gx.tag = px.tag "
			  + "AND gx.modPeptide = px.modPeptide "
//			  + "AND gx.charge = px.charge "
			  + "AND px.iniProb >= " + globals.iniProbTH + " "
			  + "GROUP BY gx.tag, gx.geneid, gx.modPeptide, gx.wt " //gx.charge, gx.wt "
			  + "ORDER BY gx.tag, gx.geneid, gx.modPeptide "
			  + ") ";
		stmt.executeUpdate(query);

		stmt.executeUpdate("CREATE INDEX g2pep_idx1 ON g2pep_(tag)");
		stmt.executeUpdate("CREATE INDEX g2pep_idx2 ON g2pep_(geneid)");
//		stmt.executeUpdate("CREATE INDEX g2pep_idx3 ON g2pep_(modPeptide,charge)");
		stmt.executeUpdate("CREATE INDEX g2pep_idx4 ON g2pep_(tag,geneid)");
		stmt.executeUpdate("CREATE INDEX g2pep_idx5 ON g2pep_(tag,modPeptide)");
		stmt.executeUpdate("CREATE INDEX g2pep_idx6 ON g2pep_(tag,geneid,modPeptide)");
//		stmt.executeUpdate("CREATE INDEX g2pep_idx7 ON g2pep_(tag,geneid,modPeptide,charge)");

		stmt.close();
		stmt = null;
	}



	/**************
	 *
	 * Function appends gene descriptions to the geneResults table.
	 * This function is only called if the gene2prot mapping file contained a
	 * third column with text that describes each gene.
	 *
	 * @param conn
	 * @throws SQLException
	 *
	 */
	public void appendGeneDescriptions(Connection conn) throws SQLException {
		Statement stmt = conn.createStatement();
		String query = null;

		query = "ALTER TABLE geneResults ADD COLUMN geneDescription VARCHAR(1000) ";
		stmt.executeUpdate(query);

		query = "UPDATE geneResults "
			  + "  SET geneDescription = ? "
			  + "WHERE geneid = ? ";
		PreparedStatement prep = conn.prepareStatement(query);

		query = "SELECT g2p.geneid, g2p.geneDefline "
			  + "FROM gene2prot g2p, geneResults r "
			  + "WHERE g2p.geneid = r.geneid "
			  + "AND r.isFwd = 1 "
			  + "GROUP BY g2p.geneid, g2p.geneDefline ";
		ResultSet rs = stmt.executeQuery(query);

		while(rs.next()) {
			String geneid = rs.getString(1);
			String defline = rs.getString(2);

			prep.setString(1, defline);
			prep.setString(2, geneid);
			prep.addBatch();
		}
		conn.setAutoCommit(false);
		prep.executeBatch();
		conn.setAutoCommit(true);
		rs.close(); rs = null;

		query = "UPDATE geneResults "
			  + " SET geneDescription = 'DECOY MATCH' "
			  + "WHERE isFwd = 0 ";
		stmt.executeUpdate(query);

		stmt.close();
		stmt = null;
	}

	public void getNSAF_values_gene(Connection conn, abacus_textArea console) throws SQLException {
		String query = null;
		String msg = null;
		Statement stmt = conn.createStatement();
		Statement stmt2 = conn.createStatement();
		Statement stmt3 = conn.createStatement();

		ResultSet rs = null;
		ResultSet rs2 = null;

		msg = "\nCreating NSAF values table (gene-centric)\n";
		if(console != null) console.append(msg);
		else System.err.print(msg);

		stmt.executeUpdate("DROP TABLE IF EXISTS nsaf_p1");
		stmt.executeUpdate("DROP TABLE IF EXISTS nsaf");

		// create initial table
		query = "CREATE CACHED TABLE nsaf_p1 ("
			  + "  geneid "
			  + ") AS ( "
			  + "SELECT geneid "
			  + "FROM geneResults "
			  + "WHERE isFwd = 1 "
			  + "GROUP BY geneid "
			  + "ORDER BY geneid ASC "
			  + ")WITH DATA ";
		stmt.executeUpdate(query);
		stmt.executeUpdate("CREATE INDEX nsaf_p1_idx1 ON nsaf_p1(geneid)");

		// create final NSAF table
		query = "CREATE CACHED TABLE nsaf ("
			  + "  geneid "
			  + ") AS ( "
			  + "SELECT geneid "
			  + "FROM geneResults "
			  + "WHERE isFwd = 1 "
			  + "GROUP BY geneid "
			  + "ORDER BY geneid ASC "
			  + ")WITH DATA ";
		stmt.executeUpdate(query);
		stmt.executeUpdate("CREATE INDEX nsaf_idx1 ON nsaf(geneid)");

		/*
		 * All NSAF values are multiplied by a factor to raise their values.
		 * This is done because of a rounding error problem inherent to HSQLDB.
		 * I can't find the cause of the bug, but by multiplying the values by this
		 * factor I can eliminate it.
		 * Our factor is computed as 10^N where N = # of digits in representing
		 * number of proteins reported + 1
		 * example: if 377 proteins are reported then N = 3+1, 377 consists of 3 digits
		 */
		rs = stmt.executeQuery("SELECT COUNT(geneid) FROM geneResults WHERE isFwd = 1");
		rs.next();
		String numProts = Integer.toString( rs.getInt(1) );
		int factor = numProts.length() + 1;
		double NSAF_FACTOR = Math.pow(10,factor); // we multiply all NSAF values by this
		globals.NSAF_FACTOR = NSAF_FACTOR;
		rs.close();

		msg = "  NSAF_FACTOR = 10^" + factor + " = " + NSAF_FACTOR + "\n";
		if(console != null) console.append(msg);
		else System.err.print(msg);



		// now you need the add columns for spectral counts
		rs = stmt.executeQuery("SELECT DISTINCT tag FROM srcFileTags WHERE fileType = 'prot' ORDER BY tag ASC");

		while(rs.next()) {
			String tag = rs.getString(1);

			stmt2.executeUpdate("ALTER TABLE nsaf_p1 ADD COLUMN " + tag + "_specsTot DOUBLE DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE nsaf_p1 ADD COLUMN " + tag + "_specsUniq DOUBLE DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE nsaf_p1 ADD COLUMN " + tag + "_specsAdj DOUBLE DEFAULT 0");

			stmt2.executeUpdate("ALTER TABLE nsaf ADD COLUMN " + tag + "_totNSAF DOUBLE DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE nsaf ADD COLUMN " + tag + "_uniqNSAF DOUBLE DEFAULT 0");
			stmt2.executeUpdate("ALTER TABLE nsaf ADD COLUMN " + tag + "_adjNSAF DOUBLE DEFAULT 0");

			// get spectral counts divided by avg. protein length
			query = "SELECT geneid, avgProtLen, "
				  + "  " + tag + "_numSpecsTot, "
				  + "  " + tag + "_numSpecsUniq, "
				  + "  " + tag + "_numSpecsAdj "
				  + "FROM geneResults "
				  + "ORDER BY geneid ";

			rs2 = stmt2.executeQuery(query);

			// assign spectral-count/length to nsaf_p1 table
			while(rs2.next()) {
				String geneid = rs2.getString(1);
				double protLen = rs2.getDouble(2);

				double tot = rs2.getDouble(3) / protLen;
				double uniq = rs2.getDouble(4) / protLen;
				double adj = rs2.getDouble(5) / protLen;

				query = "UPDATE nsaf_p1 "
					  + "  SET " + tag + "_specsTot = " + tot + ", "
					  + "      " + tag + "_specsUniq = " + uniq + ", "
					  + "      " + tag + "_specsAdj = " + adj + " "
					  + "WHERE geneid = '" + geneid + "' ";

				stmt3.executeUpdate(query);
			}
			rs2.close();

			// now add up columns of current tag
			query = "SELECT SUM(" + tag + "_specsTot), "
				  + "       SUM(" + tag + "_specsUniq), "
				  + "       SUM(" + tag + "_specsAdj) "
				  + "FROM nsaf_p1 ";
			rs2 = stmt2.executeQuery(query);
			rs2.next();

			double totSum  = rs2.getDouble(1);
			double uniqSum = rs2.getDouble(2);
			double adjSum  = rs2.getDouble(3);

			rs2.close();


			query = "SELECT geneid, "
				  + "  " + tag + "_specsTot, "
				  + "  " + tag + "_specsUniq, "
				  + "  " + tag + "_specsAdj "
				  + "FROM nsaf_p1 "
				  + "GROUP BY geneid, "
				  + "  " + tag + "_specsTot, "
				  + "  " + tag + "_specsUniq, "
				  + "  " + tag + "_specsAdj "
				  + "ORDER BY geneid ASC ";

			rs2 = stmt2.executeQuery(query);



			while(rs2.next()) {
				String geneid = rs2.getString(1);

				double x_t = rs2.getDouble(2);
				double x_u = rs2.getDouble(3);
				double x_a = rs2.getDouble(4);

				double nsaf_t = (x_t/totSum) * NSAF_FACTOR;
				double nsaf_u = (x_u/uniqSum) * NSAF_FACTOR;
				double nsaf_a = (x_a/adjSum) * NSAF_FACTOR;


				query = "UPDATE nsaf "
					  + "  SET " + tag + "_totNSAF =  " + nsaf_t + ","
					  + "      " + tag + "_uniqNSAF = " + nsaf_u + ","
					  + "      " + tag + "_adjNSAF =  " + nsaf_a
					  + "WHERE geneid = '" + geneid + "' ";

				stmt3.executeUpdate(query);
			}
			rs2.close();
		}
		rs.close();

		stmt.executeUpdate("DROP INDEX nsaf_p1_idx1");
		stmt.executeUpdate("DROP TABLE nsaf_p1");

		stmt.close();
		stmt2.close();
		stmt3.close();

		/*
		 * Now we need to add NSAF values to results table
		 */
		reformat_results(conn, console);
	}



} // end hyperSQLObject_gene Class

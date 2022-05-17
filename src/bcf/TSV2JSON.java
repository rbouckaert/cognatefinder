package bcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import beast.app.beauti.BeautiDoc;
import beast.app.util.Application;
import beast.core.Input;
import beast.core.util.Log;

/**
 * Creates JSON file from word list
 * with aligned phoneme sequences
 *
 */
public class TSV2JSON extends TSV2Nexus {

	final public Input<File> mappingInput = new Input<>("map", "phoneme mapping used to pre-process phoneme strings to get rid of infrequently occurring phonemes. "
			+ "This is a tab delimited file with first column source phoneme and second column the target phoneme. "
			+ "All occurrances of source phonemes will be replaced by target phonemes. "
			+ "Ignored when not specified.");

	final public Input<File> encodingInput = new Input<>("encoding", "phoneme mapping used to post-process phoneme strings to get rid of infrequently non-ascii characters. "
			+ "Like the mapping file, this is tab delimited with two columns");
	
	
	final public Input<File> conceptInput = new Input<>("concept", "a list of concepts to remove (each meaning class on a different line in the file)");

	final public Input<Boolean> suppresWordGapInput = new Input<>("suppresWordGap", "removes word gaps (encoded as '+' or '_') from segments", true);
	final public Input<Boolean> includeLoanInput = new Input<>("loan", "include words marked as loan=true?", true);
	final public Input<Boolean> removeConstantCognatesInput = new Input<>("variableOnly", "filter out constant cognates", false);
	
	final public Input<Integer> minTaxaPerCognateInput = new Input<>("min", "Minimum number of taxa per cognate for the cognate tree to be included", 2);
	
	final static String VOWEL     = "aeoiuyɛəøæœɨɔ";
	final static String CONSONANT = "bdfghjklmnprsttvwŋɢʃʔβ";

	
	Map<String, String> phonemeMapping = new HashMap<>();
	Map<String, String> encodingMapping = new HashMap<>();
	List<String> conceptRemoval = new ArrayList<>();
	
	@Override
	public void initAndValidate() {
		matrix = new DolgoScore();
		matrix = new SCAScore();
	}

	
	String [] token;
	String [] conceptList;
	int [] cogid;
	Set<String> doculects;
	String [] doculect;
	String [] docIds;
	String [] sequence;
	String [] concept;
	boolean [] constantZeroColumn;
	int maxCognateNumber;
	int cognateCount;
	char [][] cognateCharSeqs;
	int [] startPhonemesMC;
	int [] startPhonemesCognate;
	int [] startCognates;
	List<Integer> cogidUnique;
	
	
	List<String> cognateNumberRange;
	List<List<String>> taxaPerCognate;
	
	
	@Override
	public void run() throws Exception {
		if (tsvInput.get() == null || tsvInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("A valid TSV file must be specified");
		}
		if (mappingInput.get() != null) {
			processMapping(mappingInput.get(), phonemeMapping);
		}
		if (encodingInput.get() != null) {
			processMapping(encodingInput.get(), encodingMapping);
		}
		
		TSVImporter importer = new TSVImporter(tsvInput.get(), languagesInput.get());
		
		
		// Concept list
		conceptList = importer.getColumn("CONCEPT");
		if (conceptList == null) {
			conceptList = importer.getColumn("PARAMETER_ID");
		}
		if (conceptList == null) throw new IllegalArgumentException("Please ensure there is a column named 'CONCEPT' or 'PARAMETER_ID' which contains meaning classes");
		
		
		// Remove some rows?
		boolean[] remove = new boolean[conceptList.length];
		for (int i = 0; i < conceptList.length; i ++) {
			remove[i] = false;
		}
		
		
		// Remove any concepts from the data?
		if (conceptInput.get() != null) {
			processList(conceptInput.get(), conceptRemoval);
			
			
			for (int i = 0; i < conceptList.length; i ++) {
				String concept = conceptList[i];
				boolean isBadWord = conceptRemoval.contains(concept);
				
				if (isBadWord) {
					Log.warning("Removing " + concept + " on row " + (i+1));
					remove[i] = true;
				}
			}
			
		
		}
		

		// Remove tokens marked as loan?
		if (!includeLoanInput.get()) {
			
			// Word tokens
			token = importer.getColumn("TOKENS");
			if (token == null ) {
				token = importer.getColumn("SEGMENTS");
			}
			
			
			
			String[] loanStr = importer.getColumn("LOAN");
			if (loanStr == null) {
				loanStr = importer.getColumn("Loan");
				if (loanStr == null) {
					loanStr = importer.getColumn("loan");
				}
			}
			
			if (loanStr == null) throw new IllegalArgumentException("Cannot find a column named LOAN / Loan / loan");
			if (loanStr.length != token.length) throw new IllegalArgumentException("Inconsistent number of rows " + loanStr.length + "!="  + token.length);
			
			for (int rowNum = 0; rowNum < loanStr.length; rowNum++) {
				String val = loanStr[rowNum];
				boolean loaned = val != null && !val.isEmpty() && val.toLowerCase().trim().equals("true");
				if (loaned) {
					Log.warning("Removing '" + token[rowNum] + "' (row " + (rowNum+1) + ") because it was borrowed");
					remove[rowNum] = true;
				}
			}
			
		}
		
		
		// Remove rows from the table
		importer.removeRows(remove);
		conceptList = importer.getColumn("CONCEPT");
		if (conceptList == null) {
			conceptList = importer.getColumn("PARAMETER_ID");
		}
		
		// Word tokens
		token = importer.getColumn("TOKENS");
		if (token == null ) {
			token = importer.getColumn("SEGMENTS");
		}
		
		
		// Standardise tokens (which have not been removed)
		standardiseTokens(token);
		
		
		// Doculects
		doculect = importer.getColumn("DOCULECT");
		if (doculect == null) {
			doculect = importer.getColumn("LANGUAGE_ID");
		}
		doculects = new LinkedHashSet<>();
		for (String s : doculect) {
			doculects.add(s);
		}
		

		
		
		cogid = importer.getColumnAsInt("COGID");
		if (cogid == null ) {
			String [] cognacy = importer.getColumn("COGNACY");
			cogid = new int[conceptList.length];
			Map<String,Integer> conceptmap = new HashMap<>();
			for (int i = 0; i < cogid.length; i++) {
				String c = cognacy[i];
				if (c != null && c.trim().length() > 0) {
					if (!conceptmap.containsKey(c)) {
						conceptmap.put(c, conceptmap.size() + 1);
					}
					cogid[i] = conceptmap.get(c);
				} else {
					cogid[i] = -1;
				}
			}
		}
		
		/*
		for (int i = 0; i < cogid.length; i ++) {
			System.out.print(cogid[i] + " ");
		}
		System.out.println();
		*/
		
		Set<String> concepts;
		concepts = new LinkedHashSet<>();
		for (String s : conceptList) {
			if (s != null) {
				concepts.add(s);
			}
		}
		concept = concepts.toArray(new String[]{});
		Arrays.sort(concept);
		
		
		// Count cognates
		maxCognateNumber = 0;
		for (int d : cogid) {
			maxCognateNumber = Math.max(d, maxCognateNumber);
		}
		maxCognateNumber++;
		
		cogidUnique = new ArrayList<>();
		cogidUnique.add(-1); // First site is dummy
		//startCognates = new int[]
		for (int s : cogid) {
			if (s > 0 && !cogidUnique.contains(s)) {
				cogidUnique.add(s);
			}
		}
		cognateCount = cogidUnique.size();

		
		boolean [][] docHasConcept = new boolean[doculects.size()][concept.length];
		
		// group by cogid
		List<C> [] seqs = new List[maxCognateNumber];
		for (int i = 0; i < cogid.length; i++) {
			int id = cogid[i];
			if (id > 0) {
				if (seqs[id] == null) {
					seqs[id] = new ArrayList<>();
				}
				C c = new C(token[i], doculect[i]);
				if (token[i] != null) {
					seqs[id].add(c);
				}
			}
		}
		
		// align individual cogids
		for (int i = 0; i < maxCognateNumber; i++) {
			alignCognate(seqs[i], extendGapPenaltyInput.get(), openGapPenaltyInput.get());
			//alignCognate(seqs[i], 4.5f, -0.9f);
		}
		
		// convert cognate alignments to sequences
		StringBuilder [] charSeqs = new StringBuilder[doculects.size()];
		for (int i = 0; i < charSeqs.length; i++) {
			charSeqs[i] = new StringBuilder();
		}
		docIds = doculects.toArray(new String[]{});
		startPhonemesMC = new int[concepts.size()];
		startCognates = new int[concepts.size()];
	
		
		
		boolean [] done = new boolean[maxCognateNumber];
		
		
		// Cognate start/finish positions
		startPhonemesCognate = new int[cognateCount];
		for (int cognateNum = 1; cognateNum < cognateCount; cognateNum++) {
			int i = cogidUnique.get(cognateNum);
			int len = 0; 
			if (seqs[i].size() > 0) {
				len = seqs[i].get(0).aligned.characters.length - 2;
			}
			startPhonemesCognate[cognateNum] = len + (cognateNum > 0 ?startPhonemesCognate[cognateNum-1] : 0);
		}
		
		
		
		// Meaning class start/finish positions
		int j = 0;
		int MAX = 256;
		int [] cognatePerConceptHistogram = new int[MAX];
		List<String> [] conceptSet = new List[MAX];
		for (int i = 0; i < conceptSet.length; i++) {
			conceptSet[i] = new ArrayList<>();
		}
		int [] dataPerConceptHistogram = new int[MAX];
		List<String> [] dataSet = new List[MAX];
		for (int i = 0; i < dataSet.length; i++) {
			dataSet[i] = new ArrayList<>();
		}
		for (String c : concept) {
			int lenPhoneme = 0;
			int lenCognate = 0;
			int cognatePerConcept = 0;
			int dataPerConcept = 0;
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && conceptList[k].equals(c)) {
					int i = cogid[k];
					if (cogid[k] >= 0 && !done[i]) {
						lenPhoneme += toCharSeqs(seqs[i], docIds, charSeqs, j, docHasConcept);
						lenCognate ++;
						done[i] = true;
						cognatePerConcept++;
					}
					dataPerConcept++;
				}
			}
			cognatePerConceptHistogram[cognatePerConcept]++;
			conceptSet[cognatePerConcept].add(c);
			dataPerConceptHistogram[dataPerConcept]++;
			dataSet[dataPerConcept].add(c);
			startPhonemesMC[j] = lenPhoneme + (j>0?startPhonemesMC[j-1] : 0);
			startCognates[j] = lenCognate + (j>0?startCognates[j-1] : 0); 
			j++;
		}
				
		
		
		sequence = new String[docIds.length];
		for (int i = 0; i < docIds.length; i++) {
			sequence[i] = cleanUp(charSeqs[i].toString());
		}

		separateVowelsAndConsonants(sequence);

		
		// create binary cognate alignment
		cognateCharSeqs = new char[doculects.size()][cognateCount];
		for (int i = 0; i < cognateCharSeqs.length; i++) {
			Arrays.fill(cognateCharSeqs[i],'0');
		}
		for (int i = 0; i < cogid.length; i++) {
			int cogid_ = cogid[i];
			if (cogid_ >= 0) {
				String doculect_ = doculect[i];
				if (doculect_ != null) {
					int k = indexOf(doculect_, docIds);
					int siteNum = cogidUnique.indexOf(cogid_);
					if (siteNum >= 0) {
						cognateCharSeqs[k][siteNum] = '1';
					}
				}
			}
		}
		// check for constant columns
		constantZeroColumn = new boolean[cognateCount];
		for (int k = 1; k < constantZeroColumn.length; k++) {
			boolean isConstant = true;
			
			char c = cognateCharSeqs[0][k];
			for (int i = 1; i < doculects.size(); i++) {
				if (cognateCharSeqs[i][k] != c) {
					isConstant = false;
					break;
				}
			}
			if (isConstant) {
				if (c == '0') {
					constantZeroColumn[k] = true;
				}
				if (removeConstantCognatesInput.get()) System.err.println("Constant column ("+cognateCharSeqs[0][k]+") in binary sequence at column " + (k+1));
			}
		}
		
		
		
		
		// For each cognate, which taxa have a 1?
		taxaPerCognate = new ArrayList<>();
		for (int cognateNum = 0; cognateNum < cognateCount; cognateNum++) {
			List<String> taxa = new ArrayList<>();
			taxaPerCognate.add(taxa);
		}
		for (int i = 0; i < cogid.length; i++) {
			int cogid_ = cogid[i];
			if (cogid_ >= 0) {
				String doculect_ = doculect[i];
				if (doculect_ != null) {
					
					// Cognate index 
					int cognateIndex = cogidUnique.indexOf(cogid_);
					if (!taxaPerCognate.get(cognateIndex).contains(doculect_)) {
						taxaPerCognate.get(cognateIndex).add(doculect_);
					}
				}
			}
		}
		
		
		// Only include cognates with more than 1 taxon
		cognateNumberRange = new ArrayList<>();
		for (int cognateNum = 1; cognateNum < cognateCount; cognateNum++) {
			if (taxaPerCognate.get(cognateNum).size() >= minTaxaPerCognateInput.get()) {
				cognateNumberRange.add("" + (cognateNum+1));
			}
		}
		
		
		output();
		
		// process stats
		System.err.println();
		double mean = 0;
		double count = 0;
		for (int i = 0; i < cognatePerConceptHistogram.length; i++) {
			if (cognatePerConceptHistogram[i] > 0) {
				mean += cognatePerConceptHistogram[i] * i;
				count += cognatePerConceptHistogram[i];
				System.err.println(cognatePerConceptHistogram[i] + " concepts with " + i + " cognates " + 
					Arrays.toString(conceptSet[i].toArray()));
			}
		}
		System.err.println("On average: " + (mean / count) + " cognates per concept\n");

		mean = 0;
		count = 0;
		for (int i = 0; i < dataPerConceptHistogram.length; i++) {
			if (dataPerConceptHistogram[i] > 0) {
				mean += dataPerConceptHistogram[i] * i;
				count += dataPerConceptHistogram[i];
				System.err.println(dataPerConceptHistogram[i] + " concepts with " + i + " data " + 
					Arrays.toString(dataSet[i].toArray()));
			}
		}
		System.err.println("On average: " + (mean / count) + " data per concept\n");
		
		mean = 0;
		for (int i = 0; i < doculects.size(); i++) {
			int k = 0;
			for (boolean b : docHasConcept[i]) {
				if (b) {
					k++;
				}
			}
			mean += k;
			System.err.println(docIds[i] + " covered by " + k + " concepts");
		}
		System.err.println("On average: " + (mean / doculects.size()) + " (out of "
				+ concept.length + ") concepts covered per doculect\n");

		Log.warning("Done!");
	}
		
	private void output() throws FileNotFoundException {
		// output results
		StringBuilder buf = new StringBuilder();
		buf.append("{");
		buf.append("\"sequences\":\"\n");

		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getName());
			out = new PrintStream(outputInput.get());
		}

		
		Set<String> phonemes = new LinkedHashSet<>();
		Set<String> vowels = new LinkedHashSet<>();
		Set<String> consonants = new LinkedHashSet<>();

		for (int i = 0; i < docIds.length; i++) {
			if (docIds[i] != null) {
				buf.append("<sequence taxon='" + docIds[i] + "' ");
				if (docIds[i].length() < 25) {
					buf.append("                         ".substring(docIds[i].length()));
				}
				buf.append("value='");
				sequence[i] = processPhonemes(sequence[i], phonemes, vowels, consonants);
				
				buf.append(sequence[i]);
				buf.append("'/>\n");
			}
		}
		buf.append("\",\n");
		
		buf.append("\"binsequences\":\"\n");
		for (int i = 0; i < docIds.length; i++) {
			if (docIds[i] != null) {
				buf.append("<sequence taxon='" + docIds[i] + "' ");
				if (docIds[i].length() < 25) {
					buf.append("                         ".substring(docIds[i].length()));
				}
				buf.append("value='");
				for (int k = 0; k < cognateCount; k++) {
					if (!removeConstantCognatesInput.get() || !constantZeroColumn[k]) {
						buf.append(cognateCharSeqs[i][k]);
					}
				}
				buf.append("'/>\n");
			}
		}
		buf.append("\",\n");

		//buf.append("// _. = gap between words, .. = other cognate, -. = gap due to alignment\n");
		appendDataType(phonemes, "phonemes", buf);
		appendDataType(vowels, "vowels", buf);
		appendDataType(consonants, "consonants", buf);

		// remove ".." representing other cognates
		appendDataTypeWithMissing(phonemes, "phonemes", buf, "..", 1);
		appendDataTypeWithMissing(vowels, "vowels", buf, "..", 1);
		appendDataTypeWithMissing(consonants, "consonants", buf, "..", 1);
		
		// remove ".-" representing missing data due to alignment
		appendDataTypeWithMissing(phonemes, "phonemes", buf, "-.", 2);
		appendDataTypeWithMissing(vowels, "vowels", buf, "-.", 2);
		appendDataTypeWithMissing(consonants, "consonants", buf, "-.", 2);

		buf.append("\"words\":\"");
		for (int i = 0; i < concept.length; i++) {
			buf.append(concept[i].replaceAll("[ ,]", "_"));
			buf.append(concept[i].replaceAll("/", "-"));
			if (i < concept.length - 1) {
				buf.append(",");
			}
		}
		buf.append("\",\n");
		buf.append("\"words-1\":\"");
		for (int i = 1; i < concept.length; i++) {
			buf.append(concept[i].replaceAll("[ ,]", "_"));
			buf.append(concept[i].replaceAll("/", "-"));
			if (i < concept.length - 1) {
				buf.append(",");
			}
		}
		buf.append("\",\n");
		
		
		// Cognate numbers with more than 1 taxon
		buf.append("\"n-cognate-trees\":\"" + cognateNumberRange.size() + "\",\n");
		buf.append("\"cognate-tree-range\":\"");
		buf.append(String.join(",", cognateNumberRange));
		buf.append("\",\n");
				

		
		// Taxonset for multispecies coalescent
		buf.append("\"taxonset\":\"\n");
		for (int i = 0; i < docIds.length; i++) {
			if (docIds[i] == null || docIds[i].isEmpty()) continue;
			buf.append("<taxon id='" + docIds[i] + ".lang' spec='TaxonSet'>\n");
			buf.append("\t<taxon id='" + docIds[i] + "' spec='Taxon' />\n");
			buf.append("</taxon>\n");
		}
		buf.append("\",\n");
		
		
		// List of meaning class IDs
		buf.append("\"meaning-classes\":\"");
		for (int i = 0; i < concept.length; i++) {
			buf.append("MC_" + concept[i].replaceAll("[ ,]", "_").replaceAll("/", "_"));
			if (i < concept.length-1) buf.append(",");
		}		
		buf.append("\",\n");
		
		
		
		// List of meaning class IDs (but only ones which have a cognate tree)
		List<String> conceptsFiltered = new ArrayList<>();
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			int cogid_ = cogidUnique.get(i-1);
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && cogid[k] == cogid_) {
					String concept = conceptList[k];
					concept = "MC_" + concept.replaceAll("[ ,]", "_").replaceAll("/", "_");
					if (!conceptsFiltered.contains(concept)) conceptsFiltered.add(concept);
				}
			}
		}
		buf.append("\"n-meaning-classes-with-tree\":\"" + conceptsFiltered.size() + "\",\n");
		buf.append("\"meaning-classes-tree\":\"");
		buf.append(String.join(",", conceptsFiltered));
		buf.append("\",\n");
		
		
		
		// Cognate tree to meaning class rate mapping (mutationRate)
		buf.append("\"meaning-class-rates\":\"\n");
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			
			// Get the meaning class of this cogid
			int cogid_ = cogidUnique.get(i-1);
			
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && cogid[k] == cogid_) {
					String concept = conceptList[k];
					concept = "MC_" + concept.replaceAll("[ ,]", "_").replaceAll("/", "_");
					buf.append("<parameter id='Tree.t:" + i + "_mutationRate' spec='parameter.CompoundRealParameter' >\n");
					buf.append("\t<parameter idref='mutationRate.s:" + concept + "' />\n");
					buf.append("</parameter>\n");
					break;
				}
			}
			
		}
		buf.append("\",\n");
		
		
		// Cognate tree to meaning class rate mapping (tree)
		buf.append("\"pruned-trees\":\"\n");
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			
			// Get the meaning class of this cogid
			int cogid_ = cogidUnique.get(i-1);
			
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && cogid[k] == cogid_) {
					String concept = conceptList[k];
					concept = "MC_" + concept.replaceAll("[ ,]", "_").replaceAll("/", "_");
					buf.append("<tree id='CognateTree.t:" + i + "' spec='beast.phoneme.CognatePrunedTree' tree='@Tree.t:" + concept + "' birth='@cognateBirth.s:" + i + "' name='stateNode' estimate='true'>\n");
					buf.append("\t<data id='cognateSite." + i + "' spec='FilteredAlignment' userDataType='@CognateDataType' filter='" + i + "-" + i + "' data='@binsequences' />\n");
					buf.append("</tree>\n");
					break;
				}
			}
			
		}
		buf.append("\",\n");
		
		
		
		// Cognate tree to meaning class rate mapping (parallel operator)
		buf.append("\"parallel-distributions\":\"\n");
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			
			// Get the meaning class of this cogid
			int cogid_ = cogidUnique.get(i-1);
			
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && cogid[k] == cogid_) {
					String concept = conceptList[k];
					concept = "MC_" + concept.replaceAll("[ ,]", "_").replaceAll("/", "_");
					buf.append("<distribution spec='starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution' geneprior='@treePrior.t:" + concept + "' tree='@Tree.t:" + concept + "' treelikelihood='@CognateBirthLikelihood." + i + "'    id='@par_CognateBirthLikelihood." + i + "' />\n");
					buf.append("<distribution spec='starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution' geneprior='@treePrior.t:" + concept + "' tree='@Tree.t:" + concept + "' treelikelihood='@treeLikelihood.cognates." + i + "'   id='@par_treeLikelihood.cognates." + i + "' />\n");
					buf.append("<distribution spec='starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution' geneprior='@treePrior.t:" + concept + "' tree='@Tree.t:" + concept + "' treelikelihood='@treeLikelihood.vowels." + i + "'     id='@par_treeLikelihood.vowels." + i + "' />\n");
					buf.append("<distribution spec='starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution' geneprior='@treePrior.t:" + concept + "' tree='@Tree.t:" + concept + "' treelikelihood='@treeLikelihood.consonants." + i + "' id='@par_treeLikelihood.consonants." + i + "' />\n");
					break;
				}
			}
			
		}
		buf.append("\",\n");
		
		


		
		// Cognate trees and alignments
		buf.append("\"cognate-trees\":\"\n");
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			int j = i-1;
			buf.append("<tree id='Tree.t:" + i + "' spec='beast.evolution.tree.Tree' name='stateNode' estimate='true'>\n");
			buf.append("\t<taxonset id='taxonset-" + i + "' spec='TaxonSet'>\n");
			buf.append("\t\t<alignment id='cognate_" + i + "' spec='beast.phoneme.TaxonFilteredAlignment' filter='"+(j>0 ? startPhonemesCognate[j-1]+1 : 1)+"-"+startPhonemesCognate[j]+"' data='@data'>\n");
			buf.append("\t\t\t<taxonset spec='TaxonSet'>\n");
			for (String taxon : taxaPerCognate.get(j)) {
				buf.append("\t\t\t\t<taxon idref='" + taxon + "' />\n");
			}
			buf.append("\t\t\t</taxonset>\n");
			buf.append("\t\t</alignment>\n");
			buf.append("\t</taxonset>\n");
			buf.append("</tree>\n");
		}
		buf.append("\",\n");
		
		
		
		
		// 1 phoneme filter per cognate
		buf.append("\"phoneme-cognate-filters\":\"\n");
		for (int i = 2; i < cognateCount+1; i++) {
			if (!cognateNumberRange.contains("" + i)) continue;
			int j = i-1;
			buf.append("\t<alignment id='phoneme_cognate_" + i + "' spec='FilteredAlignment' filter='"+(j>0 ? startPhonemesCognate[j-1]+1 : 1)+"-"+startPhonemesCognate[j]+"' data='@data' />\n");
		}
		buf.append("\",\n");
		
		
		// Datatype for cognate filter
		buf.append("\"covarion-datatype\":\"\n");
		buf.append("<userDataType id='CognateDataType' spec='beast.evolution.datatype.TwoStateCovarion'/>");
		buf.append("\",\n");
		
		
		
		
		// Datatype for PD cognate filter
		buf.append("\"pd-covarion-datatype\":\"\n");
		buf.append("<userDataType id='CognateDataType' spec='beast.evolution.datatype.UserDataType' states='3' codelength='1' codeMap='A = 0,B = 1,C = 2,1 = 0 1,0 = 2,? = 0 1 2,- = 0 1 2' />\n");
		buf.append("\",\n");
		
		// One filter per meaning class (MC) - cognate
		buf.append("\"cognate-filters\":\"\n");
		for (int i = 0; i < concept.length; i++) {
				buf.append("<data id='MC_" + concept[i].replaceAll("[ ,]", "_").replaceAll("/", "_") +"_c' spec='FilteredAlignment' userDataType='@CognateDataType' filter='"+ (1+ (i>0 ? startCognates[i-1]+1 : 1))+"-"+ (1+startCognates[i])+"' data='@binsequences'/>\n");
		}		
		buf.append("\",\n");
		
		
		// One filter per meaning class - phonemes
		buf.append("\"phoneme-filters\":\"\n");
		for (int i = 0; i < concept.length; i++) {
				buf.append("<data id='MC_" + concept[i].replaceAll("[ ,]", "_").replaceAll("/", "_") +"_p' spec='FilteredAlignment' filter='"+(i>0 ? startPhonemesMC[i-1]+1 : 1)+"-"+startPhonemesMC[i]+"' data='@data'/>\n");
		}		
		buf.append("\"\n}\n");
		
		out.println(buf.toString());

	}
		

	private void standardiseTokens(String[] token) {
		boolean suppresWordGap = suppresWordGapInput.get();
		if (phonemeMapping.size() > 0 || suppresWordGap) {
			for (int k = 0; k < token.length; k++) {
				String string = token[k];
				if (string != null) {
					String [] strs = string.split("\\s");
					for (int i = 0; i < strs.length; i++) {
						if (phonemeMapping.containsKey(strs[i])) {
							strs[i] = phonemeMapping.get(strs[i]);
						}
					}
					
					StringBuilder b = new StringBuilder();
					b.append(strs[0]);
					for (int i = 1; i < strs.length; i++) {
						if (suppresWordGap && !strs[i].equals("+") && !strs[i].equals("_")) {
							b.append(' ');
							b.append(strs[i]);
						}
					}
					string = b.toString();
					token[k] = string;
				}
			}
		}
	}

	
	private void processList(File file, List<String> list) throws IOException {
		
		String s = BeautiDoc.load(file);
		String [] strs = s.split("\n");
		for (String str : strs) {
			if (!str.startsWith("#")) {
				list.add(str);
			}
		}
		
	}
	

	private void processMapping(File file, Map<String,String> phonemeMapping) throws IOException {
		String s = BeautiDoc.load(file);
		String [] strs = s.split("\n");
		for (String str : strs) {
			if (!str.startsWith("#")) {
				String [] strs2 = str.split("\t");
				if (strs2.length == 2) {
					phonemeMapping.put(strs2[0], strs2[1]);
				} else {
					Log.warning("found line with " + strs2.length + " columns in mapping file, where 2 are expected: " + str);
					Log.warning("the line is ignored");
				}
			}
		}
		
	}


	private void separateVowelsAndConsonants(String[] sequences) {
		int unhappyColumns = 0;
		for (int i = 0; i < sequences[0].length(); i+=2) {
			int isVowel = 0;
			int isCons = 0;
			for (int j = 0; j < sequences.length; j++) {
				char c = sequences[j].charAt(i);
				if (VOWEL.indexOf(c)>=0) {
					isVowel++;
				}
				if (CONSONANT.indexOf(c)>=0) {
					isCons++;
				}
			}
			if (isVowel>0 && isCons>0) {
				Map<String, Integer> set = new HashMap<>();
				for (int j = 0; j < sequences.length; j++) {
					String s = sequences[j].charAt(i) + "";
					if (!set.containsKey(s)) {
						set.put(s, 0);
					}
					set.put(s, set.get(s) + 1);
				}
				String str = "[";
				for (String s : set.keySet()) {
					str += s + "(" + set.get(s) + ") ";
				}
				str += "]";				
				Log.warning("Unhappy at column " + i + " " + isVowel + " " + isCons + " " + str);
				unhappyColumns++;
				
				// repair
				if (isVowel < isCons) {
					for (int j = 0; j < sequences.length; j++) {
						char c = sequences[j].charAt(i);
						if (VOWEL.indexOf(c)>=0) {
							sequences[j] = sequences[j].substring(0, i) + "-." + sequences[j].substring(i+2); 
						}
					}					
				} else {
					for (int j = 0; j < sequences.length; j++) {
						char c = sequences[j].charAt(i);
						if (CONSONANT.indexOf(c)>=0) {
							sequences[j] = sequences[j].substring(0, i) + "-." + sequences[j].substring(i+2); 
						}
					}					
				}
			}
		}
		
		Log.warning(unhappyColumns + " unhappy columns -- columns containing both consonants and vowels");
		if (unhappyColumns > 0) {
			Log.warning("Repaired unhappy columns");
			// sanity check: should have 0 unhappy columns
			separateVowelsAndConsonants(sequences);
		}
	}


	private void appendDataType(Set<String> phonemes, String id, StringBuilder buf) {
		String [] phonemes_ = phonemes.toArray(new String[]{});
		buf.append("\"datatype_"+id+"\":\"<userDataType id='" + id + "' spec='beast.phoneme.UserPhonemeDataType' states='" + phonemes_.length + "' codelength='2' codeMap='");
		Arrays.sort(phonemes_);
		int i = 0;
		for (String s : phonemes_) {
			buf.append(s.toUpperCase() + "=" + i);
			if (i < phonemes_.length-1) {
				buf.append(",");
			}
			i++;
		}
 		buf.append("'/>\",\n");
 		buf.append("\"gtrSymRates_" + id + "\":\"" + (phonemes.size() * (phonemes.size()-1)/2) + "\",\n");
 		buf.append("\"gtrAsymRates_" + id + "\":\"" + (phonemes.size() * (phonemes.size()-1)) + "\",\n");
		
	}

	private void appendDataTypeWithMissing(Set<String> phonemes, String id, StringBuilder buf, String codeToRemove, int countID) {
		phonemes.remove(codeToRemove);
		String [] phonemes_ = phonemes.toArray(new String[]{});
		buf.append("\"datatype_"+id+"_M" + countID + "\":\"<userDataType id='" + id + "WithMissing' spec='beast.phoneme.UserPhonemeDataType' states='" + phonemes_.length + "' codelength='2' codeMap='");
		Arrays.sort(phonemes_);
		int i = 0;
		for (String s : phonemes_) {
			buf.append(s.toUpperCase() + "=" + i);
			buf.append(",");
			i++;
		}
		buf.append("..=");
		for (i = 0; i < phonemes.size(); i++) {
			buf.append(i + " ");		
		}
		if (!phonemes.contains("-.")) {
			buf.append(",-.=");
			for (i = 0; i < phonemes.size(); i++) {
				buf.append(i + " ");			
			}
		}
		
 		buf.append("'/>\",\n");
 		buf.append("\"gtrSymRatesM" + countID + "_" + id + "\":\"" + (phonemes.size() * (phonemes.size()-1)/2) + "\",\n");
 		buf.append("\"gtrAsymRatesM" + countID + "_" + id + "\":\"" + (phonemes.size() * (phonemes.size()-1)) + "\",\n");		
	}
	
	private String processPhonemes(String sequence, Set<String> phonemes, Set<String> vowels, Set<String> consonants) {
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < sequence.length(); i += 2) {
			String phoneme = sequence.substring(i, i+2);
			if (encodingMapping.containsKey(phoneme)) {
				phoneme = encodingMapping.get(phoneme);
			}
			phonemes.add(phoneme);
			char c = phoneme.charAt(0);
			if (VOWEL.indexOf(c)>=0) {
				vowels.add(phoneme);
			} else {
				consonants.add(phoneme);
				if (phoneme.equals("..") || phoneme.equals("-.") || phoneme.equals("_.")) {
					vowels.add(phoneme);
				}
			}
			b.append(phoneme);
		}
		return b.toString();
	}

	
	private String cleanUp(String string) {
		string = string.replaceAll(" ",".");
		string = string.replaceAll("\\+","_");
		if (true) {
			return string;
		}
		// replace infrequently (<10) occurring phonemes by nearest phoneme
		string = string.replaceAll("ʰn","n.");
		string = string.replaceAll("ʲk","k.");		
		string = string.replaceAll("ʰs","s.");
		string = string.replaceAll("lʰ","l.");
		string = string.replaceAll("pʰ","p.");
		string = string.replaceAll("ᵐb","b.");

		string = string.replaceAll("tʃ","s.");
		string = string.replaceAll("ʰl","l.");
		string = string.replaceAll("kʰ","k.");
		string = string.replaceAll("ɣ","g");
		string = string.replaceAll("ʰm","m.");

		return string;
	}


	private int toCharSeqs(List<C> list, String[] docIds, StringBuilder [] charSeqs, int concept, boolean [][] docHasConcept) {
		if (list == null || list.size() == 0) {
			return 0;
		}
		int len = list.get(0).aligned.characters.length - 2;
		if (len < 0) {
			
		}
		String [] characters = list.get(0).aligned.characters;
		if (!characters[0].equals("x") || !characters[characters.length - 1].equals("x")) {
			System.err.println("Looks like the sequences i alignment");
		}
		
		// sanity check
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i).aligned.characters.length != len + 2) {
				throw new IllegalArgumentException("alignments do not have same length");
			}
		}
		
		
		int newLen = charSeqs[0].length() + 2 * len;
		boolean [] done = new boolean[charSeqs.length];
		for (C c : list) {
			int docId = indexOf(c.doculect, docIds);
			docHasConcept[docId][concept] = true;
			if (!done[docId]) {
				done[docId] = true;
				for (int i = 1; i < c.aligned.characters.length - 1; i++) {
					String s = c.aligned.characters[i];
					charSeqs[docId].append(s);
					if (s.length() < 2) {
						charSeqs[docId].append(' ');
					}
				}
			}
		}

		for (int i = 0; i < charSeqs.length; i++) {
			while (charSeqs[i].length() < newLen) {
				charSeqs[i].append(". ");
			}
		}
		return len;
	}


	public static void main(String[] args) throws Exception {
		new Application(new TSV2JSON(), "TSV to JSON converter", args);
	}
}

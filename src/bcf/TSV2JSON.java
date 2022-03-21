package bcf;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import beast.app.beauti.BeautiDoc;
import beast.app.util.Application;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.phoneme.UserPhonemeDataType;
import beast.util.ClusterTree;
import beast.util.ClusterTree.Type;

/**
 * Creates JSON file from word list
 * with aligned phoneme sequences
 *
 */
public class TSV2JSON extends TSV2Nexus {

	final public Input<File> mappingInput = new Input<>("mapping", "phoneme mapping used to preprocess phoneme strings to get rid of infrequently occurring phonemes. "
			+ "This is a tab delimited file with first column source phoneme and second column the target phoneme. "
			+ "All occurrances of source phonemes will be replaced by target phonemes. "
			+ "Ignored when not specified.");
	
	
	final public Input<Integer> ntreesInput = new Input<>("ntrees", "The number of word trees to assign meaning classes to ", 0);
	
	Map<String, String> phonemeMapping = new HashMap<>();
	
	@Override
	public void initAndValidate() {
		matrix = new DolgoScore();
		matrix = new SCAScore();
	}

	
	@Override
	public void run() throws Exception {
		if (tsvInput.get() == null || tsvInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("A valid TSV file must be specified");
		}
		if (mappingInput.get() != null) {
			processMapping();
		}
		
		TSVImporter importer = new TSVImporter(tsvInput.get(), languagesInput.get());
		
		String [] token = importer.getColumn("TOKENS");
		if (token == null ) {
			token = importer.getColumn("SEGMENTS");
		}
		standardiseTokens(token);
		int [] cogid = importer.getColumnAsInt("COGID");
		if (cogid == null ) {
			cogid = importer.getColumnAsInt("COGNACY");
		}
		String [] conceptList = importer.getColumn("CONCEPT");
		if (conceptList == null) {
			conceptList = importer.getColumn("PARAMETER_ID");
		}
		String [] doculect = importer.getColumn("DOCULECT");
		if (doculect == null) {
			doculect = importer.getColumn("LANGUAGE_ID");
		}
		Set<String> doculects = new LinkedHashSet<>();
		for (String s : doculect) {
			doculects.add(s);
		}
		Set<String> concepts = new LinkedHashSet<>();
		for (String s : conceptList) {
			if (s != null) {
				concepts.add(s);
			}
		}
		String [] concept = concepts.toArray(new String[]{});
		Arrays.sort(concept);
		
		int cognateCount = 0;
		for (int d : cogid) {
			cognateCount = Math.max(d, cognateCount);
		}
		cognateCount++;
		
		boolean [][] docHasConcept = new boolean[doculects.size()][concept.length];
		
		// group by cogid
		List<C> [] seqs = new List[cognateCount];
		for (int i = 0; i < cogid.length; i++) {
			int id = cogid[i];
			if (id > 0) {
				if (seqs[id] == null) {
					seqs[id] = new ArrayList<>();
				}
				seqs[id].add(new C(token[i], doculect[i]));
			}
		}
		
		// align individual cogids
		for (int i = 0; i < cognateCount; i++) {
			alignCognate(seqs[i], extendGapPenaltyInput.get(), openGapPenaltyInput.get());
			//alignCognate(seqs[i], 4.5f, -0.9f);
		}
		
		// convert cognate alignments to sequences
		StringBuilder [] charSeqs = new StringBuilder[doculects.size()];
		for (int i = 0; i < charSeqs.length; i++) {
			charSeqs[i] = new StringBuilder();
		}
		String [] docIds = doculects.toArray(new String[]{});
		int [] start = new int[concepts.size()];
		
		boolean [] done = new boolean[cognateCount];
		
		
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
			int len = 0;
			int cognatePerConcept = 0;
			int dataPerConcept = 0;
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k] != null && conceptList[k].equals(c)) {
					int i = cogid[k];
					if (!done[i]) {
						len +=  toCharSeqs(seqs[i], docIds, charSeqs, j, docHasConcept);
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
			start[j] = len + (j>0?start[j-1] : 0);
			j++;
		}
				
		
		Set<String> phonemes = new LinkedHashSet<>();
		Set<String> vowels = new LinkedHashSet<>();
		Set<String> consonants = new LinkedHashSet<>();
		
		String [] sequence = new String[docIds.length];
		for (int i = 0; i < docIds.length; i++) {
			sequence[i] = cleanUp(charSeqs[i].toString());
		}

		separateVowelsAndConsonants(sequence);

		
		// create binary cognate alignment
		char [][] cognateCharSeqs = new char[doculects.size()][cognateCount];
		for (int i = 0; i < cognateCharSeqs.length; i++) {
			Arrays.fill(cognateCharSeqs[i],'0');
		}
		for (int i = 0; i < cogid.length; i++) {
			int cogid_ = cogid[i];
			String doculect_ = doculect[i];
			if (doculect_ != null) {
				int k = indexOf(doculect_, docIds);
				cognateCharSeqs[k][cogid_] = '1';
			}
		}
		// check for constant columns
		boolean [] constantZeroColumn = new boolean[cognateCount];
		for (int k = 1; k < cognateCount; k++) {
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
				System.err.println("Constant column ("+cognateCharSeqs[0][k]+") in binary sequence at column " + (k+1));
			}
		}
		
		
	
		
		
		// output results
		StringBuilder buf = new StringBuilder();
		buf.append("{");
		buf.append("\"sequences\":\"\n");

		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getName());
			out = new PrintStream(outputInput.get());
		}

		
		
		for (int i = 0; i < docIds.length; i++) {
			if (docIds[i] != null) {
				buf.append("<sequence taxon='" + docIds[i] + "' ");
				if (docIds[i].length() < 25) {
					buf.append("                         ".substring(docIds[i].length()));
				}
				buf.append("value='");
				processPhonemes(sequence[i], phonemes, vowels, consonants);
				
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
					if (!constantZeroColumn[k]) {
						buf.append(cognateCharSeqs[i][k]);
					}
				}
				buf.append("'/>\n");
			}
		}
		buf.append("\",\n");

		buf.append("// _. = gap between words, .. = other cognate, -. = gap due to alignment\n");
		appendDataType(phonemes, "phonemes", buf);
		appendDataType(vowels, "vowels", buf);
		appendDataType(consonants, "consonants", buf);

		buf.append("\"words\":\"");
		for (int i = 0; i < concepts.size(); i++) {
			buf.append(concept[i].replaceAll("[ ,]", "_"));
			buf.append(concept[i].replaceAll("/", "-"));
			if (i < concepts.size() - 1) {
				buf.append(",");
			}
		}
		buf.append("\",\n");
		buf.append("\"words-1\":\"");
		for (int i = 1; i < concepts.size(); i++) {
			buf.append(concept[i].replaceAll("[ ,]", "_"));
			buf.append(concept[i].replaceAll("/", "-"));
			if (i < concepts.size() - 1) {
				buf.append(",");
			}
		}
		buf.append("\",\n");

		
		// Taxonset for multispecies coalescent
		buf.append("\"taxonset\":\"\n");
		for (int i = 0; i < docIds.length; i++) {
			buf.append("<taxon id='" + docIds[i] + ".lang' spec='TaxonSet'>\n");
			buf.append("\t<taxon id='" + docIds[i] + "' spec='Taxon' />\n");
			buf.append("</taxon>\n");
		}
		buf.append("\",\n");
		
		// List of meaning class IDs
		buf.append("\"meaning-classes\":\"");
		for (int i = 0; i < concepts.size(); i++) {
			buf.append("MC_" + concept[i].replaceAll("[ ,]", "_").replaceAll("/", "_"));
			if (i < concepts.size()-1) buf.append(",");
		}		
		buf.append("\",\n");
		
		
		
		
		
		
		// Meaning classes grouped by histories
		int ntrees = ntreesInput.get();
		if (ntrees > 0) {
			List<List<String>> historyClasses = getHistoryClasses(concept, docIds, sequence, start, phonemes, ntrees);
			
			
			// List of meaning class groups
			buf.append("\"grouped-meaning-classes\":\"");
			for (int i = 0; i < historyClasses.size(); i++) {
				buf.append("MC_" + (i+1));
				if (i < historyClasses.size()-1) buf.append(",");
			}		
			buf.append("\",\n");
			
			// One filter per meaning class (MC)
			buf.append("\"grouped-filters\":\"\n");
			for (int i = 0; i < historyClasses.size(); i++) {
				String filter = "";
				for (int c = 0; c < historyClasses.get(i).size(); c ++) {
					String concep = historyClasses.get(i).get(c);
					int conceptNum = 0;
					for (String concep2 : concepts) {
						if (concep2.equals(concep)) {
							break;
						}
						conceptNum++;
					}
					filter += (conceptNum>0 ? (start[conceptNum-1]+1) : 1)+"-"+start[conceptNum];
					if (c < historyClasses.get(i).size()-1) filter += ",";
				}
				buf.append("<data id='MC_"+ (i+1) + "' spec='FilteredAlignment' filter='"+ filter +"' data='@data'/>\n");
				
			}		
			buf.append("\"\n}\n");
			
		}
		
		
		// One filter per meaning class (MC)
		buf.append("\"filters\":\"\n");
		for (int i = 0; i < concepts.size(); i++) {
				buf.append("<data id='MC_" + concept[i].replaceAll("[ ,]", "_").replaceAll("/", "_") +"' spec='FilteredAlignment' filter='"+(i>0 ? (start[i-1]+1) : 1)+"-"+start[i]+"' data='@data'/>\n");
		}		
		buf.append("\"\n}\n");
		
		
		
		out.println(buf.toString());

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
			System.err.println(doculect[i] + " covered by " + k + " concepts");
		}
		System.err.println("On average: " + (mean / doculects.size()) + " (out of "
				+ concepts.size() + ") concepts covered per doculect\n");
		
		
		

		
		
		
		Log.warning("Done!");

	}



	private void standardiseTokens(String[] token) {
		if (phonemeMapping.size() > 0) {
			for (int k = 0; k < token.length; k++) {
				String string = token[k];
				String [] strs = string.split("\\s");
				for (int i = 0; i < strs.length; i++) {
					if (phonemeMapping.containsKey(strs[i])) {
						strs[i] = phonemeMapping.get(strs[i]);
					}
				}
				
				StringBuilder b = new StringBuilder();
				b.append(strs[0]);
				for (int i = 1; i < strs.length; i++) {
					b.append(' ');
					b.append(strs[i]);
				}
				string = b.toString();
				token[k] = string;
			}
		}
	}


	private void processMapping() throws IOException {
		String s = BeautiDoc.load(mappingInput.get());
		String [] strs = s.split("\n");
		for (String str : strs) {
			if (!str.startsWith("#")) {
				String [] strs2 = str.split("\t");
				if (strs2.length == 2) {
					if (strs2[0].length() == strs2[1].length()) {
						phonemeMapping.put(strs2[0], strs2[1]);
					} else {
						Log.warning("found line with source phoneme length not equal to target phoneme length " + str);
						Log.warning("the line is ignored");
					}
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
				if ("aeoiuy".indexOf(c)>=0) {
					isVowel++;
				}
				if ("bdfghjklmnprsttvwŋɢʃʔβ".indexOf(c)>=0) {
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
						if ("aeoiuy".indexOf(c)>=0) {
							sequences[j] = sequences[j].substring(0, i) + "-." + sequences[j].substring(i+2); 
						}
					}					
				} else {
					for (int j = 0; j < sequences.length; j++) {
						char c = sequences[j].charAt(i);
						if ("bdfghjklmnprsttvwŋɢʃʔβ".indexOf(c)>=0) {
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
		
	}


	private void processPhonemes(String sequence, Set<String> phonemes, Set<String> vowels, Set<String> consonants) {
		for (int i = 0; i < sequence.length(); i+= 2) {
			String phoneme = sequence.substring(i, i+2);
			phonemes.add(phoneme);
			char c = phoneme.charAt(0);
			if ("aeoiuy".indexOf(c)>=0) {
				vowels.add(phoneme);
			} else {
				consonants.add(phoneme);
				if (phoneme.equals("..") || phoneme.equals("-.") || phoneme.equals("_.")) {
					vowels.add(phoneme);
				}
			}
		}
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

	
	

	
	/*
	 * Return an assignment of meaning-classes into groups, based on their NJ tree similarities
	 */
	private List<List<String>> getHistoryClasses(String[] concepts, String [] docIds, String[] sequences, int[] starts, 
													Set<String> phonemes, int ntrees){
		
		
		int codelength = 2;
		
		
		// Datatype
		String codeMap = "";
		int i = 0;
		for (String s : phonemes) {
			codeMap = codeMap + s.toUpperCase() + "=" + i;
			if (i < phonemes.size()-1) {
				codeMap += ",";
			}
			i++;
		}
		UserPhonemeDataType dataType = new UserPhonemeDataType();
		//System.out.println(phonemes.size() + " " + codeMap);
		dataType.initByName("codelength", codelength, "states", phonemes.size(), "codeMap", codeMap);
		

		
		HashMap<String, String> conceptTopology = new LinkedHashMap<>();
		for (i = 0; i < concepts.length; i ++) {
			
			
			String concept = concepts[i];
			int start = (i>0 ? starts[i-1] : 0);
			int stop = starts[i] - 1;
			
			start = start*codelength;
			stop = (stop+1)*codelength;
			
			
			
			// Get alignment of subsequences
			List<Sequence> seqs = new ArrayList<>();
			for (int j = 0; j < sequences.length; j ++) {
				String label = docIds[j];
				String sequence = sequences[j].substring(start, stop);
				//Log.warning(concept + ": " + label + " has a sequence " + sequence);
				Sequence seq = new Sequence(label, sequence);
				seqs.add(seq);
			}
			Alignment alignment = new Alignment();
			alignment.initByName("sequence", seqs, "userDataType", dataType);
			
			
			// Build NJ tree
			ClusterTree tree = new ClusterTree();
			tree.initByName("clusterType", Type.upgma, "taxa", alignment);
			tree.getRoot().sort();
			String newick = tree.getRoot().toNewick(true);
			//Log.warning(concept + " : " + newick);
			conceptTopology.put(concept, newick);
		 	
		}
		
		
		// Count the topologies
		HashMap<String, Integer> topologies = new HashMap<>();
		for (String topology : conceptTopology.values()) {
			
			int count = 1;
			if (topologies.containsKey(topology)) {
				count = topologies.get(topology) + 1;
			}
			topologies.put(topology, count);
			
		}
		topologies = sortByValue(topologies);
		
		
		// Take the top ntrees topologies
		List<String> bestTopologies = new ArrayList<>();
		i = 0;
		for (String topology : topologies.keySet()) {
			if (i < ntrees) {
				bestTopologies.add(topology);
			}
			
			//Log.warning(topologies.get(topology) + " counts of " + topology);
			i ++;
			//if (i >= ntrees-1) break;
		}
		
		
		
		
		
		// Assign each concept to a class based off its topology. If it does not belong to a topology, then assign it at random
		List<List<String>> groups = new ArrayList<>();
		int randomAddTo = 0;
		for (i = 0; i < ntrees; i ++) groups.add(new ArrayList<>());
		for (String concept : concepts) {
			
			String topology = conceptTopology.get(concept);
			int which = bestTopologies.indexOf(topology);
			if (which == -1) {
				
				//Log.warning(concept + " being assigned to " + randomAddTo + " because it doesn't belong anywhere else");
				
				// Disperse it 
				groups.get(randomAddTo).add(concept);
				
				randomAddTo++;
				if (randomAddTo >= ntrees) randomAddTo = 0;
				
			}else {
				
				Log.warning(concept + " matched to " + which);
				groups.get(which).add(concept);
			}
			
		}
		
		
		return groups;
		
		
	}
	
	
	// Sort hashmap by value (decreasing order)
    private HashMap<String, Integer> sortByValue(HashMap<String, Integer> hm) {
        // Create a list from elements of HashMap
        List<Map.Entry<String, Integer> > list =
               new LinkedList<Map.Entry<String, Integer> >(hm.entrySet());
 
        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<String, Integer> >() {
            public int compare(Map.Entry<String, Integer> o1,
                               Map.Entry<String, Integer> o2)
            {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
         
        // put data from sorted list to hashmap
        HashMap<String, Integer> temp = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }
	
	
	
	

	private int toCharSeqs(List<C> list, String[] docIds, StringBuilder [] charSeqs, int concept, boolean [][] docHasConcept) {
		if (list == null) {
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

package bcf;

import java.io.PrintStream;
import java.util.*;

import beast.app.util.Application;
import beast.core.util.Log;

/**
 * Creates JSON file from word list
 * with aligned phoneme sequences
 *
 */
public class TSV2JSON extends TSV2Nexus {

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
		TSVImporter importer = new TSVImporter(tsvInput.get());
		
		String [] token = importer.getColumn("TOKENS");
		int [] cogid = importer.getColumnAsInt("COGID");
		String [] conceptList = importer.getColumn("CONCEPT");
		String [] doculect = importer.getColumn("DOCULECT");
		Set<String> doculects = new LinkedHashSet<>();
		for (String s : doculect) {
			doculects.add(s);
		}
		Set<String> concepts = new LinkedHashSet<>();
		for (String s : conceptList) {
			concepts.add(s);
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
			if (seqs[id] == null) {
				seqs[id] = new ArrayList<>();
			}
			seqs[id].add(new C(token[i], doculect[i]));			
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
				if (conceptList[k].equals(c)) {
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
			int k = indexOf(doculect_, docIds);
			cognateCharSeqs[k][cogid_] = '1';
		}
		// check for constant columns
		
		for (int k = 0; k < cognateCount; k++) {
			boolean isConstant = true;
			
			for (int i = 1; i < doculects.size(); i++) {
				if (cognateCharSeqs[i][k] != cognateCharSeqs[0][k]) {
					isConstant = false;
					break;
				}
			}
			if (isConstant) {
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
			buf.append("<sequence taxon='" + docIds[i] + "' ");
			if (docIds[i].length() < 25) {
				buf.append("                         ".substring(docIds[i].length()));
			}
			buf.append("value='");
			processPhonemes(sequence[i], phonemes, vowels, consonants);
			
			buf.append(sequence[i]);
			buf.append("'/>\n");
		}
		buf.append("\",\n");
		
		buf.append("\"binsequences\":\"\n");
		for (int i = 0; i < docIds.length; i++) {
			buf.append("<sequence taxon='" + docIds[i] + "' ");
			if (docIds[i].length() < 25) {
				buf.append("                         ".substring(docIds[i].length()));
			}
			buf.append("value='");
			for (int k = 0; k < cognateCount; k++) {
				buf.append(cognateCharSeqs[i][k]);
			}
			buf.append("'/>\n");
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

		buf.append("\"filters\":\"\n");
		for (int i = 0; i < concepts.size(); i++) {
				buf.append("<data id='" + concept[i].replaceAll("[ ,]", "_") +"' spec='FilteredAlignment' filter='"+(i>0 ? start[i-1] : 1)+"-"+start[i]+"' data='@data'/>\n");
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
		
		Log.warning(unhappyColumns + " unhappy columns");
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

	// replace infrequently (<10) occurring phonemes by nearest phoneme
	private String cleanUp(String string) {
		string = string.replaceAll("ʰn","n.");
		string = string.replaceAll("ʰs","s.");
		string = string.replaceAll("lʰ","l.");
		string = string.replaceAll("pʰ","p.");
		string = string.replaceAll("ᵐb","b.");

		string = string.replaceAll("tʃ","s.");
		string = string.replaceAll("ʰl","l.");
		string = string.replaceAll("kʰ","k.");
		string = string.replaceAll("ɣ","g");
		string = string.replaceAll("ʰm","m.");

		string = string.replaceAll(" ",".");
		string = string.replaceAll("\\+","_");

		return string;
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

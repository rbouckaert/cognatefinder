package bcf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.app.beauti.BeautiDoc;
import beast.core.Description;
import beast.core.util.Log;

@Description("import data from a tab separated file")
public class TSVImporter {
	String [] labels;
	String [][] datacolumns;
	
	public TSVImporter(File file) throws IOException {
		String data = BeautiDoc.load(file);
		process(data, null);
	}
	
	public TSVImporter(File file, File languageFile) throws IOException {
		String data = BeautiDoc.load(file);
		if (languageFile == null) {
			process(data, null);
			return;
		}
		Set<String> languages = new HashSet<>(); 
		for (String language: BeautiDoc.load(languageFile).split("\n")) {
			if (language != null && !language.startsWith("#") && language.trim().length() > 0) {
				languages.add(language);
			}
		}
		process(data, languages);
	}

	public TSVImporter(String data) {
		process(data, null);
	}

	private void process(String data, Set<String> languages) {
		String [] strs = data.split("\n");
		String separator = "\t";
		labels = strs[0].split(separator);
		if (labels.length == 1) {
			// try comma separates
			separator = ",";
			labels = strs[0].split(separator);
		}
		for (int i = 0; i < labels.length; i++) {
			labels[i] = labels[i].toUpperCase();
		}

		int n = 0;
		for (int i = 0; i < strs.length - 1; i++) {
			if (n == 0 || (!strs[i+1].startsWith("#") && matches(strs[i+1], languages, separator))) {
				n++;
			}
		}
		
		datacolumns = new String[labels.length][n];
		int k = 0;
		for (int i = 0; i < strs.length - 1; i++) {
			if (!strs[i+1].startsWith("#") && matches(strs[i+1], languages, separator)) {
				String [] s = strs[i+1].split(separator);
				if (s.length != labels.length) {
					if (s.length +1 == labels.length &&
						strs[i+1].length() > 0 && strs[i+1].charAt(strs[i+1].length()-1) == separator.charAt(0)) {
						// string ends in separator, so last entry is empty and not added to strs[i+1]
					} else {
					Log.warning("Line " + (i+2) + " does not have the same number of columns as the header: "
							+ "expected " + labels.length + " but got " + s.length);
					}
				}
				for (int j = 0; j < Math.min(labels.length, s.length); j++) {
					datacolumns[j][k] = s[j];
				}
				for (int j = s.length; j < labels.length; j++) {
					datacolumns[j][k] = null;
				}
				k++;
			}
		}		
	}
	
	
	private boolean matches(String str, Set<String> languages, String separator) {
		if (languages == null) {
			return true;
		}
		String [] strs = str.split(separator);
		for (String s : strs) {
			if (languages.contains(s)) {
				return true;
			}
		}
		return false;
	}

	String [] getColumn(int c) {
		return datacolumns[c];
	}
	
	String [] getColumn(String label) {
		for (int i = 0; i < labels.length; i++) {
			if (labels[i].equals(label)) {
				return getColumn(i);
			}
		}
		return null;
	}
	
	int [] getColumnAsInt(int c) {
		String [] data = datacolumns[c];
		int [] intData = new int[data.length];
		for (int i = 0; i < data.length; i++) {
			String str = data[i];
			if (str == null || str.trim().length() == 0) {
				intData[i] = Integer.MIN_VALUE;
			} else {
				try {
					intData[i] = Integer.parseInt(str);
				} catch (NumberFormatException e) {
					int k = str.length()-1;
					while (Character.isDigit(str.charAt(k))) {
						k--;
					}
					try {
						intData[i] = Integer.parseInt(str.substring(k+1));
					} catch (NumberFormatException e2) {
						Log.warning("Cannot determine int value of \"" + str + "\"");
						intData[i] = Integer.MIN_VALUE;
					}
				}
			}
		}
		return intData;
	}
	
	int [] getColumnAsInt(String label) {
		for (int i = 0; i < labels.length; i++) {
			if (labels[i].equals(label)) {
				return getColumnAsInt(i);
			}
		}
		return null;
	}
	
	
	/**
	 * Remove rows from the table
	 * @param remove
	 */
	public void removeRows(boolean[] remove) {

		
		// How many rows needed?
		int numRows = 0;
		for (int i = 0; i < remove.length; i ++) {
			if (!remove[i]) numRows++;
		}
		
		// Take subset
		String [][] datacolumns2 = new String[this.datacolumns.length][numRows];
		int newRowNum = 0;
		for (int oldRowNum = 0; oldRowNum < remove.length; oldRowNum ++) {
			if (!remove[oldRowNum]) {
				
				// Keep the row
				for (int colNum = 0; colNum < this.datacolumns.length; colNum ++) {
					datacolumns2[colNum][newRowNum] = datacolumns[colNum][oldRowNum];
				}
				newRowNum++;
			}
		}
		
		
		this.datacolumns = datacolumns2;
		
		
	}
	
	public static void main(String[] args) throws IOException {
//		TSVImporter io = new TSVImporter(new File("examples/mikronesian.tsv"));
//		
//		String [] value = io.getColumn("VALUE");
//		String [] form = io.getColumn("FORM");
//		for (int i = 0; i < value.length; i++) {
//			if (!value[i].equals(form[i])) {
//				System.out.println(value[i] + " " + form[i]);
//			}
//		}
//		int [] cognates = io.getColumnAsInt("COGID");

		//TSVImporter io = new TSVImporter(new File("examples/DravLex-2017-04-23.tsv"));
		//TSVImporter io = new TSVImporter(new File("examples/mikronesian.tsv"));
		TSVImporter io = new TSVImporter(new File("examples/polySmithWatermanGotoh.javanesian.tsv"));
		String [] CONCEPTS = io.getColumn("CONCEPT");
		String [] TOKENS = io.getColumn("TOKENS");
		String [] DOCULECTS = io.getColumn("DOCULECT");
		int [] COGIDS = io.getColumnAsInt("COGID");
		int n = COGIDS.length;
		int k = 1;
		String WORD = CONCEPTS[k];
		
		double sum = 0;
		int wordcount = 0;

		do {
			List<String> langs = new ArrayList<>();
			List<String[]> tokens = new ArrayList<>();
			List<Integer> labels = new ArrayList<>();
			Map<Integer,Set<String>> map = new LinkedHashMap<>();
	
			while (k < n && CONCEPTS[k].equals(WORD)) {
				String [] tokenlist = TOKENS[k].split(" ");
				for (int i = 0; i < tokenlist.length; i++) {
					tokenlist[i] = tokenlist[i].substring(0, 1);
					if (tokenlist[i].equals("ã")||tokenlist[i].equals("ā")) {
						tokenlist[i] = "a";
					}
					if (tokenlist[i].equals("d")) {
						tokenlist[i] = "t";
					}
					if (tokenlist[i].equals("ʋ")||tokenlist[i].equals("f")) {
						tokenlist[i] = "v";
					}
					if (tokenlist[i].equals("ē")) {
						tokenlist[i] = "e";
					}
					if (tokenlist[i].equals("ṣ")||tokenlist[i].equals("ʒ")||tokenlist[i].equals("ʃ")) {
						tokenlist[i] = "s";
					}
					if (tokenlist[i].equals("ṅ")||tokenlist[i].equals("ṉ")) {
						tokenlist[i] = "n";
					}
					if (tokenlist[i].equals("ṛ")||tokenlist[i].equals("ṟ")) {
						tokenlist[i] = "r";
					}
				}
				tokens.add(tokenlist);
				String lang = DOCULECTS[k];
				langs.add(lang);
				int label = COGIDS[k];
				labels.add(label);
				if (map.containsKey(label)) {
					map.get(label).add(lang);
				} else {
					Set<String> set = new LinkedHashSet<>();
					set.add(lang);
					map.put(label, set);
				}
				k++;
			}
			List<Set<String>> origSet = new ArrayList<>();
			for (int i = 0; i < labels.size(); i++) {
				Set<String> set = map.get(labels.get(i)); 
				origSet.add(set);
			}
	
			
			ThresholdCF cognateFinder = new ThresholdCF(3);
			List<Set<String>> set = cognateFinder.classifyTokens(langs, tokens);
			for (int i = 0; i < origSet.size(); i++) {
//				System.out.println(origSet.get(i) + " " + set.get(i));
			}
			double similarity = ThresholdCF.similarity(origSet, set);
			int w = 1;
			sum += similarity * w;
			wordcount+= w;
			System.out.println("similarity ("+WORD+")= " + similarity);
			if (k < n) {
				WORD = CONCEPTS[k];
			}
		} while (k < n);
		
		System.out.println("Mean: " + sum / wordcount);
	}

	
	

}

package bcf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.app.beauti.BeautiDoc;

public class ThresholdCF {

	int threshold = 2;
	
	ThresholdCF(int threshold) {
		this.threshold = threshold;
	}
	
	public List<Set<String>> classify(List<String> languages, List<String> words) {
		int n = words.size();
		int[][] dist = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				dist[i][j] = distance(words.get(i), words.get(j));
			}
		}
		
		int [] labels = new int[n];
		for (int i = 0; i < n; i++) {
			labels[i] = i;
		}
		
		boolean changed = false;
		do {
			changed = false;
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n; j++) {
					if (dist[i][j] < threshold && labels[i] != labels[j]) {
						if (labels[j] > labels[i]) {
							labels[j] = labels[i];
						} else {
							labels[i] = labels[j];
						}
						changed = true;
					}
				}
			}
		} while (changed);
		
		Map<Integer,Set<String>> map = new LinkedHashMap<>();
		for (int i = 0; i < n; i++) {
			String lang = languages.get(i);
			int label = labels[i];
			if (map.containsKey(label)) {
				map.get(label).add(lang);
			} else {
				Set<String> set = new LinkedHashSet<>();
				set.add(lang);
				map.put(label, set);
			}
		}
		
		List<Set<String>> lexemeMap = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			Set<String> set = map.get(labels[i]); 
			lexemeMap.add(set);
		}
		return lexemeMap;
	}

	// LevenshteinDistance from https://commons.apache.org/proper/commons-text/jacoco/org.apache.commons.text.similarity/LevenshteinDistance.java.html
	private static int distance(CharSequence left, CharSequence right) {
		if (left == null || right == null) {
			throw new IllegalArgumentException("Strings must not be null");
		}

		/*
		 * This implementation use two variable to record the previous cost
		 * counts, So this implementation use less memory than previous impl.
		 */

		int n = left.length(); // length of left
		int m = right.length(); // length of right

		if (n == 0) {
			return m;
		} else if (m == 0) {
			return n;
		}

		if (n > m) {
			// swap the input strings to consume less memory
			final CharSequence tmp = left;
			left = right;
			right = tmp;
			n = m;
			m = right.length();
		}

		int[] p = new int[n + 1];

		// indexes into strings left and right
		int i; // iterates through left
		int j; // iterates through right
		int upperLeft;
		int upper;

		char rightJ; // jth character of right
		int cost; // cost

		for (i = 0; i <= n; i++) {
			p[i] = i;
		}

		for (j = 1; j <= m; j++) {
			upperLeft = p[0];
			rightJ = right.charAt(j - 1);
			p[0] = j;

			for (i = 1; i <= n; i++) {
				upper = p[i];
				cost = left.charAt(i - 1) == rightJ ? 0 : 1;
				// minimum of cell to the left+1, to the top+1, diagonally left
				// and up +cost
				p[i] = Math.min(Math.min(p[i - 1] + 1, p[i] + 1), upperLeft + cost);
				upperLeft = upper;
			}
		}

		return p[n];
	}

	public static void main(String[] args) throws IOException {
		String WORD = "belly";
		//String data = BeautiDoc.load("/tmp/DravLex-2017-04-23.csv");
		String data = BeautiDoc.load("/tmp/mikronesian.tsv");
		String[] strs = data.split("\n");

		// column for {cognate class, language, word} 
		//int [] col = new int[]{1,2,3};
		int [] col = new int[]{12,1,6};
		
		WORD = "I";
		int k = 1;
		
		
		do {
			List<String> langs = new ArrayList<>();
			List<String> words = new ArrayList<>();
			List<Integer> labels = new ArrayList<>();
			Map<Integer,Set<String>> map = new LinkedHashMap<>();
	
			while (k < strs.length && strs[k].startsWith(WORD)) {
				String[] x = strs[k].split("\\s+");
				String word = x[col[2]];
				words.add(word);
				langs.add(x[col[1]]);
				int label = Integer.parseInt(x[col[0]]);
				labels.add(label);
				if (map.containsKey(label)) {
					map.get(label).add(word);
				} else {
					Set<String> set = new LinkedHashSet<>();
					set.add(word);
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
			List<Set<String>> set = cognateFinder.classify(langs, words);
			for (int i = 0; i < origSet.size(); i++) {
//				System.out.println(origSet.get(i) + " " + set.get(i));
			}
			double similarity = similarity(origSet, set);
			System.out.println("similarity ("+WORD+")= " + similarity);
			WORD = strs[k].substring(0, strs[k].indexOf('-'));
		} while (k < strs.length);
		
	}

	private static double similarity(List<Set<String>> origSet, List<Set<String>> altSet) {
		double sim = 0;
		for (int i = 0; i < origSet.size(); i++) {
			Set<String> org = origSet.get(i);
			Set<String> alt = altSet.get(i);
			Set<String> union = new HashSet<>();
			union.addAll(org);
			union.addAll(alt);
			Set<String> diff = new HashSet<>();
			diff.addAll(org);
			diff.retainAll(alt);
			sim += ((diff.size() + 0.0) / union.size());
		}
		return sim/origSet.size();
	}

}

package bcf;

import java.util.HashMap;
import java.util.Map;

import beast.core.Description;

/* Score matrices explained here:
 * https://doi.org/10.1038/nbt0804-1035
 * 
 * @article{eddy2004did,
  title={Where did the BLOSUM62 alignment score matrix come from?},
  author={Eddy, Sean R},
  journal={Nature biotechnology},
  volume={22},
  number={8},
  pages={1035--1036},
  year={2004},
  publisher={Nature Publishing Group}
}
 */
@Description("Score matrix for aligning phoneme strings ")
public abstract class Score {
	Map<String,Integer> charMap = new HashMap<>();
	double [][] score;
	
	protected void setScore(double [][] score) {
		this.score = score;
	}
	
	public double [][] getScore() {
		return score;
	}
	
	protected void process(int i, String string) {
		for (String s : string.split(",")) {
			s = s.trim().toLowerCase();
			charMap.put(s,  i);
			charMap.put(s+".",  i);
			charMap.put(s+"ː",  i);
			charMap.put(s+"ʰ",  i);
			charMap.put(s+"ⁿ",  i);
			charMap.put("ʰ"+s,  i);
			charMap.put("ⁿ"+s,  i);
			charMap.put("ʲ"+s,  i);
			s = s.toUpperCase();
			charMap.put(s,  i);
			charMap.put(s+"ː",  i);
			charMap.put(s+".",  i);
			charMap.put(s+"ʰ",  i);
			charMap.put(s+"ⁿ",  i);
			charMap.put("ʰ"+s,  i);
			charMap.put("ⁿ"+s,  i);
			charMap.put("ʲ"+s,  i);
		}		
	}


	public float score(String a, String b) {
		Integer i = charMapGet(a);
		Integer j = charMapGet(b);
		if (i == j) {
			if (a.equals(b)) {
				return (float) score[i][j];
			} else {
				return (float) score[i][j] - 0.01f;
			}
		}
		return (float) score[i][j];
	}
	
	private Integer charMapGet(String a) {
		Integer i = charMap.get(a);
		if (i == null) {
			charMap.put(a, 0);
			System.err.println("Character " + a + " was not mapped yet");
			return 0;
			//throw new IllegalArgumentException("Character " + a + " was not mapped yet");
		}
		return i;
	}

	public int getCode(String character) {
		Integer i = charMap.get(character);
		if (i == null) {
			return -1;
		}
		return i;
	}

	public float score(String [] a, String [] b, float gapPenalty) {
		if (a.length != b.length) {
			throw new IllegalArgumentException("lenght of sequences should be the same");
		}
		float score = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i].equals(SmithWatermanGotoh.Alignment.GAP)) {
				if (!b[i].equals(SmithWatermanGotoh.Alignment.GAP)) {
					score += gapPenalty;
				}
			} else {
				if (b[i].equals(SmithWatermanGotoh.Alignment.GAP)) {
					score += gapPenalty;
				} else {
					score += score(a[i], b[i]);
				}
			}
		}
		return score;
	}

}

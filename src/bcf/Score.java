package bcf;

public abstract class Score {

	abstract float score(String a, String b);

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

	public int getCode(String character) {
		return character.charAt(0);
	}
}

/*
 * $Id: SmithWatermanGotoh.java,v 1.18 2005/04/03 19:38:21 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package bcf;


import java.util.Arrays;

/**
 * An implementation of the Smith-Waterman algorithm with Gotoh's improvement
 * for biological local pairwise sequence alignment.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class SmithWatermanGotoh {
	public SmithWatermanGotoh() {
		super();
	}

	private String [] aligned1, aligned2;
	String [] getAligned1() {return aligned1;}
	String [] getAligned2() {return aligned2;}
	
	public abstract class Directions {
		/**
		 * Traceback direction stop
		 */
		public static final byte STOP = 0;
		/**
		 * Traceback direction left
		 */
		public static final byte LEFT = 1;
	    /**
		 * Traceback direction diagonal
		 */
		public static final byte DIAGONAL = 2;
		/**
		 * Traceback direction up
		 */
		public static final byte UP = 3;
	}

	public abstract class Markups {
		/**
		 * Markup line identity character
		 */
		public static final String IDENTITY	= "|";
		
		/**
		 * Markup line similarity character
		 */
		public static final String SIMILARITY	= ":";
		
		/**
		 * Markup line gap character
		 */
		public static final String GAP		= "-";
		
		/**
		 * Markup line mismatch character
		 */
		public static final String MISMATCH	= "x";
	}
	
	public class Alignment {

		/**
		 * Gap character
		 */
		public static final String GAP = "-";
	}
	

	/**
	 * Logger
	 */
	//private static final Logger logger = Logger
	//		.getLogger(SmithWatermanGotoh.class.getName());

	/**
	 * Aligns two sequences by Smith-Waterman algorithm
	 * 
	 * @param s1
	 *            sequene #1 ({@link Sequence})
	 * @param s2
	 *            sequene #2 ({@link Sequence})
	 * @param matrix
	 *            scoring matrix ({@link Matrix})
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return alignment object contains the two aligned sequences, the
	 *         alignment score and alignment statistics
	 * @see Sequence
	 * @see Matrix
	 */
	public String [] align(String [] s1, String [] s2, Score matrix,
			float o, float e) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();
		Score scores = matrix;

//		SmithWatermanGotoh sw = new SmithWatermanGotoh();

		int m = s1.length + 1;
		int n = s2.length + 1;

		byte[] pointers = new byte[m * n];

		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = Directions.STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = Directions.STOP;
		}

		short[] sizesOfVerticalGaps = new short[m * n];
		short[] sizesOfHorizontalGaps = new short[m * n];
		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}

		Cell cell = construct(s1, s2, scores, o, e, pointers,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);
		String [] markupLine = traceback(s1, s2, matrix, pointers, cell,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);
//		alignment.setName1(s1.getId());
//		alignment.setName2(s2.getId());
//		alignment.setMatrix(matrix);
//		alignment.setOpen(o);
//		alignment.setExtend(e);
//		logger.info("Finished in " + (System.currentTimeMillis() - start)
//				+ " milliseconds");
		return markupLine;
	}

	
	/**
	 * Calculate the score of the alignment, not using the score field (the
	 * function only uses sequence1, sequence2, matrix and gap penalties). (By:
	 * Bram Minnaert)
	 * 
	 * @return the calculated score
	 */
	public float calculateScore(int [] sequence1, int[] sequence2, float [][] matrix, float extend, float open) {
		final char GAP = '-';
		float calcScore = 0; // the calculated score
		boolean previous1wasGap = false; // in the previous step there was a gap
										 // in the first sequence
		boolean previous2wasGap = false; // in the previous step there was a gap
										 // in the second sequence
		int c1, c2; // the next character
		for (int i = 0, n = sequence1.length; i < n; i++) {
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == GAP) {
				if (previous1wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == GAP) {
				if (previous2wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else {
				calcScore += matrix[c1][c2];
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		return calcScore;
	}

	/**
	 * Constructs directions matrix for the traceback
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param matrix
	 *            scoring matrix
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return The cell where the traceback starts.
	 */
	private Cell construct(String [] s1, String [] s2, Score matrix, float o,
			float e, byte[] pointers, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();
		
		String[] a1 = s1;
		String[] a2 = s2;

		int m = s1.length + 1;
		int n = s2.length + 1;

		float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		float[] g = new float[n]; // score if xi aligns to a gap after yi
		Arrays.fill(g, Float.NEGATIVE_INFINITY);

		float h =  Float.NEGATIVE_INFINITY; // score if yi aligns to a gap after xi

		float[] v = new float[n]; // best score of alignment x1...xi to y1...yi
		Arrays.fill(v, 0f);

		float vDiagonal;


		float similarityScore, g1, g2, h1, h2;

		Cell cell = new Cell();

		int k = n;
		for (int i = 1; i < m; i++) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			int l = k + 1;
			for (int j = 1; j < n; j++) {
				similarityScore = matrix.score(a1[i - 1],a2[j - 1]);

				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - e;
				g2 = v[j] - o;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - e;
				h2 = v[j - 1] - o;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = Directions.STOP;
				} else if (v[j] == f) {
					pointers[l] = Directions.DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = Directions.UP;
				} else {
					pointers[l] = Directions.LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cell.getScore()) {
					cell.set(i, j, v[j]);
				}
				l++;
			}
			k += n;
		}
	//	logger.info("Finished in " + (System.currentTimeMillis() - start)
	//			+ " milliseconds");
		return cell;
	}

	/**
	 * Returns the alignment of two sequences based on the passed array of
	 * pointers
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param m
	 *            scoring matrix
	 * @param cell
	 *            The cell where the traceback starts.
	 * @return {@link Alignment}with the two aligned sequences and alignment
	 *         score.
	 * @see Cell
	 * @see Alignment
	 */
	private String [] traceback(String [] s1, String [] s2, Score m,
			byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();
		
		String[] a1 = s1;
		String[] a2 = s2;
		
		Score scores = m;

		int n = s2.length + 1;

		//Alignment alignment = new Alignment();
		//alignment.setScore(cell.getScore());

		int maxlen = s1.length + s2.length; // maximum length after the
												// aligned sequences

		String[] reversed1 = new String[maxlen]; // reversed sequence #1
		String[] reversed2 = new String[maxlen]; // reversed sequence #2
		String[] reversed3 = new String[maxlen]; // reversed markup

		int len1 = 0; // length of sequence #1 after alignment
		int len2 = 0; // length of sequence #2 after alignment
		int len3 = 0; // length of the markup line

		int identity = 0; // count of identitcal pairs
		int similarity = 0; // count of similar pairs
		int gaps = 0; // count of gaps

		String c1, c2;

		int i = cell.getRow(); // traceback start row
		int j = cell.getCol(); // traceback start col
		int k = i * n;

		boolean stillGoing = true; // traceback flag: true -> continue & false
								   // -> stop

		while (stillGoing) {
			switch (pointers[k + j]) {
			case Directions.UP:
				for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = a1[--i];
					reversed2[len2++] = Alignment.GAP;
					reversed3[len3++] = Markups.GAP;
					k -= n;
					gaps++;
				}
				break;
			case Directions.DIAGONAL:
				c1 = a1[--i];
				c2 = a2[--j];
				k -= n;
				reversed1[len1++] = c1;
				reversed2[len2++] = c2;
				if (c1 == c2) {
					reversed3[len3++] = Markups.IDENTITY;
					identity++;
					similarity++;
				} else if (scores.score(c1,c2) > 0) {
					reversed3[len3++] = Markups.SIMILARITY;
					similarity++;
				} else {
					reversed3[len3++] = Markups.MISMATCH;
				}
				break;
			case Directions.LEFT:
				for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = Alignment.GAP;
					reversed2[len2++] = a2[--j];
					reversed3[len3++] = Markups.GAP;
					gaps++;
				}
				break;
			case Directions.STOP:
				stillGoing = false;
			}
		}

//		alignment.setSequence1(reverse(reversed1, len1));
//		alignment.setStart1(i);
//		alignment.setSequence2(reverse(reversed2, len2));
//		alignment.setStart2(j);
//		alignment.setMarkupLine(reverse(reversed3, len3));
//		alignment.setIdentity(identity);
//		alignment.setGaps(gaps);
//		alignment.setSimilarity(similarity);

		//logger.info("Finished in " + (System.currentTimeMillis() - start)
		//		+ " milliseconds");
		aligned1 = reverse(reversed1, len1);
		aligned2 = reverse(reversed2, len2);
		//System.out.println(Arrays.toString(aligned1));
		//System.out.println(Arrays.toString(aligned2));
		return reverse(reversed3, len3);
	}
	

	/**
	 * Returns the maximum of 4 float numbers.
	 * 
	 * @param a
	 *            float #1
	 * @param b
	 *            float #2
	 * @param c
	 *            float #3
	 * @param d
	 *            float #4
	 * @return The maximum of a, b, c and d.
	 */
	private static float maximum(float a, float b, float c, float d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}

	/**
	 * Reverses an array of chars
	 * 
	 * @param a
	 * @param len
	 * @return the input array of char reserved
	 */
	private static String[] reverse(String[] a, int len) {
		String[] b = new String[len];
		for (int i = len - 1, j = 0; i >= 0; i--, j++) {
			b[j] = a[i];
		}
		return b;
	}
	
	
	
	
	public static void main(String[] args) {
		String [] w1 = "x t ʃ a: ɾ a m x".split(" ");
		String [] w2 = "x s a: m b a l x".split(" ");

		w2 = "_ a a\" a: a:ʰ b bʰ c d dʰ e e: f g gʰ h i i: j k kʰ l m n o o: p p: r s t tʰ u u: v x y z ã ñ ā ē ī ō ū Ɛ ɖ ə ɭ ɲ ɳ ɻ ɽ ɾ ʂ ʃ ʈ ʋ ʒ ḍ ḻ ṅ ṇ ṉ ṛ ṟ ṣ ṭ ə: ʃ: ɖʰ ".split(" ");
		
		w1 = "x, b, a, ŋ, x".split(", ");
		w2 = "x, e, β, a, x".split(", ");
		w1 = "x, x, x, b, a, ŋ, x, x, x".split(", ");
		w2 = "x, x, x, e, β, a, x, x, x".split(", ");
		
		// SmithWaterman align = new SmithWaterman(sub, d)
		SCAScore matrix = new SCAScore();
		SmithWatermanGotoh sw = new SmithWatermanGotoh();
		String [] w3 = sw.align(w1, w2, matrix, 5f, -0.8f);
		System.out.println(Arrays.toString(w3));
		System.out.println(Arrays.toString(w1));
		System.out.println(Arrays.toString(w2));
	}

}
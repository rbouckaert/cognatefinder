package bcf;

import java.util.Arrays;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.Node;
import beast.util.ClusterTree;

public class MultiAligner {
	SmithWatermanGotoh pairAligner;
	float o, e;
	Score matrix;
	
	public Seq [] align(String [] ids, String [][] seqs, Score matrix,
			float openGapPenaly, float extendGapPenalty) {
		int n = seqs.length;
		this.o = openGapPenaly;
		this.e = extendGapPenalty;
		this.matrix = matrix;
		
		// calculate distance matrix
		double [][] distances = new double[n][n];
		pairAligner = new SmithWatermanGotoh();
		for (int i = 0; i < n ; i++) {
			for (int j = i + 1; j < n; j++) {
				pairAligner.align(seqs[i], seqs[j], matrix, openGapPenaly, extendGapPenalty);
				distances[i][j] = 1.0/matrix.score(pairAligner.getAligned1(), pairAligner.getAligned2(), o);
				distances[j][i] = distances[i][j];
			}
		}
		
		Distance dist = new Distance() {			
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				return distances[taxon1][taxon2];
			}
		};
		
		// create tree		
		ClusterTree tree = new ClusterTree();
		TaxonSet taxonset = new TaxonSet();
		for (String str : ids) {
			taxonset.taxonsetInput.get().add(new Taxon(str));
		}
		taxonset.initAndValidate();
		tree.initByName("distance", dist, "taxonset", taxonset);
		
		//System.out.println("Guiding tree:" + tree.getRoot().toNewick());
		
		// traverse tree
		Seq [] alignedSeqs = new Seq[tree.getNodeCount()];
		for (int i = 0; i < n; i++) {
			alignedSeqs[i] = new Seq(seqs[i], matrix);
		}
		traverseUp(alignedSeqs, tree.getRoot());
		//System.err.println(Arrays.toString(alignedSeqs));
		traverseDown(alignedSeqs, tree.getRoot());
		//System.err.println(Arrays.toString(alignedSeqs));

		int len = 0;
		for (int i = 0; i < alignedSeqs.length; i++) {
			if (alignedSeqs[i].characters.length > len) {
				len = alignedSeqs[i].characters.length;
			}
		}
		for (int i = 0; i < alignedSeqs.length; i++) {
			alignedSeqs[i].makeLenght(len);
		}

		
		
		//System.out.println("Alignment:");
		//for (int i = 0; i < n; i++) {
		//	System.out.println(alignedSeqs[i].toString());
		//}
		return alignedSeqs;
	}
	
	
	private boolean unequalLengths(Seq[] alignedSeqs) {
		int len = alignedSeqs[0].characters.length;
		for (int i = 1; i < alignedSeqs.length; i++) {
			if (alignedSeqs[i].characters.length != len) {
				return true;
			}
		}
		return false;
	}


	private void traverseUp(Seq[] alignedSeqs, Node node) {
		if (!node.isLeaf()) {
			traverseUp(alignedSeqs, node.getLeft());
			traverseUp(alignedSeqs, node.getRight());
			pairAligner.align(alignedSeqs[node.getLeft().getNr()].characters, alignedSeqs[node.getRight().getNr()].characters, matrix, o, e);
			
			String [] result = pairAligner.getAligned1();
			String [] result2 = pairAligner.getAligned2();
			for (int i = 0; i < result.length; i++) {
				if (result[i].equals("-")) {
					result[i] = result2[i];
				}
			}
			
			alignedSeqs[node.getNr()] = new Seq(result, matrix);
		}		
	}

	private void traverseDown(Seq[] alignedSeqs, Node node) {
		if (!node.isLeaf()) {
			pairAligner.align(alignedSeqs[node.getLeft().getNr()].characters, alignedSeqs[node.getNr()].characters, matrix, o, e);
			alignedSeqs[node.getLeft().getNr()] = new Seq(pairAligner.getAligned1(), matrix);
			
			pairAligner.align(alignedSeqs[node.getRight().getNr()].characters, alignedSeqs[node.getNr()].characters, matrix, o, e);
			alignedSeqs[node.getRight().getNr()] = new Seq(pairAligner.getAligned1(), matrix);

			traverseDown(alignedSeqs, node.getLeft());
			traverseDown(alignedSeqs, node.getRight());
		}
	}



	public static void main(String[] args) {
			
			
			String [] s = new String[]{"p a ʈ ʈ e", 
					"a ʈ ʈ e",
					"p a ṭ ṭ a", 
					"p a ʈ ʈ e",
					"p a ṭ ṭ a"};
			
			s = new String[]{"ʋ a j a ɾ ə ",
					"v a y i ṟ u ",
					"p i ɾ ",
					"b a n d ʒ i ",
					"b a i ɾ u ",
					"b a: ɾ ə "
					
			};
s = new String[]{
		"k a ɖ i",
		"k a a t _ k a n n i n g ",
		"k a ḍ i ",
		"k a ɖ i ",
		"k a ɖ ə ",
		"k a d i k k u g a ",
		"k a ɖ i ",
		"k a ɖ i "
	};

s = new String[] {"a u",
		"a u",
		"a u",
		"a u",
		"a ʔ u",
		"a u",
		"v a u",
		"o u",
		"a u",
		"k a u",
		"a u",
		"a v a u",
		"v a u",
		"a u",
		"k a u",
		"k a u",
		"a u",
		"a u",
		"a ŋ a u",
		"a u",
		"a u",
		"a u",
		"k a u",
		"k a u",
		"k o a u",
		"a u",
		"a u",
		"k o a u",
		"a u",
		"a k u",
		"k o u",
		"a β a u",
		"n a u",
		"v a u"};
s = new String[] {
//		"v a ʔ u",
//		"w a l u",
//		"v a r u",
//		"w a r u",
//		"v a l u",
//		"v a ʔ u",
//		"v a l u",
//		"v a l u",
//		"v a r u",
//		"v a r u",
//		"v a ʔ u",
		"b a ŋ",
		"e β a",
		"v a l u",
		"w a l u",
		"v a r u",
		"v a l u",
		"v a r u",
		"v a l u",
		"v a r u",
		"v a l u",
		"v a r u",
		"v a l u",
		"v a ʔ u",
		"v a l u",
		"v a l u",
		"β a r u",
		"v a r u",
		"v a ɢ u"
		};

	s = new String[]{"r i m a _ ŋ a u r u",
			"l i m a g a f u l u",
			"l i m a s e f u l u",
			"n i m a ŋ o f u l u",
			"l i m a ŋ a u l u",
			"r i m a _ r o ŋ o ʔ u r u"
	};

	s = new String[]{
			"a _ ŋ a u",
			"a g a f u",
//			"a s e f u",
//			"a ŋ o f u",
//			"a ŋ a u l",
			"a _ r o ŋ o ʔ u"
	};


	SCAScore matrix = new SCAScore();
	
	System.err.println(matrix.score("ŋ","ŋ"));
	System.err.println(matrix.score("ŋ","a"));
	System.err.println(matrix.score("ŋ","g"));
	System.err.println(matrix.score("g","a"));
	
	Seq [] seqs = multiAlign(s, matrix, 5f, -1f);
	for (Seq seq : seqs) {
		System.out.println(Arrays.toString(seq.characters));
	}
			
//			String [] w1 = "x t ʃ a: ɾ a m x".split(" ");
//			String [] w2 = "x s a: m b a l x".split(" ");
//			String [] w3 = "x t ʃ a: m b a l x".split(" ");
//			// SmithWaterman align = new SmithWaterman(sub, d)
//			SCAScore matrix = new SCAScore();
//			MultiAligner sw = new MultiAligner();
//			sw.align(new String[]{"lang1", "lang2", "lang3"},
//					new String[][] {w3, w2, w1}, matrix, 8f, -1f);
	}


	public static Seq [] multiAlign(String[] s, Score matrix, float openGapPenalty, float extendGapPenalty) {
		String[][] seqs = new String[s.length][];
		String [] ids = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			seqs[i] = ("x " + s[i].trim() + " x").split(" ");
			ids[i] = "id" + (i+1);
		}
		
		MultiAligner sw = new MultiAligner();
		return sw.align(ids, seqs, matrix, openGapPenalty, extendGapPenalty);
	}

}

package bcf;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;

public class TSV2Nexus extends Runnable {
	public Input<File> tsvInput = new Input<>("tsv", "TSV file with cognate data",
			new File("file.tsv"));
	public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("[[none]]"));
	public Input<Float> extendGapPenaltyInput = new Input<>("egp","extend gap penalty used in aligning sequences. "
			+ "These can be tuned to minimise the nr of unhappy (containing both vowels and consonants) columns", 4.0f); 
	public Input<Float> openGapPenaltyInput = new Input<>("ogp","open gap penalty used in aligning sequences. "
			+ "These can be tuned as togetehr with egp input", -0.75f); 
	public Input<File> languagesInput = new Input<>("languages", "text file with languages to include, one per line. If not specified, all languages are included");
	
	
	Score matrix;

	@Override
	public void initAndValidate() {
		matrix = new SCAScore();
	}

	class C {
		String token;
		String doculect;
		Seq aligned;
		
		C(String token, String doculect) {
			this.token = token;
			this.doculect = doculect;
		}
		@Override
		public String toString() {
			if (aligned != null) {
				return aligned.toString();
			}
			return token;
		}
	}
	
	@Override
	public void run() throws Exception {
		if (tsvInput.get() == null || tsvInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("A valid TSV file must be specified");
		}
		TSVImporter importer = new TSVImporter(tsvInput.get(), languagesInput.get());
		
		String [] token = importer.getColumn("TOKENS");
		int [] cogid = importer.getColumnAsInt("COGID");
		String [] doculect = importer.getColumn("DOCULECT");
		Set<String> doculects = new LinkedHashSet<>();
		for (String s : doculect) {
			doculects.add(s);
		}
		int cognateCount = 0;
		for (int d : cogid) {
			cognateCount = Math.max(d, cognateCount);
		}
		cognateCount++;
		
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
		}
		
		// convert cognate alignments to sequences
		StringBuilder [] charSeqs = new StringBuilder[doculects.size()];
		for (int i = 0; i < charSeqs.length; i++) {
			charSeqs[i] = new StringBuilder();
		}
		String [] docIds = doculects.toArray(new String[]{});
		int len = 0;
		for (int i = 0; i < cognateCount; i++) {
			len +=  toCharSeqs(seqs[i], docIds, charSeqs);
		}
		
		// convert character sequences to standard sequences
		charSeqs = standardiseCharSeqs(charSeqs);
		
		
		// output results
		StringBuilder buf = new StringBuilder();
		buf.append("#NEXUS\n");
		buf.append("BEGIN DATA;\n");
		buf.append("DIMENSIONS NTAX=" + doculects.size() + " NCHAR=" + len + ";\n");
		buf.append("FORMAT DATATYPE=STANDARD MISSING=? GAP=-  SYMBOLS=\"0123456789\";\n");
		buf.append("MATRIX\n");

		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getName());
			out = new PrintStream(outputInput.get());
		}

		for (int i = 0; i < docIds.length; i++) {
			buf.append(docIds[i] + " ");
			if (docIds[i].length() < 25) {
				buf.append("                         ".substring(docIds[i].length()));
			}
			buf.append(charSeqs[i].toString());
			buf.append("\n");
		}
		buf.append(";\nend;\n");
		
		out.println(buf.toString());

		Log.warning("Done!");

	}

	private StringBuilder[] standardiseCharSeqs(StringBuilder[] charSeqs) {
		StringBuilder[] newSeqs= new StringBuilder[charSeqs.length];
		for (int i = 0; i < newSeqs.length; i++) {
			newSeqs[i] = new StringBuilder();
		}
		int n = charSeqs[0].length()/2;
		for (int i = 0; i < n; i++) {
			Map<String, Integer> map = new HashMap<>();
			map.put("- ", 0);
			for (int j = 0; j < charSeqs.length; j++) {
				String s = charSeqs[j].substring(i*2, i*2+2);
				if (!s.equals("? ")) {
					if (!map.containsKey(s)) {
						map.put(s, map.size());
					}
				}
			}
			if (map.size() > 1) {
	 			for (int j = 0; j < charSeqs.length; j++) {
					String s = charSeqs[j].substring(i*2, i*2+2);
					if (s.equals("? ")) {
						newSeqs[j].append("0");
					} else {
						newSeqs[j].append(1+map.get(s));
					}
					//newSeqs[j].append(" ");
				}
			}
		}

		return newSeqs;
	}

	private int toCharSeqs(List<C> list, String[] docIds, StringBuilder [] charSeqs) {
		if (list == null) {
			return 0;
		}
		int len = list.get(0).aligned.characters.length - 2;
		
		// sanity check
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i).aligned.characters.length != len + 2) {
				throw new IllegalArgumentException("alignments do not have same length");
			}
		}
		
		
		int newLen = charSeqs[0].length() + 2 * len;
		for (C c : list) {
			int docId = indexOf(c.doculect, docIds);
			for (int i = 1; i < c.aligned.characters.length - 1; i++) {
				String s = c.aligned.characters[i];
				charSeqs[docId].append(s);
				if (s.length() < 2) {
					charSeqs[docId].append(' ');
				}
			}
		}

		for (int i = 0; i < charSeqs.length; i++) {
			while (charSeqs[i].length() < newLen) {
				charSeqs[i].append("? ");
			}
		}
		return len;
	}

	protected int indexOf(String doculect, String[] docIds) {
		for (int i = 0; i < docIds.length; i++) {
			if (docIds[i] != null && docIds[i].equals(doculect)) {
				return i;
			}
		}
		throw new IllegalArgumentException("Cannot find doculect " + doculect + " in docIds");
	}

	protected void alignCognate(List<C> list, float openGapPenalty, float extendGapPenalty) {
		if (list == null || list.size() == 0) {
			return;
		}
		if (list.size() == 1) {
			C c = list.get(0);
			c.aligned = new Seq(c.token.split(" "), matrix);
		}
		String [] s = new String[list.size()];
		for (int i = 0; i < s.length; i++) {
			s[i] = list.get(i).token;
		};
		Seq [] seqs = MultiAligner.multiAlign(s, matrix, openGapPenalty, extendGapPenalty);
		for (int i = 0; i < s.length; i++) {
			list.get(i).aligned = seqs[i];
		}
	}

	public static void main(String[] args) throws Exception {
		new Application(new TSV2Nexus(), "TSV to Nexus converter", args);
	}
}

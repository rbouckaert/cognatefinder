package bcf;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;

/**
 * Creates JSON file from word list
 * with aligned phoneme sequences
 *
 */
public class TSV2JSON extends TSV2Nexus {

	@Override
	public void initAndValidate() {
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
			alignCognate(seqs[i]);
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
		for (String c : concept) {
			int len = 0;
			for (int k = 0; k < conceptList.length; k++) {
				if (conceptList[k].equals(c)) {
					int i = cogid[k];
					if (!done[i]) {
						len +=  toCharSeqs(seqs[i], docIds, charSeqs);
						done[i] = true;
					}
				}
			}
			start[j] = len + (j>0?start[j-1] : 0);
			j++;
		}
//		int len = 0;
//		for (int i = 0; i < cognateCount; i++) {
//			len +=  toCharSeqs(seqs[i], docIds, charSeqs);
//		}
				
		
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
			buf.append(cleanUp(charSeqs[i].toString()));
			buf.append("'/>\n");
		}
		buf.append("\",\n");
		buf.append("\"filters\":\"\n");
		for (int i = 0; i < concepts.size(); i++) {
				buf.append("<data id='" + concept[i] +"' spec='FilteredAlignment' filter='"+(i>0 ? start[i-1] : 1)+"-"+start[i]+"'/>\n");
		}		
		buf.append("\"\n}\n");
		
		out.println(buf.toString());

		Log.warning("Done!");

	}

	private Object cleanUp(String string) {
		string = string.replaceAll("ʰn","n_");
		string = string.replaceAll("ʰs","s_");
		string = string.replaceAll("lʰ","l_");
		string = string.replaceAll("pʰ","p_");
		string = string.replaceAll("ᵐb","b_");

		string = string.replaceAll("tʃ","s_");
		string = string.replaceAll("ʰl","l_");
		string = string.replaceAll("kʰ","k_");
		string = string.replaceAll("ɣ","g");
		string = string.replaceAll("ʰm","m_");

		string = string.replaceAll(" ","_");
		string = string.replaceAll("\\+","_");

		return string;
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
				charSeqs[i].append("O ");
			}
		}
		return len;
	}


	public static void main(String[] args) throws Exception {
		new Application(new TSV2JSON(), "TSV to JSON converter", args);
	}
}

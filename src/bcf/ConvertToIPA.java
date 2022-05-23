package bcf;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;

@Description("Convert a list of phonemes separated by '/', or with modifiers (e.g. superscripts, long sounds) into IPA, or give a warning when not possible")
public class ConvertToIPA extends Runnable {
	
	public Input<File> encodingInput = new Input<>("encoding", "phoneme map to something resembling IPA. This is tab delimited with two columns");
	public Input<File> tsvInput = new Input<>("tsv", "TSV file with cognate data", Input.Validate.REQUIRED);
	public Input<File> languagesInput = new Input<>("languages", "text file with languages to include, one per line. If not specified, all languages are included");
	public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified", new OutFile("[[none]]"));
	
	Map<String, String> encodingMapping = new HashMap<>();
	
	// https://www.internationalphoneticassociation.org/IPAcharts/IPA_chart_orig/IPA_charts_E.html
	// 2020
	// I converted some of the additional symbols into 2 letter forms (eg. ʦ -> ts), and I replaced the non ascii 'ɡ' with the ascii 'g'
	final static String[] IPAs =
			("p b t d ʈ ɖ c ɟ k g q ɢ ʔ m ɱ n ɳ ɲ ŋ ɴ ʙ r ʀ ⱱ ɾ ɽ ɸ β f v θ ð s z ʃ ʒ ʂ ʐ ç ʝ x ɣ χ ʁ ħ ʕ h ɦ ɬ ɮ ʋ ɹ ɻ j ɰ l ɭ ʎ ʟ ɫ " + // Pulmonic
			"ʘ ǀ ǃ ǂ ǁ ɓ ɗ ʄ ɠ ʛ ʼ pʼ tʼ kʼ sʼ " + // Non-pulmonic
			"ʍ w ɥ ʜ ʢ ʡ ɕ ʑ ɺ ɧ " + // Other
			"i y ɪ ʏ e ø ɛ œ æ a ɶ ɨ ʉ ɘ ɵ ə ɜ ɞ ɐ ɯ u ʊ ɤ o ʌ ɔ ɑ ɒ ɚ " + // Vowels
			"ts dz tʃ dʒ ʨ ʥ ɝ".toLowerCase()).split(" ") ; // Additional symbols
	
	
	@Override
	public void initAndValidate() {
		
		
		
		
	}

	@Override
	public void run() throws Exception {
		
		
		Log.warning("Converting to IPA...");
		
		
		// Load in data
		if (encodingInput.get() != null) TSV2JSON.processMapping(encodingInput.get(), encodingMapping, true);
		TSVImporter importer = new TSVImporter(tsvInput.get(), languagesInput.get());
		
		
		Log.warning("Loaded data");
		
		// Word tokens
		String [] tokens = importer.getColumn("TOKENS");
		if (tokens == null ) {
			tokens = importer.getColumn("SEGMENTS");
		}
		
		
		Log.warning("Found tokens");
		
		// Get list of phonemes and their counts
		HashMap<String, Integer> phonemes = new HashMap<>();
		for (String token : tokens) {
			
			if (token == null || token.isEmpty()) continue;
			
			String[] phonemesToken = token.trim().split(" ");
			for (String p : phonemesToken) {
				
				
				if (p.isEmpty()) continue;
				

				// Skip special characters
				if (p.equals("+") || p.equals("_") || p.equals("-") || p.equals(".")) continue;
				
				
				
				int count = 0;
				if (phonemes.containsKey(p)) {
					count = phonemes.get(p) + 1;
				}else {
					count = 1;
				}
				phonemes.put(p, count);
			}
			
		}
		phonemes = sortTableByCount(phonemes);
		
		
		Log.warning("Detecting phonemes:");
		for (String phoneme : phonemes.keySet()) {
			Log.warning("Phoneme '" + phoneme + "' (" + phonemes.get(phoneme) + " counts)");
		}
		Log.warning("\n");
		
		
		// Conversion
		Log.warning("Mapping phonemes to IPA:");
		HashMap<String, String> phonemeMapping = new HashMap<>();
		int errorCount = 0;
		List<String> unMapped = new ArrayList<>();
		for (String phonemeDataset : phonemes.keySet()) {
			
			
			String phoneme = phonemeDataset;
			if (encodingMapping != null && encodingMapping.containsKey(phonemeDataset)) {
				phoneme = encodingMapping.get(phonemeDataset);
				Log.warning("User mapping from '" + phonemeDataset + "' to '" + phoneme + "'");
			}
			
			
			//Log.warning(phoneme);
			// Length check
			//if (phoneme.length() > 2) {
				//Log.warning("Error: '" + phoneme + "' is more than 2 digits! All phonemes should be 1 or 2 digits!");
				//errorCount++;
			//}
			
			
			if (isIPA(phoneme, true)) {
				Log.warning("Success: '" + phoneme + "' is already IPA");
				phonemeMapping.put(phonemeDataset, phoneme);
			}else {
				
				
				// Split by '/'
				String[] bits = phoneme.split("/");
				if (phoneme.length() == 2 && bits.length == 1) {
					//bits = phoneme.split("");
				}
				if (bits.length == 2) {
					
					boolean firstIsIPA = isIPA(bits[0], false);
					boolean secondIsIPA = isIPA(bits[1], false);
					
					// Map to first 
					if (firstIsIPA && !secondIsIPA) {
						phonemeMapping.put(phoneme, bits[0]);
						Log.warning("Success: '" + phonemeDataset + "' is being mapped to '" + bits[0] + "'");
						continue;
					}
					
					// Map to second
					else if (!firstIsIPA && secondIsIPA) {
						phonemeMapping.put(phoneme, bits[1]);
						Log.warning("Success: '" + phonemeDataset + "' is being mapped to '" + bits[1] + "'");
						continue;
					}
					
					// Map to most probable
					else if (firstIsIPA && secondIsIPA) {
						
						int count1 = phonemes.containsKey(bits[0]) ? phonemes.get(bits[0]) : 0;
						int count2 = phonemes.containsKey(bits[1]) ? phonemes.get(bits[1]) : 0;
						String mapped = count1 >= count2 ? bits[0] : bits[1];
						
						phonemeMapping.put(phoneme, mapped);
						Log.warning("Success: '" + phonemeDataset + "' is being mapped to '" + mapped + "' because it has a larger frequency (" + count1 + "/" + count2 + ")");
						continue;
					}
					
					// Neither are IPA.
					else {
						
						
						
						
						
					}
					
				}
				Log.warning("Error: '" + phoneme + "' is not IPA!");
				unMapped.add(phonemeDataset);
				errorCount++;
			}
			
			
			
			
		}
		Log.warning("\n");
		
		
		// Quality check: make sure that either colons OR diamonds are used, but not both
		
		
		

		
		// Sort by vaue
		phonemeMapping = sortTableByValue(phonemeMapping);
		
		// Output
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getName());
			out = new PrintStream(outputInput.get());
		}
		for (String phoneme : phonemeMapping.keySet()) {
			out.println(phoneme + "\t" + phonemeMapping.get(phoneme));
		}
		out.close();
		
		
		List<String> uniqueMapped = new ArrayList<>();
		for (String phoneme : phonemeMapping.keySet()) {
			String value =  phonemeMapping.get(phoneme);
			if (!uniqueMapped.contains(value)) uniqueMapped.add(value);
		}
		
		
		Log.warning(uniqueMapped.size() + " phonemes down from " + phonemes.size() + ". Values: " + uniqueMapped);
		Log.warning("Detected " + errorCount + " errors");
		if (errorCount > 0) {
			Log.warning("Unable to map the following" + unMapped);
		}
		Log.warning("\n");
		
		
	}
	
	
	
	
	public static boolean isIPA(String phoneme, boolean showLongSymbolWarning) {
		
		String phoneme2 = phoneme.trim();
		
		for (String p : ConvertToIPA.IPAs) {
			
			
			p = p.trim();
			
			if (p.equals(phoneme2)) return true;
			
			// Allow a long-sound (colons or diamonds) after the phoneme
			if ((p+":").equals(phoneme2)) {
				if (showLongSymbolWarning) {
					Log.warning("'" + phoneme + "' is considered a long-sound (colon ':')");
				}
				return true;
			}
			if ((p+"ː").equals(phoneme2)) {
				if (showLongSymbolWarning) {
					Log.warning("'" + phoneme + "' is considered a long-sound (diamonds 'ː')");
				}
				return true;
			}
			
			
			if (p.equals(phoneme2.toLowerCase())) {
				if (showLongSymbolWarning) {
					Log.warning("Warning: '" + phoneme + "' is uppercase. Express in lower case form '" + phoneme2.toLowerCase() + "' to match to a phoneme");
				}
			}
			
		}
		return false;
	}
	
	
	/**
	 * Sort a table by count
	 * @param table
	 * @return
	 */
	private static HashMap<String, Integer> sortTableByCount(HashMap<String, Integer> table) {
		
		 // Create a list from elements of HashMap
        List<Map.Entry<String, Integer> > list = new LinkedList<Map.Entry<String, Integer>>(table.entrySet());
 
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
	
	
	
	/**
	 * Sort a table by string value
	 * @param table
	 * @return
	 */
	private static HashMap<String, String> sortTableByValue(HashMap<String, String> table) {
		
		 // Create a list from elements of HashMap
        List<Map.Entry<String, String> > list = new LinkedList<Map.Entry<String, String>>(table.entrySet());
 
        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<String, String> >() {
            public int compare(Map.Entry<String, String> o1,
                               Map.Entry<String, String> o2)
            {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });
         
        // put data from sorted list to hashmap
        HashMap<String, String> temp = new LinkedHashMap<String, String>();
        for (Map.Entry<String, String> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
		
		
	}
	
	
	

	public static void main(String[] args) throws Exception {
		new Application(new ConvertToIPA(), "Attempt to convert all phonemes into IPA", args);
	}

}

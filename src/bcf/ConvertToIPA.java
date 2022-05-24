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
	public Input<Boolean> splitDipthongsInput = new Input<>("dipthong", "set true to split dipthongs into 2 phonemes", true);
	public Input<Boolean> removeLongConsonantsInput = new Input<>("removeLongConsonants", "set true to remove the long modifier from consonants", true);
	
	
	
	
	Map<String, String> encodingMapping = new HashMap<>();
	
	// https://www.internationalphoneticassociation.org/IPAcharts/IPA_chart_orig/IPA_charts_E.html
	// 2020
	// I converted some of the 'Additional symbols' into 2 letter forms (eg. ʦ -> ts), and I replaced the non ascii 'ɡ' with the ascii 'g'
	final static String[] CONSONANTS =
			("p b t d ʈ ɖ c ɟ k g q ɢ ʔ m ɱ n ɳ ɲ ŋ ɴ ʙ r ʀ ⱱ ɾ ɽ ɸ β f v θ ð s z ʃ ʒ ʂ ʐ ç ʝ x ɣ χ ʁ ħ ʕ h ɦ ɬ ɮ ʋ ɹ ɻ j ɰ l ɭ ʎ ʟ ɫ " + // Pulmonic
			"ʘ ǀ ǃ ǂ ǁ ɓ ɗ ʄ ɠ ʛ ʼ pʼ tʼ kʼ sʼ " + // Non-pulmonic
			"ʍ w ɥ ʜ ʢ ʡ ɕ ʑ ɺ ɧ " + // Other
			"ts dz tʃ dʒ ʨ ʥ ɝ".toLowerCase()).split(" ") ; // Additional symbols
	
	
	final static String[] VOWELS =
			("i y ɪ ʏ e ø ɛ œ æ a ɶ ɨ ʉ ɘ ɵ ə ɜ ɞ ɐ ɯ u ʊ ɤ o ʌ ɔ ɑ ɒ ɚ".toLowerCase()).split(" ") ; 

	
	
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
		List<String> uniqueMapped = new ArrayList<>();
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
			
			
			boolean isConsonant = isConsonant(phoneme, true);
			boolean isVowel = isVowel(phoneme, true);
			
			if (isConsonant || isVowel) {
				String phonemeTidy = phoneme;
				if (isConsonant) {
					if (this.removeLongConsonantsInput.get()) phonemeTidy = removeLong(phonemeTidy, true);
					Log.warning("Success: '" + phonemeTidy + "' is already an IPA consonant");
				}else {
					Log.warning("Success: '" + phonemeTidy + "' is already an IPA vowel");
				}
				if (!uniqueMapped.contains(phonemeTidy)) uniqueMapped.add(phonemeTidy);
				phonemeMapping.put(phonemeDataset, phonemeTidy);
			}else {
				
				
				
				
					
				if (isDipthong(phoneme)) {
				
					String[] dipthongArr = splitDipthong(phoneme);
					String dipthong = dipthongArr[0] + " " + dipthongArr[1];
					Log.warning("'" + phonemeDataset + "' is a dipthong '" + dipthong + "'");
					
					// Split dipthongs?
					if (splitDipthongsInput.get()) {
						phonemeMapping.put(phoneme, dipthong);
						if (!uniqueMapped.contains(dipthongArr[0])) uniqueMapped.add(dipthongArr[0]);
						if (!uniqueMapped.contains(dipthongArr[1])) uniqueMapped.add(dipthongArr[1]);
						Log.warning("Success: '" + phonemeDataset + "' is being split into dipthong '" + dipthong + "'");
						continue;
					}
					
				}
				
				
				
				// Split by '/' . Typically the second part is IPA and the first part is not
				String[] bits = phoneme.split("/");
				if (bits.length == 2) {
					
					boolean firstIsIPA = isConsonant(bits[0], false) || isVowel(bits[0], false);
					boolean secondIsIPA = isConsonant(bits[1], false) || isVowel(bits[1], false);
					boolean secondIsDipthong = isDipthong(bits[1]);
					
					// Split dipthongs? Just consider the second part
					if (splitDipthongsInput.get() && secondIsDipthong) {
						String[] dipthongArr = splitDipthong(bits[1]);
						String dipthong = dipthongArr[0] + " " + dipthongArr[1];
						phonemeMapping.put(phoneme, dipthong);
						if (!uniqueMapped.contains(dipthongArr[0])) uniqueMapped.add(dipthongArr[0]);
						if (!uniqueMapped.contains(dipthongArr[1])) uniqueMapped.add(dipthongArr[1]);
						Log.warning("Success: '" + phonemeDataset + "' is being split into dipthong '" + dipthong + "'");
						continue;
					}
					
					
					// Map to first 
					if (firstIsIPA && !secondIsIPA) {
						
						String phonemeTidy = bits[0];
						if (this.removeLongConsonantsInput.get() && isConsonant(phonemeTidy, false)) phonemeTidy = removeLong(phonemeTidy, true);
						
						phonemeMapping.put(phoneme, phonemeTidy);
						if (!uniqueMapped.contains(phonemeTidy)) uniqueMapped.add(phonemeTidy);
						Log.warning("Success: '" + phonemeDataset + "' is being mapped to '" + phonemeTidy + "'");
						continue;
					}
					
					// Map to second
					else if (!firstIsIPA && secondIsIPA) {
						
						String phonemeTidy = bits[1];
						if (this.removeLongConsonantsInput.get() && isConsonant(phonemeTidy, false)) phonemeTidy = removeLong(phonemeTidy, true);
						
						phonemeMapping.put(phoneme, phonemeTidy);
						if (!uniqueMapped.contains(phonemeTidy)) uniqueMapped.add(phonemeTidy);
						Log.warning("Success: '" + phonemeDataset + "' is being mapped to '" + phonemeTidy + "'");
						continue;
					}
					
					// Map to most probable. Take second as a tie breaker
					else if (firstIsIPA && secondIsIPA) {
						
						String phonemeTidy1 = bits[0];
						String phonemeTidy2 = bits[1];
						if (this.removeLongConsonantsInput.get() && isConsonant(phonemeTidy1, false)) phonemeTidy1 = removeLong(phonemeTidy1, true);
						if (this.removeLongConsonantsInput.get() && isConsonant(phonemeTidy2, false)) phonemeTidy2 = removeLong(phonemeTidy2, true);
						
						
						int count1 = phonemes.containsKey(phonemeTidy1) ? phonemes.get(phonemeTidy1) : 0;
						int count2 = phonemes.containsKey(phonemeTidy2) ? phonemes.get(phonemeTidy2) : 0;
						String mapped = count1 > count2 ? phonemeTidy1 : phonemeTidy2;
						
						phonemeMapping.put(phoneme, mapped);
						if (!uniqueMapped.contains(mapped)) uniqueMapped.add(mapped);
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
		
		// Print encoding and ensure that everything is 2 characters. Add . if necessary
		for (String phoneme : phonemeMapping.keySet()) {
			String value = phonemeMapping.get(phoneme);
			
			String[] bits = value.split(" ");
			if (bits.length == 1) {
				if (value.length() == 1) value = value + ".";
				if (value.length() > 2) {
					Log.warning("Error: '" + value + "' is too long. All phonemes must be 1-2 characters.");
					errorCount++;
				}
			}else if (bits.length == 2) {
				String value1 = bits[0];
				String value2 = bits[1];
				if (value1.length() == 1) value1 = value1 + ".";
				if (value2.length() == 1) value2 = value2 + ".";
				value = value1 + " " + value2;
			}else {
				Log.warning("Error: '" + value + "' is too long. All phonemes must be 1-2 characters.");
				errorCount++;
			}
			
			out.println(phoneme + "\t" + value);
		}
		out.close();
		

		
		Collections.sort(uniqueMapped);
		
		
		Log.warning(uniqueMapped.size() + " phonemes down from " + phonemes.size() + ". Values: " + uniqueMapped);
		Log.warning("Detected " + errorCount + " errors");
		if (errorCount > 0) {
			Log.warning("Unable to map the following" + unMapped);
		}
		Log.warning("\n");
		
		
	}
	
	
	/**
	 * Remove long modifier from phoneme
	 * @param phoneme
	 * @return
	 */
	public static String removeLong(String phoneme, boolean showWarning) {
		String phoneme2 = phoneme.replaceAll("(:|ː)", "");
		if (showWarning && !phoneme.equals(phoneme2)) {
			Log.warning("Removing long sound from consonant '" + phoneme + "' to give '" + phoneme2 + "'");
		}
		return phoneme2;
	}
	
	
	/**
	 * Check if the phoneme is a consonant
	 * @param phoneme
	 * @param showLongSymbolWarning
	 * @return
	 */
	public static boolean isConsonant(String phoneme, boolean showLongSymbolWarning) {
		return phonemeOnList(phoneme, ConvertToIPA.CONSONANTS, showLongSymbolWarning);
	}
	
	
	/**
	 * Check if the phoneme is a vowel
	 * @param phoneme
	 * @param showLongSymbolWarning
	 * @return
	 */
	public static boolean isVowel(String phoneme, boolean showLongSymbolWarning) {
		return phonemeOnList(phoneme, ConvertToIPA.VOWELS, showLongSymbolWarning);
	}
	
	
	/**
	 * It is a dipthong if it made of 2 vowels. eg. oːy
	 * @param phoneme
	 * @param showLongSymbolWarning
	 * @return
	 */
	public static boolean isDipthong(String phoneme) {
		return splitDipthong(phoneme) != null;
	}
	
	
	/**
	 * Attempt to split a string into 2 vowels, or return null if not possible. eg. oːy
	 * @param phoneme
	 * @return
	 */
	public static String[] splitDipthong(String phoneme) {
		
		// Split by diamond or :
		//String[] bits = phoneme.split("(:|ː)");
		
		
		// Can contain 2,3, or 4 characters
		if (phoneme.length() < 2 || phoneme.length() > 4) return null;
		
		// Remove long sounds
		String phoneme2 = phoneme.replaceAll("(:|ː)", "");
		
		// Must be made of 2 vowels
		if (phoneme2.length() != 2) return null;
		String part1 = phoneme2.substring(0, 1);
		String part2 = phoneme2.substring(1, 2);
		if (!isVowel(part1, false) || !isVowel(part2, false)) return null;
		
		
		// This is a dipthong. Add back the long sounds
		String part1b = phoneme.substring(1,2);
		boolean part1IsLong = false;
		if (part1b.equals(":") || part1b.equals("ː")) {
			part1IsLong = true;
			part1 = part1 + part1b;
		}
		
		
		String part2b = "";
		if (part1IsLong && phoneme.length() == 4) {
			part2b = phoneme.substring(3,4);
		}
		else if (!part1IsLong && phoneme.length() == 3) {
			part2b = phoneme.substring(2,3);
		}
		if (part2b.equals(":") || part2b.equals("ː")) part2 = part2 + part2b;
		
		
		
		
		String[] parts = new String[] {part1, part2};
		return parts;
		
	}
	
	
	
	private static boolean phonemeOnList(String phoneme, String[] arr, boolean showLongSymbolWarning) {
		
		String phoneme2 = phoneme.trim();
		for (String p : arr) {
			
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

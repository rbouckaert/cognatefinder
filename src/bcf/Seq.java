package bcf;

import java.util.Arrays;

public class Seq {
	String [] characters;
	int [] codes;
	
	Seq(String [] chars, Score score) {
		characters = chars;
		codes = new int[chars.length];
		for (int i = 0; i < chars.length; i++) {
			codes[i] = score.getCode(chars[i]);
		}
	}
	
	@Override
	public String toString() {
		return Arrays.toString(characters);
	}

	public void makeLenght(int len) {
		if (characters.length < len) {
			String [] newChars = new String[len];
			for (int i = 0; i < characters.length - 1; i++) {
				newChars[i] = characters[i];
			}
			for (int i = characters.length - 1; i < len -1; i++) {
				newChars[i] = "-";
			}
			newChars[len - 1] = characters[characters.length - 1];
			characters = newChars;
			// todo: update Seq.codes
		}
		
	}

}

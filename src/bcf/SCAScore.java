package bcf;

import java.util.HashMap;
import java.util.Map;

public class SCAScore extends Score {
	double [][] score = new double[][] { 
			/* + */ {0.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -100.00 , -5.00 , -100.00 , -100.00},
			/* 0 */ {-100.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00},
			/* 1 */ {-100.00 , 0.00 , 2.00 , 1.00 , 1.00 , 1.00 , 1.00 , 1.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 2 */ {-100.00 , 0.00 , 1.00 , 2.00 , 1.00 , 1.00 , 1.00 , 1.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 3 */ {-100.00 , 0.00 , 1.00 , 1.00 , 2.00 , 1.00 , 1.00 , 1.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 4 */ {-100.00 , 0.00 , 1.00 , 1.00 , 1.00 , 2.00 , 1.00 , 1.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 5 */ {-100.00 , 0.00 , 1.00 , 1.00 , 1.00 , 1.00 , 2.00 , 1.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 6 */ {-100.00 , 0.00 , 1.00 , 1.00 , 1.00 , 1.00 , 1.00 , 2.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -20.00 , -20.00},
			/* 9 */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 10.00 , -10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 9.00 , 9.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* A */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 5.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , 4.00 , -6.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -6.00 , 0.00 , 4.00 , -10.00},
			/* B */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 6.00 , 0.00 , 0.00 , 0.00 , -10.00 , 6.00 , 0.00 , -9.00 , -10.00},
			/* C */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 10.00 , 2.00 , -10.00 , 2.00 , 2.00 , -10.00 , 0.00 , 6.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 6.00 , 6.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* D */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 2.00 , 10.00 , -10.00 , 0.00 , 2.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 6.00 , 6.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* E */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , 5.00 , -10.00 , -10.00 , 4.00 , -6.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -6.00 , 0.00 , 4.00 , -10.00},
			/* G */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 2.00 , 0.00 , -10.00 , 10.00 , 2.00 , -10.00 , 0.00 , 6.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 6.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* H */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 2.00 , 2.00 , -10.00 , 2.00 , 10.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 6.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* I */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , 5.00 , -5.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , -10.00 , 3.00 , -6.00 , 0.00 , 4.00 , -10.00},
			/* J */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -6.00 , 0.00 , 0.00 , 0.00 , -6.00 , 0.00 , 0.00 , -5.00 , 10.00 , 0.00 , 0.00 , 0.00 , 0.00 , -6.00 , 0.00 , 0.00 , 0.00 , 0.00 , -6.00 , 0.00 , 0.00 , -6.00 , -10.00},
			/* K */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 6.00 , 0.00 , -10.00 , 6.00 , 0.00 , -10.00 , 0.00 , 10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 2.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* L */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 4.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* M */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 9.00 , -10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 10.00 , 1.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* N */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 9.00 , -10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 1.00 , 10.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* O */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , 4.00 , -6.00 , -10.00 , -10.00 , -10.00 , -10.00 , 5.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -6.00 , 0.00 , 4.00 , -10.00},
			/* P */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 6.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 2.00 , 0.00 , -10.00 , -10.00},
			/* R */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 4.00 , 0.00 , 0.00 , -10.00 , 0.00 , 10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* S */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 6.00 , 6.00 , -10.00 , 6.00 , 6.00 , -10.00 , 0.00 , 2.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 10.00 , 2.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* T */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -10.00 , 0.00 , 6.00 , 6.00 , -10.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -10.00 , 0.00 , 0.00 , 2.00 , 10.00 , -10.00 , 0.00 , 0.00 , -10.00 , -10.00},
			/* U */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , 3.00 , -6.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , -10.00 , 5.00 , -6.00 , 0.00 , 4.00 , -10.00},
			/* W */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , 0.00 , -6.00 , 6.00 , 0.00 , 0.00 , -6.00 , 0.00 , 0.00 , -6.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , -6.00 , 2.00 , 0.00 , 0.00 , 0.00 , -6.00 , 10.00 , 0.00 , -5.00 , -10.00},
			/* X */ {-5.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00},
			/* Y */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , 4.00 , -9.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , 4.00 , -6.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -10.00 , -10.00 , -10.00 , -10.00 , 4.00 , -5.00 , 0.00 , 5.00 , -10.00},
			/* _ */ {-100.00 , 0.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -20.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , -10.00 , 0.00 , -10.00 , 10.00}
			};

			
	Map<String,Integer> charMap;
	public SCAScore() {
		charMap = new HashMap<>();		
		/* 1 */ process(2, "₁₁,₂₂,¹¹,²²");
		/* 2 */ process(3, "₁₂,₁₃,₁₄,₁₅,₂₃,₂₄,₂₅,₃₄,₃₅,₄₅,¹²,¹³,¹⁴,¹⁵,²³,²⁴,²⁵,³⁵,³⁴,⁴⁵");
		/* 3 */ process(4, "₅₁,₅₂,₅₃,₅₄,₄₃,₄₂,₄₁,₃₂,₃₁,₂₁,⁵¹,⁵²,⁵³,⁵⁴,⁴¹,⁴²,⁴³,³¹,³²,²¹");
		/* 4 */ process(5, "₃₃,³³");
		/* 5 */ process(6, "₄₄,₅₅,⁵⁵,⁴⁴");
		/* 6 */ process(7, "⁰,¹,²,³,⁴,⁵,⁻,₁,₂,₃,₄,˥,˦,˧,˨,₅,₆,₀");
		/* 9 */ process(8, "∼");
		/* A */ process(9, "a,ᴀ,ã,ɑ,á,á,à,ā,ǎ,â,ä,ă,ă,ạ,а,å,ạ");
		/* B */ process(10, "ɸ,β,f,p͡f,p͜f,ƀ,p͡f,b͡v,pf,bv,v,ʙ,ḇ");
		/* C */ process(11, "t͡s,t͜s,d͡z,d͜z,ʦ,ʣ,t͡ɕ,t͜ɕ,d͡ʑ,d͜ʑ,ʨ,ʥ,t͡ʃ,t͜ʃ,d͡ʒ,d͜ʒ,ʧ,ʤ,c,ɟ,t͡ʂ,t͜ʂ,d͡ʐ,d͜ʐ,č,t͡θ,t͜θ,ʄ,ǰ,ĵ,Ɉ,ʈʂ,ɖʐ,ʈʂʰ,tɕ,tɕʰ,dʑ,ts,dz,tsʰ");
		/* D */ process(12, "θ,ð,ŧ,þ,đ,Ɵ");
		/* E */ process(13, "ɛ,æ,ɜ,ɐ,ʌ,e,ᴇ,ə,ɘ,ɤ,è,é,ē,ě,ê,ɚ,ǝ,ẽ,ĕ,ḛ,ε,е,ę,ḛ,ȇ,ë,ε");
		/* G */ process(14, "x,ɣ,χ");
		/* H */ process(15, "ʔ,ħ,ʕ,h,ɦ,ḥ,Ɂ,ʡ,'ʷ");
		/* I */ process(16, "i,ɪ,ɨ,ɿ,ʅ,ɯ,ĩ,í,ǐ,ì,î,ī,ı,ĭ,ḭ,ɩ,ï,ị,ๅ,ḭ");
		/* J */ process(17, "j,ɥ,ɰ");
		/* K */ process(18, "k,g,q,ɢ,ɡ,ḳ,ǥ,ǵ,ḡ");
		/* L */ process(19, "l,ȴ,l,ɭ,ʎ,ʟ,ɬ,ɮ,ł,ɫ,ḷ,ļ");
		/* M */ process(20, "m,ɱ,ʍ,ṃ");
		/* N */ process(21, "n,ȵ,ɳ,ŋ,ɴ,ň,ń,ɲ,ñ,ṇ,ῃ,ņ,ṋ");
		/* O */ process(22, "Œ,ɒ");
		/* P */ process(23, "p,b,ɓ,ᵐb,ᵐp,р");
		/* R */ process(24, "ɹ,ɻ,ʀ,ɾ,r,ʁ,ɽ,ɐ̯,ɺ,ṛ,ᵲ,ř,ȓ,ṛ́,ṙ");
		/* S */ process(25, "s,z,ʃ,∫,ʒ,ʂ,ʐ,ç,ʝ,š,ž,ɕ,ɧ,ʑ,ś,ṣ,ß,ŝ,ż,ẓ");
		/* T */ process(26, "t,d,ȶ,ȡ,ɗ,ʈ,ɖ,ţ,т,ṱ,ṭ,ḍ,ḏ");
		/* U */ process(27, "œ,ɞ,ɔ,ø,ɵ,o,õ,ó,ò,ō,ɶ,ô,ɷ,ǒ,ö,ŏ,ʮ,ọ,ȯ,ố,ǫ,ṍ");
		/* W */ process(28, "w,ʋ,ⱱ,ṿ,υ,ṽ");
		/* Y */ process(30, "y,ỹ,ỹ,ṳ,ṵ,ʏ,ʉ,u,ᴜ,ʊ,ú,ù,ũ,ü,ŭ,ǔ,ụ,ū,ỳ,û,û,ý,ў,ȗ,ṹ");
		/* _ */ process(31, "_,#,+,◦,·");
	}
	
	private void process(int i, String string) {
		for (String s : string.split(",")) {
			charMap.put(s,  i);
		}		
	}

	@Override
	float score(String a, String b) {
		int i = charMap.get(a);
		int j = charMap.get(b);
		return (float) score[i][j];
	}

}

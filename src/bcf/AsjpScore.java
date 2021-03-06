package bcf;

import beast.core.Description;

@Description("ASJP Sound-Class Model "
		+ "Sound-Class model following Brown et al. (2008) and Brown et al. (2011) "
		+ "Sourced from lingpy/data/models/asjp.")
public class AsjpScore extends Score {

	public AsjpScore() {	
		/* ! */ process(47," ǃ, ǂ, ǁ, ǀ, ʘ, gǃ, gǂ, gǁ, gǀ, gʘ, ǃŋ, ǂŋ, ǁŋ, ǀŋ, ʘŋ, !, |, g!, g|, !ŋ, |ŋ");
		/* p */ process(36," p, ᵐp, р, p͡f, pf, p͜f");
		/* b */ process(22," b, ɓ, ᵐb, β, b͡v, ḇ, bv, ʙ");
		/* f */ process(26," f, ɸ");
		/* v */ process(42," v, ƀ, ⱱ, ṿ, ṽ");
		/* m */ process(33," m, ɱ, ṃ");
		/* w */ process(43," w, ɰ, ʋ, υ, ʍ");
		/* 8 */ process(8," θ, Ɵ, ð, đ, þ");
		/* t */ process(40," t, ʈ, ţ, т, ṱ, ṭ, ŧ, ȶ, ƫ");
		/* T */ process(17," c, Ɉ");
		/* d */ process(24," d, ɖ, ḍ, ḏ, ȡ, ɗ, ᶁ");
		/* s */ process(39," s, ś, ṣ, ß, ᶊ");
		/* z */ process(46," z, ż, ẓ, ᶎ");
		/* c */ process(23," ʦ, ʣ, t͡s, d͡z");
		/* n */ process(34," n, ɳ, ṋ, ṇ, ñ, ῃ, ņ, ň, ń, ∼");
		/* r */ process(38," ɾ, r, ʀ, ɽ, ɐ̯, ɺ, ɹ, ṛ, ř, ȓ, ṛ́, ṙ, ᵲ, ᶉ");
		/* l */ process(32," l, ʎ, ł, ɫ, ḷ, ļ, ᶅ");
		/* S */ process(16," ʃ, ∫, ʂ, ɕ, ç, ɧ, ʝ, š, ʆ, ŝ");
		/* Z */ process(19," ʒ, ʑ, ʐ, ž");
		/* C */ process(9," ʧ, ʨ, ɟ, t͡ʃ, t͡ɕ, t͡ʂ, t͡θ, t͜θ, t͜ʃ, d͜ʐ, ʄ, č, t͜ɕ, t͜ʂ, ʈʂ, ʈʂʰ, tɕ, tɕʰ, ts, tsʰ");
		/* j */ process(30," ʤ, ʥ, d͡ʒ, d͡ʑ, d͡ʐ, ǰ, ĵ, d͜z, d͜ʒ, t͜s, d͜ʑ, ɖʐ, dz, dʑ");
		/* H */ process(12," ɲ, ȵ, ᶇ");
		/* y */ process(45," j, ɥ");
		/* k */ process(31," k, ḳ");
		/* g */ process(27," g, ɡ, ḡ, ǵ, ɠ");
		/* x */ process(44," x, ɣ, ǥ");
		/* N */ process(15," ŋ, ɴ");
		/* q */ process(37," q, ɢ, ʛ");
		/* G */ process(11," χ, ʁ, ħ, ʕ, ʡ");
		/* h */ process(28," h, ḥ, ɦ");
		/* 7 */ process(7," Ɂ, ʔ, 'ʷ");
		/* L */ process(14," ʟ, ɭ, ȴ, ɭ, ʟ, ɬ, ɮ, ɻ");
		/* i */ process(29," i, y, ỹ, ṳ, ʏ, ɪ, ı, ɩ, ɿ, ʅ, ï, ĩ, í, ǐ, ì, î, ī, ü, ĭ, ỳ, ý, ў, ḭ, ๅ, ị, ŷ, ʯ");
		/* e */ process(25," e, ø, ẽ, è, é, ē, ě, ê, ĕ, е, ḛ, è, é, ê, ë, ē, ĕ, ę, ě, ȇ");
		/* E */ process(10," æ, ɛ, œ, Œ, ᴇ, ɶ, ö, ε, ε");
		/* I */ process(13," ɨ, ǝ, ɘ, ə, ɚ, ɜ, ɵ, ʉ");
		/* a */ process(21," a, à, á, â, ã, ä, å, ā, ă, ǎ, ɐ, а, ᴀ, ạ, ạ");
		/* u */ process(41," u, ɯ, ʊ, ɷ, ʮ, ᴜ, ɞ, ŭ, ú, ṵ, ũ, ǔ, ụ, ū, û, ù, ȗ, ṹ");
		/* o */ process(35," ɤ, ʌ, ɑ, o, ɔ, ɒ, õ, ŏ, ŏ, ọ, ȯ, ǫ, ō, ò, ô, ố, ǒ, ó, ṍ");
		/* + */ process(0," +");
		/* 1 */ process(1," ₁₁, ₂₂, ¹¹, ²²");
		/* 2 */ process(2," ₁₂, ₁₃, ₁₄, ₁₅, ₂₃, ₂₄, ₂₅, ₃₄, ₃₅, ₄₅, ¹², ¹³, ¹⁴, ¹⁵, ²³, ²⁴, ²⁵, ³⁵, ³⁴, ⁴⁵");
		/* 3 */ process(3," ₅₁, ₅₂, ₅₃, ₅₄, ₄₃, ₄₂, ₄₁, ₃₂, ₃₁, ₂₁, ⁵¹, ⁵², ⁵³, ⁵⁴, ⁴¹, ⁴², ⁴³, ³¹, ³², ²¹");
		/* 4 */ process(4," ₃₃, ³³");
		/* 5 */ process(5," ₄₄, ₅₅, ⁵⁵, ⁴⁴");
		/* 6 */ process(6," ⁰, ¹, ², ³, ⁴, ⁵, ⁻, ₁, ₂, ₃, ₄, ₅, ˥, ˦, ˧, ˨, ₆, ₀");
		/* _ */ process(20," _, ◦, #, ·");
	setScore(new double[][] {
		/* +	*/ {0.00, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -2.19, -8.75, 0.00, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 1	*/ {-8.75, 4.38, 2.19, 2.19, 2.19, 2.19, 2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 2	*/ {-8.75, 2.19, 4.38, 2.19, 2.19, 2.19, 2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 3	*/ {-8.75, 2.19, 2.19, 4.38, 2.19, 2.19, 2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 4	*/ {-8.75, 2.19, 2.19, 2.19, 4.38, 2.19, 2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 5	*/ {-8.75, 2.19, 2.19, 2.19, 2.19, 4.38, 2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 6	*/ {-8.75, 2.19, 2.19, 2.19, 2.19, 2.19, 4.38, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -21.88, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75},
		/* 7	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 10.00, 0.00, 0.00, -8.31, 1.31, 0.00, -8.31, 0.88, 0.44, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.44, 0.00, 0.44, -8.31, 0.00, 1.31, 2.62, -8.31, 0.00, 3.94, 0.88, 0.00, 0.44, -8.31, 0.00, 3.50, 0.88, 0.44, 1.75, -8.31, 0.00, 0.88, 1.75, 0.88, 0.00, -8.75},
		/* 8	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 10.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 1.31, 1.31, 0.00, 0.00, -21.88, -8.31, 0.00, 1.31, 4.38, -8.31, 0.00, 0.00, 0.88, -8.31, 0.00, 0.00, 1.75, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, 3.06, 4.81, -8.31, 0.00, 0.00, 0.00, 1.75, 1.31, -8.75},
		/* C	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 10.00, -8.31, 1.31, 0.00, -8.31, 0.00, 0.00, 1.31, 8.31, 0.00, 0.88, -21.88, -8.31, 0.00, 7.44, 0.88, -8.31, 0.00, 0.00, 1.31, -8.31, 0.44, 0.88, 0.44, 0.00, 0.44, -8.31, 0.00, 0.00, 0.88, 2.62, 2.62, -8.31, 0.00, 0.00, 0.44, 0.44, 0.44, 0.00},
		/* E	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 8.75, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 7.88, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, -8.31, -8.31, 5.69, -8.31, -8.31, -8.31, -8.31, -8.31, 5.25, -8.31, -8.31, -8.31, -8.31, -8.31, -8.75},
		/* G	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 1.31, 0.00, 1.31, -8.31, 10.00, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.00, 0.00, 0.00, -8.31, 0.00, 1.31, 3.50, -8.31, 0.00, 3.94, 0.00, 0.00, 0.00, -8.31, 0.00, 2.62, 6.12, 0.88, 0.00, -8.31, 0.00, 0.00, 9.19, 0.00, 0.00, 0.00},
		/* H	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, -8.31, 0.00, 10.00, -8.31, 0.00, 2.19, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, 0.44, 0.44, 3.94, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, 0.44, 0.44, 0.00, 0.00},
		/* I	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 6.56, -8.31, -8.31, 8.75, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 7.00, -8.31, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, 7.00, -8.31, -8.31, -8.31, -8.31, -8.31, 6.12, -8.31, -8.31, -8.31, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, -8.31, -8.31, -8.75},
		/* L	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.88, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 10.00, 0.00, 0.00, 1.31, 0.00, 0.00, -21.88, -8.31, 0.00, 0.88, 0.88, -8.31, 0.00, 0.00, 1.31, -8.31, 0.00, 0.00, 8.31, 0.00, 0.44, -8.31, 0.00, 0.00, 3.50, 0.00, 1.31, -8.31, 0.00, 0.44, 0.88, 1.31, 1.31, -8.75},
		/* N	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.44, 0.00, 0.00, -8.31, 0.00, 2.19, -8.31, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.00, 0.00, 0.44, -8.31, 0.00, 0.88, 0.44, -8.31, 0.00, 0.88, 0.44, 0.88, 5.69, -8.31, 0.00, 0.00, 0.44, 0.00, 0.44, -8.31, 0.00, 0.00, 0.88, 0.88, 0.00, 0.00},
		/* S	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 1.31, 1.31, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 10.00, 0.88, 0.00, 0.88, -21.88, -8.31, 0.00, 0.44, 0.00, -8.31, 0.00, 0.00, 1.31, -8.31, 0.44, 0.44, 0.00, 0.00, 0.44, -8.31, 0.00, 0.00, 0.44, 7.88, 0.44, -8.31, 0.00, 0.00, 0.44, 0.44, 1.31, -8.75},
		/* T	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 1.31, 8.31, -8.31, 0.00, 0.00, -8.31, 1.31, 0.00, 0.88, 10.00, 0.00, 3.94, -21.88, -8.31, 0.00, 2.62, 2.19, -8.31, 0.00, 0.88, 0.00, -8.31, 4.38, 1.31, 0.44, 0.00, 0.00, -8.31, 0.00, 0.00, 0.88, 1.31, 1.31, -8.31, 0.00, 0.44, 0.00, 1.31, 1.75, 0.00},
		/* X	*/ {-2.19, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
		/* Z	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.88, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.88, 3.94, 0.00, 10.00, -21.88, -8.31, 0.00, 0.88, 0.44, -8.31, 0.00, 0.44, 0.00, -8.31, 7.00, 0.00, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, 1.75, 0.44, 0.00, -8.31, 0.00, 0.00, 0.00, 3.06, 1.75, -8.75},
		/* _	*/ {0.00, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, 0.00, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -21.88, -8.75},
		/* a	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, 7.00, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 8.75, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, -8.31, 6.12, -8.31, -8.31, -8.31, -8.31, -8.31, 7.44, -8.31, -8.31, -8.31, -8.31, -8.31, 5.25, -8.31, -8.31, -8.31, -8.31, -8.31, -8.75},
		/* b	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.44, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 10.00, 0.00, 0.44, -8.31, 1.75, 0.44, 0.44, -8.31, 0.00, 0.00, 0.44, 0.88, 0.00, -8.31, 4.38, 0.00, 0.00, 0.00, 0.00, -8.31, 6.12, 2.62, 0.44, 0.00, 0.00, -8.75},
		/* c	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 1.31, 7.44, -8.31, 0.00, 0.00, -8.31, 0.88, 0.00, 0.44, 2.62, 0.00, 0.88, -21.88, -8.31, 0.00, 10.00, 1.31, -8.31, 0.00, 0.88, 0.44, -8.31, 3.50, 0.44, 0.88, 0.44, 0.00, -8.31, 0.00, 0.00, 0.44, 4.81, 2.19, -8.31, 0.00, 0.00, 0.00, 0.88, 1.31, 0.00},
		/* d	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.44, 4.38, 0.88, -8.31, 0.00, 0.00, -8.31, 0.88, 0.44, 0.00, 2.19, 0.00, 0.44, -21.88, -8.31, 0.44, 1.31, 10.00, -8.31, 0.44, 0.88, 0.44, -8.31, 1.75, 0.44, 1.75, 0.00, 1.75, -8.31, 0.00, 0.00, 4.38, 0.44, 5.25, -8.31, 0.44, 0.00, 0.44, 0.44, 1.75, 0.00},
		/* e	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 7.88, -8.31, -8.31, -8.31, 8.75, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, -8.31, -8.31, -8.31, 5.69, -8.31, -8.31, -8.31, -8.31, -8.31, 4.81, -8.31, -8.31, -8.31, -8.31, -8.31, -8.75},
		/* f	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 1.75, 0.00, 0.44, -8.31, 10.00, 0.00, 2.19, -8.31, 0.00, 0.44, 0.00, 0.00, 0.00, -8.31, 4.38, 0.00, 0.00, 0.44, 0.44, -2.62, 4.81, 1.31, 0.00, 0.00, 0.00, -8.75},
		/* g	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 1.31, 0.00, 0.00, -8.31, 1.31, 0.00, -8.31, 0.00, 0.88, 0.00, 0.88, 0.00, 0.44, -21.88, -8.31, 0.44, 0.88, 0.88, -8.31, 0.00, 10.00, 1.31, -8.31, 0.44, 5.69, 1.31, 0.00, 0.00, -8.31, 0.00, 0.88, 0.88, 0.00, 0.00, -8.31, 0.44, 0.00, 2.19, 0.44, 0.44, 0.00},
		/* h	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 2.62, 0.88, 1.31, -8.31, 3.50, 0.00, -8.31, 1.31, 0.44, 1.31, 0.00, 0.00, 0.00, -21.88, -8.31, 0.44, 0.44, 0.44, -8.31, 2.19, 1.31, 10.00, -8.31, 0.44, 1.75, 0.88, 0.00, 0.00, -8.31, 1.75, 0.88, 1.31, 3.06, 0.88, -8.31, 0.44, 0.88, 7.00, 0.88, 0.44, 0.00},
		/* i	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 6.56, -8.31, -8.31, 7.00, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 6.12, -8.31, -8.31, -8.31, 7.88, -8.31, -8.31, -8.31, 8.75, -8.31, -8.31, -8.31, -8.31, -8.31, 5.25, -8.31, -8.31, -8.31, -8.31, -8.31, 6.12, -8.31, -8.31, -8.31, -0.88, -8.31, -8.75},
		/* j	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.44, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.44, 4.38, 0.00, 7.00, -21.88, -8.31, 0.00, 3.50, 1.75, -8.31, 0.00, 0.44, 0.44, -8.31, 10.00, 0.00, 1.75, 0.44, 0.44, -8.31, 0.00, 0.00, 0.44, 0.00, 1.31, -8.31, 0.00, 0.44, 0.44, 5.25, 3.50, 0.00},
		/* k	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 3.94, 0.00, 0.88, -8.31, 3.94, 0.00, -8.31, 0.00, 0.88, 0.44, 1.31, 0.00, 0.00, -21.88, -8.31, 0.00, 0.44, 0.44, -8.31, 0.44, 5.69, 1.75, -8.31, 0.00, 10.00, 0.00, 0.00, 0.00, -8.31, 0.00, 5.25, 0.00, 0.44, 0.88, -8.31, 0.00, 0.44, 3.06, 0.00, 0.00, 0.00},
		/* l	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.88, 1.75, 0.44, -8.31, 0.00, 0.44, -8.31, 8.31, 0.44, 0.00, 0.44, 0.00, 0.00, -21.88, -8.31, 0.44, 0.88, 1.75, -8.31, 0.00, 1.31, 0.88, -8.31, 1.75, 0.00, 10.00, 0.00, 1.31, -8.31, 0.00, 0.00, 7.88, 0.44, 0.88, -8.31, 0.44, 0.00, 1.31, 1.75, 0.88, -8.75},
		/* m	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, -8.31, 0.00, 0.44, -8.31, 0.00, 0.88, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.88, 0.44, 0.00, -8.31, 0.00, 0.00, 0.00, -8.31, 0.44, 0.00, 0.00, 10.00, 1.31, -8.31, 0.44, 0.00, 0.00, 0.00, 0.00, -8.31, 0.44, 0.44, 0.00, 0.00, 0.00, -8.75},
		/* n	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.44, 0.00, 0.44, -8.31, 0.00, 3.94, -8.31, 0.44, 5.69, 0.44, 0.00, 0.00, 0.00, -21.88, -8.31, 0.00, 0.00, 1.75, -8.31, 0.00, 0.00, 0.00, -8.31, 0.44, 0.00, 1.31, 1.31, 10.00, -8.31, 0.00, 0.44, 1.31, 0.44, 0.44, -8.31, 0.00, 0.00, 0.44, 0.44, 0.00, -8.75},
		/* o	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 5.69, -8.31, -8.31, 6.12, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 7.44, -8.31, -8.31, -8.31, 5.69, -8.31, -8.31, -8.31, 5.25, -8.31, -8.31, -8.31, -8.31, -8.31, 8.75, -8.31, -8.31, -8.31, -8.31, -8.31, 7.44, -8.31, -8.31, -8.31, -8.31, -8.31, -8.75},
		/* p	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 4.38, 0.00, 0.00, -8.31, 4.38, 0.00, 1.75, -8.31, 0.00, 0.00, 0.00, 0.44, 0.00, -8.31, 10.00, 0.00, 0.00, 0.00, 0.00, -8.31, 1.75, 0.88, 0.00, 0.00, 0.00, -8.75},
		/* q	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 3.50, 0.00, 0.00, -8.31, 2.62, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 0.00, 0.00, 0.00, -8.31, 0.00, 0.88, 0.88, -8.31, 0.00, 5.25, 0.00, 0.00, 0.44, -8.31, 0.00, 10.00, 0.00, 0.00, 0.44, -8.31, 0.00, 0.00, 4.81, 0.00, 0.00, 0.00},
		/* r	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.88, 0.00, 0.88, -8.31, 6.12, 0.00, -8.31, 3.50, 0.44, 0.44, 0.88, 0.00, 1.75, -21.88, -8.31, 0.00, 0.44, 4.38, -8.31, 0.00, 0.88, 1.31, -8.31, 0.44, 0.00, 7.88, 0.00, 1.31, -8.31, 0.00, 0.00, 10.00, 0.88, 1.31, -8.31, 0.00, 0.44, 2.62, 1.75, 3.06, -8.75},
		/* s	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.44, 3.06, 2.62, -8.31, 0.88, 0.00, -8.31, 0.00, 0.00, 7.88, 1.31, 0.00, 0.44, -21.88, -8.31, 0.00, 4.81, 0.44, -8.31, 0.44, 0.00, 3.06, -8.31, 0.00, 0.44, 0.44, 0.00, 0.44, -8.31, 0.00, 0.00, 0.88, 10.00, 1.31, -8.31, 0.00, 0.00, 0.88, 0.44, 4.38, -8.75},
		/* t	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 1.75, 4.81, 2.62, -8.31, 0.00, 0.00, -8.31, 1.31, 0.44, 0.44, 1.31, 0.00, 0.00, -21.88, -8.31, 0.00, 2.19, 5.25, -8.31, 0.44, 0.00, 0.88, -8.31, 1.31, 0.88, 0.88, 0.00, 0.44, -8.31, 0.00, 0.44, 1.31, 1.31, 10.00, -8.31, 0.00, 0.00, 0.44, 0.00, 0.00, 0.00},
		/* u	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.31, -8.31, -8.31, 5.25, -8.31, -8.31, 6.56, -8.31, -8.31, -8.31, -8.31, 0.00, -8.31, -21.88, 5.25, -8.31, -8.31, -8.31, 4.81, -2.62, -8.31, -8.31, 6.12, -8.31, -8.31, -8.31, -8.31, -8.31, 7.44, -8.31, -8.31, -8.31, -8.31, -8.31, 8.75, -1.75, -0.88, -8.31, -8.31, -8.31, -8.75},
		/* v	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -21.88, -8.31, 6.12, 0.00, 0.44, -8.31, 4.81, 0.44, 0.44, -8.31, 0.00, 0.00, 0.44, 0.44, 0.00, -8.31, 1.75, 0.00, 0.00, 0.00, 0.00, -1.75, 10.00, 7.00, 1.75, 0.00, 0.00, -8.75},
		/* w	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.88, 0.00, 0.00, -8.31, 0.00, 0.00, -8.31, 0.44, 0.00, 0.00, 0.44, 0.00, 0.00, -21.88, -8.31, 2.62, 0.00, 0.00, -8.31, 1.31, 0.00, 0.88, -8.31, 0.44, 0.44, 0.00, 0.44, 0.00, -8.31, 0.88, 0.00, 0.44, 0.00, 0.00, -0.88, 7.00, 10.00, 0.44, 0.00, 0.44, -8.75},
		/* x	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 1.75, 0.00, 0.44, -8.31, 9.19, 0.44, -8.31, 0.88, 0.88, 0.44, 0.00, 0.00, 0.00, -21.88, -8.31, 0.44, 0.00, 0.44, -8.31, 0.00, 2.19, 7.00, -8.31, 0.44, 3.06, 1.31, 0.00, 0.44, -8.31, 0.00, 4.81, 2.62, 0.88, 0.44, -8.31, 1.75, 0.44, 10.00, 0.88, 0.00, 0.00},
		/* y	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.88, 1.75, 0.44, -8.31, 0.00, 0.44, -8.31, 1.31, 0.88, 0.44, 1.31, 0.00, 3.06, -21.88, -8.31, 0.00, 0.88, 0.44, -8.31, 0.00, 0.44, 0.88, -0.88, 5.25, 0.00, 1.75, 0.00, 0.44, -8.31, 0.00, 0.00, 1.75, 0.44, 0.00, -8.31, 0.00, 0.00, 0.88, 10.00, 1.75, -8.75},
		/* z	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, 1.31, 0.44, -8.31, 0.00, 0.00, -8.31, 1.31, 0.00, 1.31, 1.75, 0.00, 1.75, -21.88, -8.31, 0.00, 1.31, 1.75, -8.31, 0.00, 0.44, 0.44, -8.31, 3.50, 0.00, 0.88, 0.00, 0.00, -8.31, 0.00, 0.00, 3.06, 4.38, 0.00, -8.31, 0.00, 0.44, 0.00, 1.75, 10.00, -8.75},
		/* !	*/ {-8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, -8.75, 0.00, 0.00, -8.75, -8.75, 0.00, -8.75, 0.00, 0.00, -8.75, -8.75, -8.75, -8.75, 0.00, 0.00, -8.75, -8.75, 0.00, 0.00, -8.75, 0.00, 0.00, -8.75, -8.75, -8.75, -8.75, -8.75, 0.00, -8.75, -8.75, 0.00, -8.75, -8.75, -8.75, 0.00, -8.75, -8.75, 10.00},
		});
	}

}

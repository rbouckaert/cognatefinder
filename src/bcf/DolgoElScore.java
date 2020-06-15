package bcf;

import beast.core.Description;

@Description("DOLGO Sound-Class Model for use with ASJP Data "
		+ "Sound-Class model based on Dolgopolsky (1986) "
		+ "Sourced from lingpy/data/models/dolgo_el.")
public class DolgoElScore extends Score {

	public DolgoElScore() {	
		/* V */ process(12," a, e, i, o, u, E, 3");
		/* K */ process(5," k, g, x, X, C, c, q, j");
		/* P */ process(8," p, b, f, v");
		/* H */ process(3," H, h, 7, !");
		/* J */ process(4," y");
		/* M */ process(6," m");
		/* N */ process(7," n, 5, N, 4");
		/* S */ process(10," s, z, S, Z");
		/* R */ process(9," r, l, R, L");
		/* T */ process(11," t, d, 8, T");
		/* W */ process(13," w");
		/* + */ process(0," +");
		/* 1 */ process(2," 1");
		/* _ */ process(15," _");
	setScore(new double[][] {
		/* +	*/ {0.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -100.00, -5.00, -100.00},
		/* 0	*/ {-100.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
		/* 1	*/ {-100.00, 0.00, 2.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, 0.00, -20.00},
		/* H	*/ {-100.00, 0.00, -20.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* J	*/ {-100.00, 0.00, -20.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* K	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* M	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* N	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* P	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* R	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* S	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00, -10.00, 0.00, 0.00, -20.00},
		/* T	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, -10.00, 0.00, 0.00, -20.00},
		/* V	*/ {-100.00, 0.00, -20.00, -10.00, -10.00, -10.00, -10.00, -10.00, -10.00, -10.00, -10.00, -10.00, 5.00, -10.00, 0.00, -20.00},
		/* W	*/ {-100.00, 0.00, -20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -10.00, 10.00, 0.00, -20.00},
		/* X	*/ {-5.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
		/* _	*/ {-100.00, 0.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, -20.00, 0.00, 0.00},
		});
	}

}
package main.java.edu.usc.epigenome.uecgatk.bissnp;

import org.broadinstitute.gatk.utils.BaseUtils;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated massively
 * parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform Copyright (C) <2011>
 * <Yaping Liu: lyping1986@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program. If
 * not, see <http://www.gnu.org/licenses/>.
 */

public class BaseUtilsMore {
	// list of IUPAC code used for DNA nuleotide
	public final static byte A = (byte) 'A';
	public final static byte B = (byte) 'B'; // C,G,T
	public final static byte C = (byte) 'C';
	public final static byte D = (byte) 'D'; // A,G,T
	public final static byte G = (byte) 'G';
	public final static byte H = (byte) 'H'; // A,C,T
	public final static byte K = (byte) 'K'; // G,T
	public final static byte M = (byte) 'M'; // A,C
	public final static byte N = (byte) 'N'; // A,C,G,T
	public final static byte R = (byte) 'R'; // A,G
	public final static byte S = (byte) 'S'; // G,C
	public final static byte T = (byte) 'T';
	public final static byte V = (byte) 'V'; // A,C,G

	public final static byte W = (byte) 'W'; // A,T
	public final static byte X = (byte) 'X'; // to represent converted cytosine:
												// reference cytosine, but read
												// base is 'T'

	public final static byte Y = (byte) 'Y'; // C,T

	/**
	 * 
	 * @param contig
	 *            give contig string
	 * @return contigIndex return a number represent this contig
	 */
	static public int contigToIndex(String contig) {
		byte contigs = (byte) contig.charAt(3); // a bug here
		switch (contigs) {
		case 'X':
		case 'x':
			return 23;
		case 'Y':
		case 'y':
			return 24;
		case 'M':
		case 'm':
			return 25;
		default:
			return Integer.parseInt(Byte.toString(contigs));
		}
	}

	/**
	 * Creating a byte array and passing it to the String constructor
	 * 
	 * @param b
	 *            give byte
	 * @return String return a String
	 */
	static public String convertByteToString(byte b) {
		return new String(new byte[] { b });

	}

	static public byte convertTo2ndEnd(byte base) {
		switch (base) { // the wildcard character counts as an A
		case 'A':
		case 'a':
			return 'B';

		case 'C':
		case 'c':
			return 'D';

		case 'G':
		case 'g':
			return 'H';

		case 'T':
		case 't':
			return 'U';

		case 'X':
		case 'x':
			return 'Y';

		default:
			return base;
		}
	}

	// todo: better to get a bissnp engine based mismatched judges rather than
	// reference compared in the future
	/**
	 * It is a mismatch in Bisulfite space or not
	 * 
	 * @param refBase
	 * @param base
	 * @param negStrand
	 * @param secondPair
	 * @return mismatchOrNotInBisulfiteSpace
	 */
	static public boolean isBisulfiteMismatch(byte refBase, byte base, boolean negStrand, boolean secondPair) {
		if (!negStrand) {
			if (secondPair) {
				if (BaseUtils.basesAreEqual(refBase, G) && BaseUtils.basesAreEqual(base, A)) {
					return false;
				} else {
					return true;
				}
			} else {
				if (BaseUtils.basesAreEqual(refBase, C) && BaseUtils.basesAreEqual(base, T)) {
					return false;
				} else {
					return true;
				}
			}

		} else {
			if (secondPair) {
				if (BaseUtils.basesAreEqual(refBase, C) && BaseUtils.basesAreEqual(base, T)) {
					return false;
				} else {
					return true;
				}
			} else {
				if (BaseUtils.basesAreEqual(refBase, G) && BaseUtils.basesAreEqual(base, A)) {
					return false;
				} else {
					return true;
				}
			}

		}
	}

	static public boolean isRegularBase(byte base) {
		return simpleBaseToBaseIndex(base) != -1;
	}

	// convert IUPAC code pattern to its complement code (in reverse strand)
	static public byte iupacCodeComplement(byte base) {
		base = toUpperCase(base);
		switch (base) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'R':
			return 'Y';
		case 'Y':
			return 'R';
		case 'S':
			return 'S';
		case 'W':
			return 'W';
		case 'K':
			return 'M';
		case 'M':
			return 'K';
		case 'B':
			return 'V';
		case 'H':
			return 'D';
		case 'D':
			return 'H';
		case 'V':
			return 'B';
		default:
			return base;
		}
	}

	// find the observe base is matched IUPAC code pattern or not, consider
	// cytosine methylation, so T will be equal to C..
	static public boolean iupacCodeEqual(byte pattern, byte observe, boolean negStrand) {
		pattern = toUpperCase(pattern);
		observe = toUpperCase(observe);
		switch (observe) {
		case 'A':
			switch (pattern) {
			case 'A':
			case 'R':
			case 'W':
			case 'M':
			case 'H':
			case 'V':
			case 'D':
			case 'N':
				return true;
			case 'G':
			case 'S':
			case 'K':
			case 'B':
				if (negStrand) {
					return true;
				} else {
					return false;
				}
			default:
				return false;
			}

		case 'C':
			switch (pattern) {
			case 'C':
			case 'Y':
			case 'S':
			case 'M':
			case 'H':
			case 'V':
			case 'B':
			case 'N':
				return true;
			default:
				return false;
			}

		case 'G':
			switch (pattern) {
			case 'G':
			case 'R':
			case 'S':
			case 'K':
			case 'D':
			case 'V':
			case 'B':
			case 'N':
				return true;
			default:
				return false;
			}

		case 'T':
			switch (pattern) {
			case 'T':
			case 'W':
			case 'K':
			case 'D':
			case 'Y':
			case 'H':
			case 'B':
			case 'N':
				return true;
			case 'C':
			case 'S':
			case 'M':
			case 'V':
				if (negStrand) {
					return false;
				} else {
					return true;
				}
			default:
				return false;
			}

		case 'N':
			return true;
		default:
			System.err.println("error! wrong observed base!");
			return false;
		}

	}

	// find the observe base is matched IUPAC code pattern or not, but not
	// consider cytosine methylation. this is the one i use
	static public boolean iupacCodeEqualNotConsiderMethyStatus(byte pattern, byte observe) {
		pattern = toUpperCase(pattern);
		observe = toUpperCase(observe);
		if (BaseUtils.basesAreEqual(observe, A)) {
			switch (pattern) {
			case 'A':
			case 'R':
			case 'W':
			case 'M':
			case 'H':
			case 'V':
			case 'D':
			case 'N':
				return true;

			default:
				return false;
			}
		} else if (BaseUtils.basesAreEqual(observe, C)) {
			switch (pattern) {
			case 'C':
			case 'Y':
			case 'S':
			case 'M':
			case 'H':
			case 'V':
			case 'B':
			case 'N':
				return true;
			default:
				return false;
			}
		} else if (BaseUtils.basesAreEqual(observe, G)) {
			switch (pattern) {
			case 'G':
			case 'R':
			case 'S':
			case 'K':
			case 'D':
			case 'V':
			case 'B':
			case 'N':
				return true;
			default:
				return false;
			}
		} else if (BaseUtils.basesAreEqual(observe, T)) {
			switch (pattern) {
			case 'T':
			case 'W':
			case 'K':
			case 'D':
			case 'Y':
			case 'H':
			case 'B':
			case 'N':
				return true;

			default:
				return false;
			}
		} else {
			return true;
		}
	}

	static public boolean iupacCodeEqualNotConsiderMethyStatus(String pattern, String observe) {
		byte[] patterns = pattern.getBytes();
		byte[] observes = observe.getBytes();
		if (patterns.length != observes.length)
			return false;
		for (int i = 0; i < observes.length; i++) {
			if (!iupacCodeEqualNotConsiderMethyStatus(patterns[i], observes[i])) {
				return false;
			}
		}
		return true;
	}

	// observed pattern and search pattern has some overlap, like 'H' overlapped with 'B'
	static public boolean iupacCodeOverlap(byte pattern, byte observe) {
		pattern = toUpperCase(pattern);
		observe = toUpperCase(observe);
		if (observe == N || pattern == N) {
			return true;
		} else if (pattern == observe) {
			return true;
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.R)) {
			switch (observe) {
			case 'C':
			case 'T':
			case 'Y':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.Y)) {
			switch (observe) {
			case 'A':
			case 'G':
			case 'R':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.S)) {
			switch (observe) {
			case 'A':
			case 'T':
			case 'W':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.W)) {
			switch (observe) {
			case 'G':
			case 'C':
			case 'S':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.K)) {
			switch (observe) {
			case 'A':
			case 'C':
			case 'M':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.M)) {
			switch (observe) {
			case 'G':
			case 'T':
			case 'K':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.B)) {
			switch (observe) {
			case 'A':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.H)) {
			switch (observe) {
			case 'G':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.V)) {
			switch (observe) {
			case 'T':
				return false;
			default:
				return true;
			}
		} else if (BaseUtils.basesAreEqual(pattern, BaseUtilsMore.D)) {
			switch (observe) {
			case 'C':
				return false;
			default:
				return true;
			}
		} else {

			System.err.println("error! wrong observed base!");
			return false;
		}
	}

	static public boolean iupacCodeOverlap(String pattern, String observe) {
		byte[] patterns = pattern.getBytes();
		byte[] observes = observe.getBytes();
		if (patterns.length != observes.length)
			return false;
		for (int i = 0; i < observes.length; i++) {
			if (!iupacCodeOverlap(patterns[i], observes[i])) {
				return false;
			}
		}
		return true;
	}

	static public byte makeIupacCodeFrom2Byte(byte a, byte b) {
		a = toUpperCase(a);
		b = toUpperCase(b);
		if (BaseUtils.basesAreEqual(a, b)) {
			return a;
		} else {
			if (BaseUtils.basesAreEqual(a, A)) {
				switch (b) {

				case 'C':
					return 'M';
				case 'G':
					return 'R';
				case 'T':
					return 'W';
				default:
					return 'N';
				}
			} else if (BaseUtils.basesAreEqual(a, C)) {
				switch (b) {
				case 'A':
					return 'M';

				case 'G':
					return 'S';
				case 'T':
					return 'Y';
				default:
					return 'N';
				}
			} else if (BaseUtils.basesAreEqual(a, G)) {
				switch (b) {
				case 'A':
					return 'R';
				case 'C':
					return 'S';

				case 'T':
					return 'K';
				default:
					return 'N';
				}
			} else if (BaseUtils.basesAreEqual(a, T)) {
				switch (b) {
				case 'A':
					return 'W';
				case 'C':
					return 'Y';
				case 'G':
					return 'K';

				default:
					return 'N';
				}
			} else {
				return 'N';
			}
		}

	}

	static public String makeIupacCodeFrom2String(String a, String b) {
		byte[] as = a.getBytes();
		byte[] bs = b.getBytes();
		String c = "";
		for (int i = 0; i < Math.max(as.length, bs.length); i++) {
			if (i < as.length & i < bs.length) {
				c += (char) makeIupacCodeFrom2Byte(as[i], bs[i]);
			} else if (i < as.length) {
				c += (char) as[i];
			} else {
				c += (char) bs[i];
			}
		}
		return c;
	}

	static public boolean searchIupacPatternFromBases(byte[] searchPatterns, byte[] bases, boolean complement) {

		for (int i = 0; i < searchPatterns.length; i++) {
			if (i >= bases.length)
				return false;
			byte base = bases[i];
			if (complement)
				base = BaseUtilsMore.iupacCodeComplement(base);
			if (!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(searchPatterns[i], base))
				return false;
		}
		return true;
	}

	/**
	 * Converts a simple base to a base index
	 * 
	 * @param base
	 *            [AaCcGgTt]
	 * @return 0, 1, 2, 3, 4, or -1 if the base can't be understood
	 */
	static public int simpleBaseToBaseIndex(byte base) {
		switch (base) {
		case '*': // the wildcard character counts as an A
		case 'A':
		case 'a':
			return 0;

		case 'C':
		case 'c':
			return 1;

		case 'G':
		case 'g':
			return 2;

		case 'T':
		case 't':
			return 3;

		case 'X':
		case 'x':
			return 4;
		case 'B':
		case 'b':
			return 5;

		case 'D':
		case 'd':
			return 6;

		case 'H':
		case 'h':
			return 7;

		case 'U':
		case 'u':
			return 8;

		case 'Y':
		case 'y':
			return 9;
		default:
			return -1;
		}
	}

	// reverse the read's base order
	static public byte[] simpleReverse(byte[] bases) {
		byte[] rcbases = new byte[bases.length];

		for (int i = 0; i < bases.length; i++) {
			rcbases[i] = bases[bases.length - 1 - i];
		}

		return rcbases;
	}

	// reverse the read's base order and make it to be complementary
	static public byte[] simpleReverseIupacCodeComplement(byte[] bases) {
		byte[] rcbases = new byte[bases.length];

		for (int i = 0; i < bases.length; i++) {
			rcbases[i] = iupacCodeComplement(bases[bases.length - 1 - i]);
		}

		return rcbases;
	}
	
	static public byte[] simpleIupacCodeComplement(byte[] bases) {
        byte[] cbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
        	cbases[i] = iupacCodeComplement(bases[i]);
        }

        return cbases;
    }

	// convert to IUPAC code in NOMe-seq mode
	static public byte toIupacCodeNOMeSeqMode(byte base, int pos) {
		base = toUpperCase(base);
		if (pos == 1) {
			switch (base) {
			case 'A':
				return 'W';
			case 'T':
				return 'W';
			default:
				return base;
			}
		} else if (pos == 3) {
			switch (base) {
			case 'A':
				return 'H';
			case 'T':
				return 'H';
			case 'C':
				return 'H';
			default:
				return base;
			}
		} else {
			return base;
		}

	}

	static public byte toUpperCase(byte b) {
		return (byte) Character.toUpperCase((char) b);
	}

	// convert input bases to uppercase
	static public byte[] toUpperCase(byte[] in) {
		byte[] out = new byte[in.length];
		for (int i = 0; i < in.length; i++) {
			out[i] = toUpperCase(in[i]);
		}
		return out;
	}

}

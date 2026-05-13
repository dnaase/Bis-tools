#!/usr/bin/perl -w
## scNOMe-HiC / non-directional NOMe-seq wrapper for BisSNP.
##
## Companion to `bissnp_easy_usage.pl`, but tailored for non-directional
## bisulfite libraries (single-cell NOMe-HiC, scNMT-seq variants, etc.) where
## R2 reverse reads carry XG=CT and the BAM `is_reverse` flag does NOT
## reflect the bisulfite-conversion frame. Pairs with the XG-aware filter
## chain in BisSNP ≥1.1.
##
## Outputs (with default --prefix derived from BAM):
##   <prefix>.cyt.vcf                                  cytosine VCF (raw, per-base methylation)
##   <prefix>.snp.vcf                                  SNP VCF
##   <prefix>.GCH.txt                                  per-read GCH methylation status
##   <prefix>.HCG.txt                                  per-read HCG methylation status
##   <prefix>.cyt.filtered.sort.vcf                    filtered + sorted cytosine VCF
##   <prefix>.cyt.filtered.sort.GCH.6plus2.bed         per-base GCH methylation BED (combined strand)
##   <prefix>.cyt.filtered.sort.HCG.6plus2.bed         per-base HCG methylation BED (combined strand)
##   <prefix>.cyt.filtered.sort.GCH.strand.6plus2.bed  per-base GCH methylation BED (per strand)
##   <prefix>.cyt.filtered.sort.HCG.strand.6plus2.bed  per-base HCG methylation BED (per strand)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.GCH.bedgraph             methylation bedGraph (GCH)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.HCG.bedgraph             methylation bedGraph (HCG)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.GCH.coverage.bedgraph    C+T read coverage bedGraph (GCH)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.HCG.coverage.bedgraph    C+T read coverage bedGraph (HCG)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.GCH.bw                   methylation bigWig (GCH)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.HCG.bw                   methylation bigWig (HCG)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.GCH.coverage.bw          coverage bigWig (GCH)
##   <prefix>.cyt.filtered.sort.BisSNP-<v>.HCG.coverage.bw          coverage bigWig (HCG)
##   (and GCG / HCH variants if --allC)

use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

sub usage {
    print <<'USAGE';

Usage:
  perl bissnp_nomehic_usage.pl [Options] BISSNP_JAR INPUT_BAM REF_GENOME DBSNP_VCF

NON-DIRECTIONAL NOMe-seq / scNOMe-HiC wrapper. Runs BisulfiteGenotyper with
the flag set required for libraries where the XG tag (not BAM is_reverse)
defines the bisulfite-conversion frame, then post-processes the outputs
into per-base BED files. Per-read GCH/HCG methylation tables are emitted
inline by the genotyper.

PREREQUISITES
  - BisSNP jar version >= 1.1 (for XG-aware filter chain).
  - Java 8 on PATH or supplied via --java8.
  - INPUT_BAM coordinate-sorted, indexed, with MD tag. For scNOMe-HiC,
    BAM should already have Fix-D preprocessing applied (run
    scripts/strip_paired_all.py before this wrapper) — without it, the
    GATK 3.8 adaptor-clip filter silently drops the reverse mate of
    chimeric Hi-C ligations.
  - Helper scripts found alongside this wrapper or in ../utils/:
      sortByRefAndCor.pl   vcf2bed6plus2.pl   vcf2bed6plus2.strand.pl
      vcf2bedGraph.pl      vcf2coverage.pl
  - bedGraphToBigWig (UCSC) on PATH or supplied via --bgtobw, for bigWig
    output. Use --noBigwig to skip bigWig conversion.
  - REF .fai sidecar (used to auto-derive chrom.sizes for bigWig).

POSITIONAL ARGUMENTS
  BISSNP_JAR    Path to BisSNP jar (e.g. Bis-SNP.latest.jar or BisSNP-1.1.jar).
  INPUT_BAM     Input BAM. Must contain ReadGroup tag.
  REF_GENOME    Reference .fa with .fai sidecar.
  DBSNP_VCF     dbSNP VCF (e.g. Homo_sapiens_assembly38.dbsnp138.vcf).

OPTIONS
  --interval FILE   restrict analysis to a BED of intervals (e.g. chr22.bed)
                    or a single chromosome name passed via -L (e.g. chr22).
  --prefix STR      output filename prefix (default: BAM basename without .bam)
  --mem NUM         Java heap GB (default: 20)
  --nt NUM          BisSNP threads (default: 1; values > 1 are NOT safe with
                    NOMESEQ_MODE — keep at 1)
  --mmq NUM         min mapping quality (default: 30)
  --mbq NUM         min base quality (default: 5)
  --qual NUM        stand_call_conf (default: 1 — lowered from BisSNP default 20
                    because single-cell coverage is sparse)
  --mm NUM          max non-bisulfite mismatch fraction (default: 0.3)
  --minPatConv NUM  min CH-pattern conversion ratio (default: 1.0 = disabled)
  --minConv NUM     min 5'-converted Cs required (default: 0)
  --allC            also produce per-base BED / bedGraph / bigWig for GCG
                    and HCH (default: only GCH and HCG)
  --skipBed         skip the BED / bedGraph / bigWig post-processing chain
                    (run BisSNP only)
  --noStrand        skip the per-strand 6plus2 BED step (combined-strand only)
  --noBigwig        skip bedGraph→bigWig conversion (bedGraphs still produced)
  --bgtobw PATH     path to bedGraphToBigWig binary (default: search PATH)
  --java8 PATH      path to Java 8 binary
                    (default: /software/java/jdk1.8.0_191/bin/java)
  --dryRun          print commands without executing
  --help            show this message

EXAMPLE (Quest, scNOMe-HiC test cell)
  perl bissnp_nomehic_usage.pl \
      --interval chr22 --prefix bs_pr --mem 20 \
      /home/jmj7858/epifluidlab/software/Bis-tools/Bis-SNP/Bis-SNP.latest.jar \
      batch2.scD24.TAGCTT.dedup.calmd.allunpaired.bam \
      /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
      /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf

USAGE
    exit 1;
}

my $interval   = "";
my $prefix     = "";
my $mem        = 20;
my $nt         = 1;
my $mmq        = 30;
my $mbq        = 5;
my $qual       = 1;
my $mm         = 0.3;
my $minPatConv = 1.0;
my $minConv    = 0;
my $allC       = 0;
my $skipBed    = 0;
my $noStrand   = 0;
my $noBigwig   = 0;
my $bgtobw     = "";
my $java8      = "/software/java/jdk1.8.0_191/bin/java";
my $dryRun     = 0;
my $help       = 0;

GetOptions(
    "interval=s"   => \$interval,
    "prefix=s"     => \$prefix,
    "mem=i"        => \$mem,
    "nt=i"         => \$nt,
    "mmq=i"        => \$mmq,
    "mbq=i"        => \$mbq,
    "qual=i"       => \$qual,
    "mm=f"         => \$mm,
    "minPatConv=f" => \$minPatConv,
    "minConv=i"    => \$minConv,
    "allC"         => \$allC,
    "skipBed"      => \$skipBed,
    "noStrand"     => \$noStrand,
    "noBigwig"     => \$noBigwig,
    "bgtobw=s"     => \$bgtobw,
    "java8=s"      => \$java8,
    "dryRun"       => \$dryRun,
    "help"         => \$help,
) or usage();

usage() if $help || scalar(@ARGV) != 4;

my ($BISSNP, $BAM, $REF, $DBSNP) = @ARGV;

foreach my $pair (["BISSNP_JAR",$BISSNP], ["INPUT_BAM",$BAM],
                  ["REF_GENOME",$REF], ["DBSNP_VCF",$DBSNP]) {
    my ($name, $f) = @$pair;
    die "$name not found: $f\n" unless -e $f;
}
die "REF .fai sidecar not found: ${REF}.fai\n" unless -e "${REF}.fai";

unless ($prefix) {
    $prefix = $BAM;
    $prefix =~ s|\.bam$||;
}

my $script_dir = dirname(abs_path(__FILE__));
# Helpers live in upstream Bis-tools/utils/; some forks copy them next to
# the wrapper. Search both locations.
sub find_helper {
    my $name = shift;
    foreach my $d ("$script_dir", "$script_dir/../utils") {
        my $p = "$d/$name";
        return $p if -e $p;
    }
    return undef;
}
my $SORT_VCF     = find_helper("sortByRefAndCor.pl");
my $VCF2BED      = find_helper("vcf2bed6plus2.pl");
my $VCF2BED_STR  = find_helper("vcf2bed6plus2.strand.pl");
my $VCF2BG       = find_helper("vcf2bedGraph.pl");
my $VCF2COV      = find_helper("vcf2coverage.pl");

unless ($skipBed) {
    die "Helper sortByRefAndCor.pl not found near $script_dir\n"     unless $SORT_VCF;
    die "Helper vcf2bed6plus2.pl not found near $script_dir\n"       unless $VCF2BED;
    die "Helper vcf2bed6plus2.strand.pl not found near $script_dir\n"
        if !$noStrand && !$VCF2BED_STR;
    die "Helper vcf2bedGraph.pl not found near $script_dir\n"        unless $VCF2BG;
    die "Helper vcf2coverage.pl not found near $script_dir\n"        unless $VCF2COV;
}

# Discover bedGraphToBigWig binary (only matters if we'll actually do bw step).
sub which {
    my $prog = shift;
    foreach my $d (split /:/, ($ENV{PATH} || "")) {
        my $p = "$d/$prog";
        return $p if -x $p;
    }
    return undef;
}
unless ($skipBed || $noBigwig) {
    if (!$bgtobw) {
        $bgtobw = which("bedGraphToBigWig");
    }
    die "bedGraphToBigWig not found on PATH; pass --bgtobw PATH or use --noBigwig\n"
        unless $bgtobw && -x $bgtobw;
}

my $JAVA = "$java8 -Xmx${mem}G";

my $cyt_vcf      = "${prefix}.cyt.vcf";
my $snp_vcf      = "${prefix}.snp.vcf";
my $gch_per_read = "${prefix}.GCH.txt";
my $hcg_per_read = "${prefix}.HCG.txt";
my $cyt_sort     = "${prefix}.cyt.sort.vcf";
my $snp_sort     = "${prefix}.snp.sort.vcf";
my $cyt_filt     = "${prefix}.cyt.filtered.sort.vcf";

sub stamp { my $t = localtime(); return "[$t]"; }

sub run {
    my ($desc, $cmd) = @_;
    print STDERR stamp() . "  $desc\n$cmd\n";
    return if $dryRun;
    my $rc = system($cmd);
    die stamp() . "  FAILED ($desc): exit $rc\n" if $rc != 0;
}

sub step_genotyper {
    my $cmd = "$JAVA -jar $BISSNP \\\n"
            . "    -R $REF \\\n"
            . "    -I $BAM \\\n"
            . "    -D $DBSNP \\\n"
            . "    -T BisulfiteGenotyper \\\n"
            . "    -vfn1 $cyt_vcf \\\n"
            . "    -vfn2 $snp_vcf \\\n"
            . "    -gchreads $gch_per_read \\\n"
            . "    -cpgreads $hcg_per_read \\\n"
            . "    -out_modes NOMESEQ_MODE -sm GM \\\n"
            . "    -stand_call_conf $qual \\\n"
            . "    -nonDirectional -badMate \\\n"
            . "    -mm $mm -minPatConv $minPatConv -minConv $minConv \\\n"
            . "    -mmq $mmq -mbq $mbq -nt $nt";
    $cmd .= " \\\n    -L $interval" if $interval;
    run("BisulfiteGenotyper (NOMESEQ_MODE, non-directional)", $cmd);
}

sub step_sort {
    foreach my $vcf ($cyt_vcf, $snp_vcf) {
        my $out = $vcf;
        $out =~ s|\.vcf$|.sort.vcf|;
        run("Sort $vcf",
            "perl $SORT_VCF --k 1 --c 2 --tmp ./ $vcf ${REF}.fai > $out");
    }
}

sub step_filter {
    my $cmd = "$JAVA -jar $BISSNP \\\n"
            . "    -R $REF -T VCFpostprocess \\\n"
            . "    -C GCH -C HCH -C GCG -C HCG \\\n"
            . "    -oldVcf $cyt_sort -snpVcf $snp_sort \\\n"
            . "    -newVcf $cyt_filt \\\n"
            . "    -o $cyt_filt.cpgSummary.txt -minCT 1 -qual $qual";
    $cmd .= " \\\n    -L $interval" if $interval;
    run("VCFpostprocess (NOMe contexts)", $cmd);
}

sub step_to_bed {
    # Combined-strand 6plus2 BED via upstream vcf2bed6plus2.pl.
    my @ctx = $allC ? qw(GCH HCG GCG HCH) : qw(GCH HCG);
    foreach my $c (@ctx) {
        run("vcf2bed6plus2 $c (combined strand)",
            "perl $VCF2BED --only_good_call $cyt_filt $c");
    }
}

sub step_to_bed_strand {
    # Per-strand 6plus2 BED via vcf2bed6plus2.strand.pl. Output names get
    # a `.strand.` infix (e.g. *.GCH.strand.6plus2.bed) so they coexist
    # with the combined-strand BEDs from step_to_bed.
    my @ctx = $allC ? qw(GCH HCG GCG HCH) : qw(GCH HCG);
    foreach my $c (@ctx) {
        run("vcf2bed6plus2.strand $c (per strand)",
            "perl $VCF2BED_STR $cyt_filt $c");
    }
}

sub step_bedgraph {
    # Methylation bedGraph (chr, start, end, methy%, numCT) and
    # coverage bedGraph (chr, start, end, C+T reads).
    # Output filenames bake in the BisSNP version from the VCF header,
    # so we can't fully predict them here; the bigWig step globs.
    my @ctx = $allC ? qw(GCH HCG GCG HCH) : qw(GCH HCG);
    foreach my $c (@ctx) {
        run("vcf2bedGraph $c (methylation)",
            "perl $VCF2BG $cyt_filt $c");
        run("vcf2coverage $c (C+T reads)",
            "perl $VCF2COV $cyt_filt $c");
    }
}

sub step_bigwig {
    # bedGraph -> bigWig. Requires:
    #   (1) a chrom.sizes file (auto-derived from ${REF}.fai),
    #   (2) bedGraph input with no `track` header, no '.' value rows,
    #       only chromosomes present in chrom.sizes, sorted by chr+pos.
    # bedGraphToBigWig is strict and segfaults / errors on any of the above.

    my $chrom_sizes = "${prefix}.chrom.sizes";
    run("Derive chrom.sizes from ${REF}.fai",
        "awk 'BEGIN{OFS=\"\\t\"} {print \$1,\$2}' ${REF}.fai > $chrom_sizes");

    # Find every bedGraph produced from $cyt_filt (both methylation and
    # coverage variants; both have the BisSNP version baked into their
    # filenames so we glob).
    my $cyt_stem = $cyt_filt;
    $cyt_stem =~ s|\.vcf$||;
    my @bgs = glob("${cyt_stem}.*.bedgraph");
    if (!@bgs && !$dryRun) {
        warn stamp() . "  No bedGraphs found matching ${cyt_stem}.*.bedgraph — skipping bigWig step\n";
        return;
    }

    foreach my $bg (@bgs) {
        my $sorted = $bg;
        $sorted =~ s|\.bedgraph$|.sorted.bedgraph|;
        my $bw    = $bg;
        $bw    =~ s|\.bedgraph$|.bw|;

        # Strip `track` header, drop rows with chrom not in chrom.sizes
        # or value == '.' (low-CT sites in methylation bedGraph), then sort.
        my $prep = qq(awk 'NR==FNR{c[\$1]=1;next} /^track/{next} (\$1 in c) && \$4!="."' )
                 . qq($chrom_sizes $bg | sort -k1,1 -k2,2n -T ./ > $sorted);
        run("Prepare $bg for bigWig", $prep);

        run("bedGraphToBigWig $sorted",
            "$bgtobw $sorted $chrom_sizes $bw");

        # Drop the intermediate sorted bedGraph; keep the original bedGraph
        # and the .bw.
        run("rm $sorted", "rm -f $sorted") unless $dryRun;
    }
}

print STDERR stamp() . "  Output prefix: $prefix\n";

step_genotyper();

unless ($skipBed) {
    step_sort();
    step_filter();
    step_to_bed();
    step_to_bed_strand() unless $noStrand;
    step_bedgraph();
    step_bigwig()        unless $noBigwig;
}

print STDERR stamp() . "  Done. Per-read tables: $gch_per_read, $hcg_per_read\n";
unless ($skipBed) {
    print STDERR stamp() . "  Per-base BEDs:      ${prefix}.cyt.filtered.sort.{GCH,HCG}*.bed\n";
    print STDERR stamp() . "  Methy/cov bedGraph: ${prefix}.cyt.filtered.sort.*.{GCH,HCG}*.bedgraph\n";
    print STDERR stamp() . "  Methy/cov bigWig:   ${prefix}.cyt.filtered.sort.*.{GCH,HCG}*.bw\n"
        unless $noBigwig;
}

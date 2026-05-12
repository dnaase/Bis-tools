# BisSNP installation

Three paths, in order of effort.

| Path | When to use | Time |
|---|---|---|
| (A) Use the prebuilt jar | You just want to run BisSNP | ~1 min |
| (B) Rebuild from source | You modified the source, or want byte-level provenance | ~30 sec compile |
| (C) Use a release-tagged version | You need to pin to a specific tag (e.g. `v1.1`) | ~1 min |

---

## (A) Prebuilt jar (recommended for users)

The repo tracks a prebuilt fatjar at `Bis-SNP/Bis-SNP.latest.jar`. It bundles
GATK 3.8 + htsjdk + all other compile-time dependencies, so you don't need
Maven or to download anything else.

**Requirements**
- Java 8 (GATK 3.8 hard-requires it; will NOT run on Java 11+).

```bash
git clone https://github.com/dnaase/Bis-tools.git
cd Bis-tools
/path/to/java8 -jar Bis-SNP/Bis-SNP.latest.jar --help | head -3
```

That's it — `Bis-SNP.latest.jar` is the canonical artifact. There is no
version-pinned `BisSNP-X.Y.jar` in the repo; the version is recorded in the
git tag of the commit that produced the jar, and in the runtime banner the
jar prints at startup.

---

## (B) Rebuild from source

Use this if you modified the source or want to verify byte equivalence
against the tracked fatjar.

**Requirements**
- Java 8 (build + runtime)
- Maven 3.x
- A working internet connection on first build (Maven downloads `htsjdk`,
  `commons-lang3`, `crc-64` from Maven Central into `~/.m2/repository/`)
- **GATK 3.8 jar manually placed at**
  `Bis-SNP/lib/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar`

GATK 3.8 is declared in `pom.xml` as `<scope>system</scope>` with a
relative path, so it cannot be auto-downloaded. Broad's redistribution
restrictions are why it isn't bundled here.

**Get GATK 3.8 from Broad**

```bash
cd Bis-SNP
wget https://github.com/broadgsa/gatk-protected/releases/download/3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
mkdir -p lib
tar -xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 -C lib/
ls lib/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar  # should exist
```

**Build**

```bash
cd Bis-SNP
mvn -DskipTests package
# Output: target/bissnp-1.0.0-jar-with-dependencies.jar (fatjar)
# (filename uses pom version 1.0.0 — runtime banner says BisSNP-1.0.1 +
#  whatever XG patches are present; the version in the filename and the
#  runtime banner are intentionally distinct here.)
```

**Overlay-build alternative**

For incremental rebuilds where you want the same fatjar layout as the
tracked `Bis-SNP.latest.jar`, use the overlay procedure that drops compiled
classes onto a pre-extracted GATK fatjar tree. This is faster than a full
`mvn package` and matches the deploy artifact byte-for-byte. See the test
runner in `bistools_pr_test/run_pr_test.sh` (in the development repo, not
shipped here) for an example.

---

## (C) Pin to a release tag

For reproducibility, check out a specific tag:

```bash
git clone https://github.com/dnaase/Bis-tools.git
cd Bis-tools
git checkout v1.1                    # or whatever tag you want
ls Bis-SNP/Bis-SNP.latest.jar        # the jar that shipped with that tag
```

Fork releases with version-named asset downloads also exist, e.g.

```
https://github.com/yoshihiko1218/Bis-tools/releases/tag/v1.1
```

The asset `BisSNP-1.1.jar` on that page is identical to the
`Bis-SNP.latest.jar` in the repo at commit `c127b31..035051a`.

---

## Java 8 on cluster systems

If your cluster doesn't expose Java 8 by default, look for it under
`/software/java/` or load via your module system. Quest's path, for
reference:

```bash
/software/java/jdk1.8.0_191/bin/java -jar Bis-SNP/Bis-SNP.latest.jar ...
```

## Verifying the install

```bash
/path/to/java8 -jar Bis-SNP/Bis-SNP.latest.jar --help 2>&1 | head -10
```

Expected banner:

```
The BisSNP-1.0.1, Compiled 2018/02/19 05:43:50
Based on The Genome Analysis Toolkit (GATK) v3.8-1-0-gf15c1c3ef
```

(The `1.0.1` and `2018/02/19` in the banner are historical strings baked
into the source; the behavior is whatever the current commit's source
implements.)

## Common gotchas

- **`UnsupportedClassVersionError`** at startup — Java >= 11 detected. Use Java 8.
- **`NoClassDefFoundError: org/broadinstitute/...`** during `mvn package` — GATK 3.8 jar missing from `Bis-SNP/lib/`. See section (B).
- **Process killed silently with exit code 9** on a cluster — you're on a login node. Switch to a compute node via `salloc` / `sbatch`.
- **`Unknown option: seperate_strand`** from `vcf2bed6plus2.pl` — old wrappers passed a flag the upstream helper doesn't support. The current `bissnp_nomehic_usage.pl` no longer passes it.

## Test suite for barena.
## Run with: nimble test   (or: nim c -r tests/test_barena.nim)

import std/[unittest, os, osproc, tables, sets]
import "../src/taxonomy"
import "../src/krakenio"
import "../src/krakensignatures"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

const
  Data       = "tests/data"
  InspectDb  = Data / "k2_inspect_db_output.txt"
  KrakenRaw  = Data / "kraken_output.raw"
  R1         = Data / "A01_R1.fq.gz"
  R2         = Data / "A01_R2.fq.gz"

let barena =
  if fileExists("./barena"): "./barena"
  elif fileExists("./bin/barena"): "./bin/barena"
  elif fileExists("./src/barena"): "./src/barena"
  else: ""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

proc countReads(path: string): int =
  ## Count records in a 4-line FASTQ file.
  if not fileExists(path): return 0
  var n = 0
  for _ in lines(path): inc n
  n div 4

proc run(args: string): tuple[stderr: string, code: int] =
  let (outp, code) = execCmdEx(barena & " " & args & " 2>&1")
  (outp, code)

proc freshDir(path: string) =
  if dirExists(path):
    removeDir(path)
  createDir(path)

# ---------------------------------------------------------------------------
# Taxonomy tests
# ---------------------------------------------------------------------------

suite "taxonomy parser":

  test "parses inspect-db without error":
    let tree = parseInspectDb(InspectDb)
    check tree.len > 0

  test "45409 (Lachnospiraceae) is a direct child of 45405 (Lachnospirales)":
    let tree = parseInspectDb(InspectDb)
    check 45405'u32 in tree
    check 45409'u32 in tree[45405'u32]

  test "45409 and 45448 are descendants of 45405":
    let tree = parseInspectDb(InspectDb)
    let desc = getDescendants(tree, @[45405'u32])
    check 45409'u32 in desc
    check 45448'u32 in desc

  test "45326 (Christensenellaceae) is a direct child of 45325 (Christensenellales)":
    let tree = parseInspectDb(InspectDb)
    check 45325'u32 in tree
    check 45326'u32 in tree[45325'u32]

  test "45328 and 45326 are descendants of 45325":
    let tree = parseInspectDb(InspectDb)
    let desc = getDescendants(tree, @[45325'u32])
    check 45326'u32 in desc
    check 45328'u32 in desc

  test "45405, 45409, 45448, 45325, 45326, 45328 are all descendants of 3":
    let tree = parseInspectDb(InspectDb)
    let desc = getDescendants(tree, @[3'u32])
    for id in [45405'u32, 45409'u32, 45448'u32, 45325'u32, 45326'u32, 45328'u32]:
      check id in desc

# ---------------------------------------------------------------------------
# krakenio unit tests
# ---------------------------------------------------------------------------

suite "kraken line parser":

  test "classified line returns true and correct taxid":
    var taxId: uint32
    let c = parseKrakenLine("C\tread1\t46230\t301|301\t3:10", taxId)
    check c == true
    check taxId == 46230'u32

  test "unclassified line returns false and taxid 0":
    var taxId: uint32
    let c = parseKrakenLine("U\tread1\t0\t301|301\t0:267", taxId)
    check c == false
    check taxId == 0'u32

  test "stripPairSuffix removes /1 and /2":
    check stripPairSuffix("read/1") == "read"
    check stripPairSuffix("read/2") == "read"
    check stripPairSuffix("read")   == "read"
    check stripPairSuffix("read/3") == "read/3"

# ---------------------------------------------------------------------------
# Kraken signature tests
# ---------------------------------------------------------------------------

suite "kraken signature parser and comparison":

  test "parses paired-end RLE signature":
    let sig = parseKrakenSignature("3:3 46409:2 3:75 |:| 0:267")
    check sig.isPaired()
    check sig.mateCount() == 2
    check sig.mates[0].kmerCount() == 80
    check sig.mates[1].kmerCount() == 267
    check sig.toKrakenString() == "3:3 46409:2 3:75 |:| 0:267"

  test "normalizes adjacent equal runs while parsing":
    let sig = parseKrakenSignature("3:1 3:2|:|A:1 A:2")
    check sig.toKrakenString() == "3:3 |:| A:3"

  test "extracts signature from full kraken line":
    let sig = parseKrakenLineSignature("C\tread1\t46230\t301|301\t3:2 0:1 |:| A:3")
    check sig.toKrakenString() == "3:2 0:1 |:| A:3"

  test "exact positional similarity keeps 0 and A as real tokens":
    let a = parseKrakenSignature("3:2 0:1 |:| A:2")
    let b = parseKrakenSignature("3:1 4:1 0:1 |:| A:1 0:1")
    check abs(positionalSimilarity(a, b) - 0.6) < 1e-9

  test "weighted jaccard compares composition and keeps mate identity":
    let a = parseKrakenSignature("3:2 0:1 |:| A:2")
    let b = parseKrakenSignature("3:1 4:1 0:1 |:| A:1 0:1")
    check abs(weightedJaccardSimilarity(a, b) - (3.0 / 7.0)) < 1e-9

  test "taxonomy-aware similarity gives siblings partial credit":
    var tree = initTable[uint32, seq[uint32]]()
    tree[1'u32] = @[2'u32, 3'u32]
    tree[2'u32] = @[4'u32, 5'u32]
    tree[4'u32] = @[6'u32]

    let taxonomy = buildTaxonomyIndex(tree)
    var lca: uint32

    check taxonomy.sameParent(4'u32, 5'u32)
    check taxonomy.lowestCommonAncestor(6'u32, 5'u32, lca)
    check lca == 2'u32
    check abs(taxonomicTokenSimilarity(taxonToken(4'u32), taxonToken(5'u32), taxonomy) - 0.5) < 1e-9

    let a = parseKrakenSignature("4:1 0:1")
    let b = parseKrakenSignature("5:1 0:1")
    check abs(taxonomicPositionalSimilarity(a, b, taxonomy) - 0.75) < 1e-9

  test "confidence scoring follows Kraken2 denominator rules":
    var tree = initTable[uint32, seq[uint32]]()
    tree[1'u32] = @[561'u32]
    tree[561'u32] = @[562'u32]

    let taxonomy = buildTaxonomyIndex(tree)
    let sig = parseKrakenSignature("562:13 561:4 A:31 0:1 562:3")

    check sig.queriedKmerCount() == 21
    check abs(confidenceScore(sig, 562'u32, taxonomy) - (16.0 / 21.0)) < 1e-9
    check abs(confidenceScore(sig, 561'u32, taxonomy) - (20.0 / 21.0)) < 1e-9
    check taxIdAtConfidence(sig, 562'u32, taxonomy, 16.0 / 21.0) == 562'u32
    check taxIdAtConfidence(sig, 562'u32, taxonomy, (16.0 / 21.0) + 1e-6) == 561'u32
    check taxIdAtConfidence(sig, 562'u32, taxonomy, (20.0 / 21.0) + 1e-6) == 0'u32
    check taxIdAtConfidence(sig, taxonomy, 0.8) == 561'u32

  test "confidence scoring returns unclassified for no queried k-mers":
    var tree = initTable[uint32, seq[uint32]]()
    tree[1'u32] = @[562'u32]

    let taxonomy = buildTaxonomyIndex(tree)
    let sig = parseKrakenSignature("A:31")

    check sig.queriedKmerCount() == 0
    check confidenceScore(sig, 562'u32, taxonomy) == 0.0
    check taxIdAtConfidence(sig, 562'u32, taxonomy, 0.0) == 0'u32

# ---------------------------------------------------------------------------
# Integration / filtering tests (requires compiled binary)
# ---------------------------------------------------------------------------

suite "filtering (paired-end)":

  setup:
    if barena == "":
      skip()

  let tmp = getTempDir() / "barena_tests"
  createDir(tmp)

  test "taxon 3, no children: 1 read":
    let tmpl = tmp / "t3_nochildren_#.fq"
    let (_, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " -t 3 -c " & tmpl.quoteShell
    )
    check code == 0
    check countReads(tmp / "t3_nochildren_R1.fq") == 1

  test "taxon 3, with children: 5 reads":
    let tmpl = tmp / "t3_children_#.fq"
    let (_, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " -t 3 --child -d " & InspectDb &
      " -c " & tmpl.quoteShell
    )
    check code == 0
    check countReads(tmp / "t3_children_R1.fq") == 5

  test "taxon 3, children + --include-unclassified: 8 reads":
    let tmpl = tmp / "t3_unclass_#.fq"
    let (_, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " -t 3 --child --include-unclassified -d " & InspectDb &
      " -c " & tmpl.quoteShell
    )
    check code == 0
    check countReads(tmp / "t3_unclass_R1.fq") == 8

  test "taxon 46230, no children: 3 reads":
    let tmpl = tmp / "t46230_nochildren_#.fq"
    let (_, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " -t 46230 -c " & tmpl.quoteShell
    )
    check code == 0
    check countReads(tmp / "t46230_nochildren_R1.fq") == 3

  test "taxon 46230, with children: still 3 reads":
    let tmpl = tmp / "t46230_children_#.fq"
    let (_, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " -t 46230 --child -d " & InspectDb &
      " -c " & tmpl.quoteShell
    )
    check code == 0
    check countReads(tmp / "t46230_children_R1.fq") == 3

  test "split mode applies stem, explicit extension, template, and basename-dot naming":
    let outdir = tmp / "split_names"
    freshDir(outdir)
    let splitFile = tmp / "split_names.tsv"
    writeFile(splitFile,
      "TaxID\tOutputFile\n" &
      "3\tStem\n" &
      "46230\texplicit.fq\n" &
      "44098\tnested.dir/NoDot\n" &
      "0\tunclass_#.fq\n")

    let (err, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " --split " & splitFile.quoteShell &
      " --outdir " & outdir.quoteShell
    )
    checkpoint err
    check code == 0
    check countReads(outdir / "Stem_R1.fastq") == 1
    check countReads(outdir / "Stem_R2.fastq") == 1
    check countReads(outdir / "explicit_R1.fq") == 3
    check countReads(outdir / "explicit_R2.fq") == 3
    check countReads(outdir / "nested.dir" / "NoDot_R1.fastq") == 1
    check countReads(outdir / "nested.dir" / "NoDot_R2.fastq") == 1
    check countReads(outdir / "unclass_R1.fq") == 3
    check countReads(outdir / "unclass_R2.fq") == 3

  test "split child mode writes overlapping matches to all groups":
    let outdir = tmp / "split_child"
    freshDir(outdir)
    let splitFile = tmp / "split_child.tsv"
    writeFile(splitFile,
      "TaxID\tOutputFile\n" &
      "3\tRoot\n" &
      "46230\tLeaf\n")

    let (err, code) = run(
      "-1 " & R1 & " -2 " & R2 &
      " -k " & KrakenRaw &
      " --split " & splitFile.quoteShell &
      " --child -d " & InspectDb &
      " --outdir " & outdir.quoteShell
    )
    checkpoint err
    check code == 0
    check countReads(outdir / "Root_R1.fastq") == 5
    check countReads(outdir / "Leaf_R1.fastq") == 3

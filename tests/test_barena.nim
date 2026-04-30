## Test suite for barena.
## Run with: nimble test   (or: nim c -r tests/test_barena.nim)

import std/[unittest, os, osproc, tables, sets]
import "../src/taxonomy"
import "../src/krakenio"

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

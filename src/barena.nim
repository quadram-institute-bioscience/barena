## barena — filter FASTQ reads by taxon using Kraken2 classification output.
##
## Reads Kraken2 raw output (stdin or file) and a FASTQ file (or pair) in
## lock-step (positional order is assumed to match). Emits reads that match
## the target taxon(s) and/or those that don't.

import std/[sets, strutils, os, sequtils]
import argparse
import readfx
import taxonomy, krakenio, outfiles

# ---------------------------------------------------------------------------
# Output path resolution
# ---------------------------------------------------------------------------

type OutPaths = object
  classR1, classR2, unclassR1, unclassR2: string

proc applyTemplate(tmpl, tag: string): string =
  if '#' in tmpl: tmpl.replace("#", tag) else: tmpl

proc buildOutPaths(outdir, classOut, classOutR2, unclassOut, unclassOutR2: string,
                   paired: bool): OutPaths =
  if outdir.len > 0:
    createDir(outdir)
    result.classR1   = outdir / "classified_R1.fastq.gz"
    result.unclassR1 = outdir / "unclassified_R1.fastq.gz"
    if paired:
      result.classR2   = outdir / "classified_R2.fastq.gz"
      result.unclassR2 = outdir / "unclassified_R2.fastq.gz"
    return

  if classOut.len > 0:
    result.classR1 = applyTemplate(classOut, "R1")
    if paired:
      if '#' in classOut:       result.classR2 = classOut.replace("#", "R2")
      elif classOutR2.len > 0:  result.classR2 = classOutR2
  elif classOutR2.len > 0 and paired:
    result.classR2 = classOutR2

  if unclassOut.len > 0:
    result.unclassR1 = applyTemplate(unclassOut, "R1")
    if paired:
      if '#' in unclassOut:      result.unclassR2 = unclassOut.replace("#", "R2")
      elif unclassOutR2.len > 0: result.unclassR2 = unclassOutR2
  elif unclassOutR2.len > 0 and paired:
    result.unclassR2 = unclassOutR2

proc hasAnyOutput(p: OutPaths): bool =
  p.classR1.len > 0 or p.classR2.len > 0 or
  p.unclassR1.len > 0 or p.unclassR2.len > 0

# ---------------------------------------------------------------------------
# Filtering logic
# ---------------------------------------------------------------------------

proc isMatch(classified: bool, taxId: uint32,
             matchIds: HashSet[uint32], includeUnclassified: bool): bool {.inline.} =
  if classified:
    taxId in matchIds
  else:
    includeUnclassified

proc checkStrict(krakenLine, fastqName: string) {.inline.} =
  let kName = stripPairSuffix(parseKrakenReadName(krakenLine))
  let fName = stripPairSuffix(fastqName)
  if kName != fName:
    stderr.writeLine "barena: read name mismatch\n  kraken: ", kName, "\n  fastq:  ", fName
    quit(1)

const ProgressInterval = 20_000'u64

proc reportProgress(nTotal, nMatch: uint64) {.inline.} =
  let pct = if nTotal > 0: nMatch.float / nTotal.float * 100.0 else: 0.0
  stderr.write "\r  ", nMatch, " / ", nTotal, " reads matched (",
               formatFloat(pct, ffDecimal, 1), "%)"
  stderr.flushFile()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

when isMainModule:
  let p = newParser("barena"):
    help("Filter FASTQ reads by taxon using Kraken2 classification output")

    # Input
    option("-1", "--r1",         help="FASTQ file (single-end or R1 of pair)", required=true)
    option("-2", "--r2",         help="FASTQ file (R2 of pair; enables paired-end mode)")
    option("-k", "--kraken",     help="Kraken2 raw output file (default: stdin)")
    option("-d", "--inspect-db", help="kraken2-inspect output (required with --child)")

    # Output  (mutually exclusive with -o; -o always produces gzipped files)
    option("-o", "--outdir",         help="Output directory (gzipped .fastq.gz, overrides -c/-u)")
    option("-c", "--class-out",      help="Classified R1 output path (use # as R1/R2 placeholder)")
    option("-C", "--class-out-r2",   help="Classified R2 output path")
    option("-u", "--unclass-out",    help="Unclassified R1 output path (use # as R1/R2 placeholder)")
    option("-U", "--unclass-out-r2", help="Unclassified R2 output path")

    # Filter
    option("-t", "--taxon", help="Taxonomy ID to filter (repeat for multiple IDs)",
           multiple=true)
    flag("--child",
         help="Include all descendants of given taxids (requires -d/--inspect-db)")
    flag("--include-unclassified",
         help="Route kraken-unclassified reads (taxId=0) to the classified output")
    flag("--strict",
         help="Abort if read order differs between kraken output and FASTQ")
    flag("--verbose",
         help="Print taxonomy stats and progress every 20,000 reads to stderr")

  try:
    let opts = p.parse()
    let paired = opts.r2.len > 0

    # ---- validate output specification ----
    let paths = buildOutPaths(
      opts.outdir, opts.class_out, opts.class_out_r2,
      opts.unclass_out, opts.unclass_out_r2, paired
    )
    if not paths.hasAnyOutput():
      stderr.writeLine "barena: no output specified — use -o, -c/-C, or -u/-U"
      quit(1)

    # ---- build taxid match set ----
    if opts.taxon.len == 0 and not opts.include_unclassified:
      stderr.writeLine "barena: no --taxon given and --include-unclassified not set; nothing would match"
      quit(1)

    var matchIds: HashSet[uint32]
    if opts.child:
      if opts.inspect_db.len == 0:
        stderr.writeLine "barena: --child requires --inspect-db (-d)"
        quit(1)
      var nodeCount: int
      let tree = parseInspectDb(opts.inspect_db, nodeCount)
      matchIds = getDescendants(tree, opts.taxon.mapIt(parseUInt(it).uint32))
      if opts.verbose:
        stderr.writeLine "Taxonomy parsed: ", nodeCount, " nodes  |  match set: ", matchIds.len, " taxids"
    else:
      matchIds = initHashSet[uint32]()
      for t in opts.taxon:
        matchIds.incl(parseUInt(t).uint32)

    # ---- open outputs ----
    var classR1   = openOut(paths.classR1)
    var classR2   = openOut(paths.classR2)
    var unclassR1 = openOut(paths.unclassR1)
    var unclassR2 = openOut(paths.unclassR2)
    defer:
      classR1.close(); classR2.close()
      unclassR1.close(); unclassR2.close()

    # ---- open kraken stream (plain or gzip, file or stdin) ----
    var krakenStream = xopen[GzFile](if opts.kraken.len > 0: opts.kraken else: "-")
    defer: krakenStream.close()

    # ---- shared state ----
    var krakenLine: string
    var taxId: uint32
    var writeBuf = newStringOfCap(1024)
    var nTotal, nMatch, nUnmatch: uint64

    # ---- main filter loop ----
    if paired:
      var r1, r2: FQRecord
      var f1 = xopen[GzFile](opts.r1)
      var f2 = xopen[GzFile](opts.r2)
      defer: f1.close(); f2.close()

      while f1.readFastx(r1):
        if not f2.readFastx(r2):
          stderr.writeLine "barena: R2 file ended before R1"
          quit(1)
        if not krakenStream.readLine(krakenLine):
          stderr.writeLine "barena: kraken output ended before FASTQ"
          quit(1)
        if opts.strict:
          checkStrict(krakenLine, r1.name)

        let classified = parseKrakenLine(krakenLine, taxId)
        inc nTotal
        if isMatch(classified, taxId, matchIds, opts.include_unclassified):
          inc nMatch
          classR1.writeFq(r1, writeBuf)
          classR2.writeFq(r2, writeBuf)
        else:
          inc nUnmatch
          unclassR1.writeFq(r1, writeBuf)
          unclassR2.writeFq(r2, writeBuf)
        if opts.verbose and nTotal mod ProgressInterval == 0:
          reportProgress(nTotal, nMatch)

    else:
      var r1: FQRecord
      var f1 = xopen[GzFile](opts.r1)
      defer: f1.close()

      while f1.readFastx(r1):
        if not krakenStream.readLine(krakenLine):
          stderr.writeLine "barena: kraken output ended before FASTQ"
          quit(1)
        if opts.strict:
          checkStrict(krakenLine, r1.name)

        let classified = parseKrakenLine(krakenLine, taxId)
        inc nTotal
        if isMatch(classified, taxId, matchIds, opts.include_unclassified):
          inc nMatch
          classR1.writeFq(r1, writeBuf)
        else:
          inc nUnmatch
          unclassR1.writeFq(r1, writeBuf)
        if opts.verbose and nTotal mod ProgressInterval == 0:
          reportProgress(nTotal, nMatch)

    if opts.verbose:
      stderr.writeLine ""  # end the \r progress line
    stderr.writeLine "Processed: ", nTotal, "  matched: ", nMatch,
                     "  unmatched: ", nUnmatch

  except ShortCircuit as e:
    if e.flag == "argparse_help":
      echo p.help
      quit(0)
  except UsageError as e:
    stderr.writeLine "barena: ", e.msg
    quit(1)

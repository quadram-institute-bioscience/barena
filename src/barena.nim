## barena — filter FASTQ reads by taxon using Kraken2 classification output.
##
## Reads Kraken2 raw output (stdin or file) and a FASTQ file (or pair) in
## lock-step (positional order is assumed to match). Emits reads that match
## the target taxon(s) and/or those that don't.

import std/[sets, strutils, os, sequtils, tables]
import argparse
import readfx
import taxonomy, krakenio, outfiles

const NimblePkgVersion {.strdefine.} = "0.0.0"

# ---------------------------------------------------------------------------
# Output path resolution
# ---------------------------------------------------------------------------

type OutPaths = object
  classR1, classR2, unclassR1, unclassR2: string

type SplitTarget = object
  taxId: uint32
  label: string
  r1Path, r2Path: string

const
  DefaultSplitExt = ".fastq"
  MaxSplitOutputFiles = 240

proc die(msg: string) {.noreturn.} =
  stderr.writeLine "barena: ", msg
  quit(1)

proc applyTemplate(tmpl, tag: string): string =
  if '#' in tmpl: tmpl.replace("#", tag) else: tmpl

proc hasDotInBasename(path: string): bool =
  '.' in path.extractFilename

proc qualifyOutputPath(outdir, path: string): string =
  if outdir.len > 0 and not path.isAbsolute:
    outdir / path
  else:
    path

proc splitKnownExt(path: string): tuple[stem, ext: string] =
  let lower = path.toLowerAscii()
  for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
    if lower.endsWith(ext):
      return (path[0 ..< path.len - ext.len], path[path.len - ext.len .. ^1])
  let parts = splitFile(path)
  result.stem = if parts.dir.len > 0: parts.dir / parts.name else: parts.name
  result.ext = parts.ext

proc addPairTag(path, tag: string): string =
  let (stem, ext) = splitKnownExt(path)
  stem & "_" & tag & ext

proc buildSplitPaths(outdir, outputSpec: string, paired: bool,
                     lineNo: int): tuple[r1, r2: string] =
  let trimmed = outputSpec.strip()
  if trimmed.len == 0:
    die("--split line " & $lineNo & " has an empty OutputFile")

  let base = trimmed.extractFilename
  if base.len == 0:
    die("--split line " & $lineNo & " has an invalid OutputFile")

  if hasDotInBasename(trimmed):
    let path = qualifyOutputPath(outdir, trimmed)
    if '#' in path:
      result.r1 = applyTemplate(path, "R1")
      if paired:
        result.r2 = applyTemplate(path, "R2")
    else:
      if paired:
        result.r1 = addPairTag(path, "R1")
        result.r2 = addPairTag(path, "R2")
      else:
        result.r1 = path
  else:
    if '#' in trimmed:
      die("--split line " & $lineNo & " uses # without an explicit extension")
    let stem = qualifyOutputPath(outdir, trimmed)
    if paired:
      result.r1 = stem & "_R1" & DefaultSplitExt
      result.r2 = stem & "_R2" & DefaultSplitExt
    else:
      result.r1 = stem & DefaultSplitExt

proc parseTaxId(raw, context: string): uint32 =
  try:
    result = parseUInt(raw).uint32
  except ValueError:
    die(context & " is not a valid TaxID: " & raw)

proc readSplitTargets(path, outdir: string, paired: bool): seq[SplitTarget] =
  var firstContent = true
  var seenOutputs = initHashSet[string]()

  var lineNo = 0
  for line in path.lines:
    inc lineNo
    let stripped = line.strip()
    if stripped.len == 0 or stripped[0] == '#':
      continue

    let parts = stripped.splitWhitespace()
    if firstContent:
      firstContent = false
      if parts.len == 2 and parts[0].cmpIgnoreCase("TaxID") == 0 and
         parts[1].cmpIgnoreCase("OutputFile") == 0:
        continue

    if parts.len != 2:
      die("--split line " & $lineNo & " must contain exactly two columns: TaxID OutputFile")

    let paths = buildSplitPaths(outdir, parts[1], paired, lineNo)
    for outputPath in [paths.r1, paths.r2]:
      if outputPath.len == 0:
        continue
      if seenOutputs.containsOrIncl(outputPath):
        die("--split output path is duplicated: " & outputPath)

    result.add SplitTarget(
      taxId: parseTaxId(parts[0], "--split line " & $lineNo),
      label: parts[1],
      r1Path: paths.r1,
      r2Path: paths.r2
    )

  if result.len == 0:
    die("--split file contains no targets: " & path)

proc checkSplitOutputCount(targets: seq[SplitTarget], paired: bool,
                           unclass: OutPaths) =
  var outputCount = targets.len * (if paired: 2 else: 1)
  if unclass.unclassR1.len > 0: inc outputCount
  if unclass.unclassR2.len > 0: inc outputCount
  if outputCount > MaxSplitOutputFiles:
    die("--split would open " & $outputCount & " output files; limit is " &
        $MaxSplitOutputFiles)

proc checkUnclassOutputCollisions(targets: seq[SplitTarget], unclass: OutPaths) =
  for target in targets:
    if unclass.unclassR1.len > 0 and
       (unclass.unclassR1 == target.r1Path or unclass.unclassR1 == target.r2Path):
      die("--unclass-out duplicates a --split output path: " & unclass.unclassR1)
    if unclass.unclassR2.len > 0 and
       (unclass.unclassR2 == target.r1Path or unclass.unclassR2 == target.r2Path):
      die("--unclass-out-r2 duplicates a --split output path: " & unclass.unclassR2)

proc addRoute(routes: var Table[uint32, seq[int]], taxId: uint32, index: int) =
  routes.mgetOrPut(taxId, @[]).add(index)

proc buildUnclassPaths(outdir, unclassOut, unclassOutR2: string,
                       paired: bool): OutPaths =
  if unclassOut.len > 0:
    result.unclassR1 = qualifyOutputPath(outdir, applyTemplate(unclassOut, "R1"))
    if paired:
      if '#' in unclassOut:
        result.unclassR2 = qualifyOutputPath(outdir, unclassOut.replace("#", "R2"))
      elif unclassOutR2.len > 0:
        result.unclassR2 = qualifyOutputPath(outdir, unclassOutR2)
  elif unclassOutR2.len > 0 and paired:
    result.unclassR2 = qualifyOutputPath(outdir, unclassOutR2)

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

    # Output  (normal mode: -o produces gzipped defaults; split mode: -o is a root)
    option("-o", "--outdir",         help="Output directory (normal mode uses default .fastq.gz names; split mode roots OutputFile paths)")
    option("-c", "--class-out",      help="Classified R1 output path (use # as R1/R2 placeholder)")
    option("-C", "--class-out-r2",   help="Classified R2 output path")
    option("-u", "--unclass-out",    help="Unclassified R1 output path (use # as R1/R2 placeholder)")
    option("-U", "--unclass-out-r2", help="Unclassified R2 output path")

    # Filter
    option("-t", "--taxon", help="Taxonomy ID to filter (repeat for multiple IDs)",
           multiple=true)
    option("--split",
           help="TSV/whitespace file with TaxID and OutputFile columns for per-taxon outputs")
    flag("--child",
         help="Include all descendants of given taxids (requires -d/--inspect-db)")
    flag("--include-unclassified",
         help="Route kraken-unclassified reads (taxId=0) to the classified output")
    flag("--strict",
         help="Abort if read order differs between kraken output and FASTQ")
    flag("--verbose",
         help="Print taxonomy stats and progress every 20,000 reads to stderr")

  if "--version" in commandLineParams():
    echo "barena v", NimblePkgVersion
    quit(0)
  try:
    let opts = p.parse()
    let paired = opts.r2.len > 0
    let splitMode = opts.split.len > 0

    # ---- validate output specification ----
    var paths: OutPaths
    var splitTargets: seq[SplitTarget]
    if splitMode:
      if opts.taxon.len > 0:
        die("--split cannot be combined with --taxon")
      if opts.class_out.len > 0 or opts.class_out_r2.len > 0:
        die("--split cannot be combined with --class-out/--class-out-r2")
      if opts.include_unclassified:
        die("--include-unclassified cannot be combined with --split; add a TaxID 0 row instead")
      splitTargets = readSplitTargets(opts.split, opts.outdir, paired)
      paths = buildUnclassPaths(opts.outdir, opts.unclass_out, opts.unclass_out_r2, paired)
      checkSplitOutputCount(splitTargets, paired, paths)
      checkUnclassOutputCollisions(splitTargets, paths)
    else:
      paths = buildOutPaths(
        opts.outdir, opts.class_out, opts.class_out_r2,
        opts.unclass_out, opts.unclass_out_r2, paired
      )
      if not paths.hasAnyOutput():
        stderr.writeLine "barena: no output specified — use -o, -c/-C, or -u/-U"
        quit(1)

    # ---- build taxid match set ----
    if not splitMode and opts.taxon.len == 0 and not opts.include_unclassified:
      stderr.writeLine "barena: no --taxon given and --include-unclassified not set; nothing would match"
      quit(1)

    var matchIds: HashSet[uint32]
    var splitRoutes = initTable[uint32, seq[int]]()
    if splitMode and opts.child:
      if opts.inspect_db.len == 0:
        stderr.writeLine "barena: --child requires --inspect-db (-d)"
        quit(1)
      var nodeCount: int
      let tree = parseInspectDb(opts.inspect_db, nodeCount)
      for i, target in splitTargets:
        for id in getDescendants(tree, @[target.taxId]):
          splitRoutes.addRoute(id, i)
      if opts.verbose:
        stderr.writeLine "Taxonomy parsed: ", nodeCount, " nodes  |  split targets: ",
                         splitTargets.len, "  |  routed taxids: ", splitRoutes.len
    elif splitMode:
      for i, target in splitTargets:
        splitRoutes.addRoute(target.taxId, i)
    elif opts.child:
      if opts.inspect_db.len == 0:
        stderr.writeLine "barena: --child requires --inspect-db (-d)"
        quit(1)
      var nodeCount: int
      let tree = parseInspectDb(opts.inspect_db, nodeCount)
      matchIds = getDescendants(tree, opts.taxon.mapIt(parseTaxId(it, "--taxon")))
      if opts.verbose:
        stderr.writeLine "Taxonomy parsed: ", nodeCount, " nodes  |  match set: ", matchIds.len, " taxids"
    else:
      matchIds = initHashSet[uint32]()
      for t in opts.taxon:
        matchIds.incl(parseTaxId(t, "--taxon"))

    # ---- open outputs ----
    var classR1   = openOut(paths.classR1)
    var classR2   = openOut(paths.classR2)
    var unclassR1 = openOut(paths.unclassR1)
    var unclassR2 = openOut(paths.unclassR2)
    var splitR1: seq[OutFile]
    var splitR2: seq[OutFile]
    var splitCounts = newSeq[uint64](splitTargets.len)
    if splitMode:
      for target in splitTargets:
        splitR1.add(openOut(target.r1Path))
        if paired:
          splitR2.add(openOut(target.r2Path))
    defer:
      classR1.close(); classR2.close()
      unclassR1.close(); unclassR2.close()
      for i in 0 ..< splitR1.len:
        splitR1[i].close()
      for i in 0 ..< splitR2.len:
        splitR2[i].close()

    # ---- open kraken stream (plain or gzip, file or stdin) ----
    var krakenStream = xopen[GzFile](if opts.kraken.len > 0: opts.kraken else: "-")
    defer: krakenStream.close()

    # ---- shared state ----
    var krakenLine: string
    var taxId: uint32
    var writeBuf = newStringOfCap(1024)
    var nTotal, nMatch, nUnmatch, nSplitWrites: uint64

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
        if splitMode:
          if taxId in splitRoutes:
            inc nMatch
            for i in splitRoutes[taxId]:
              inc splitCounts[i]
              inc nSplitWrites
              splitR1[i].writeFq(r1, writeBuf)
              splitR2[i].writeFq(r2, writeBuf)
          else:
            inc nUnmatch
            unclassR1.writeFq(r1, writeBuf)
            unclassR2.writeFq(r2, writeBuf)
        else:
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
        if splitMode:
          if taxId in splitRoutes:
            inc nMatch
            for i in splitRoutes[taxId]:
              inc splitCounts[i]
              inc nSplitWrites
              splitR1[i].writeFq(r1, writeBuf)
          else:
            inc nUnmatch
            unclassR1.writeFq(r1, writeBuf)
        else:
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
      if splitMode:
        for i, target in splitTargets:
          stderr.writeLine "  ", target.label, ": ", splitCounts[i], " reads"
    stderr.writeLine "Processed: ", nTotal, "  matched: ", nMatch,
                     "  unmatched: ", nUnmatch,
                     (if splitMode: "  split-writes: " & $nSplitWrites else: "")

  except ShortCircuit as e:
    if e.flag == "argparse_help":
      echo p.help
      quit(0)
  except UsageError as e:
    stderr.writeLine "barena: ", e.msg
    quit(1)

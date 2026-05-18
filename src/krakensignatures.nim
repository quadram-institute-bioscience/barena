## Parse, store, and compare Kraken2 k-mer signatures.
##
## Kraken2's final raw-output field is a run-length encoded stream of
## taxonomic k-mer calls, optionally split by "|:|" for paired-end reads:
##
##   3:80 0:5 46230:20 |:| 1:11 3:96 A:2
##
## This module keeps that representation compact while exposing comparison
## routines that operate as if each run had been expanded position-by-position.

import std/[algorithm, hashes, sets, strutils, tables]

import taxonomy

type
  KrakenSignatureError* = object of ValueError
    ## Raised when a Kraken2 k-mer signature cannot be parsed.

  KmerTokenKind* = enum
    ktkUnmapped,     ## Kraken token "0": no database hit for this k-mer.
    ktkAmbiguous,    ## Kraken token "A": k-mer contains ambiguous bases.
    ktkTaxon         ## A positive NCBI/Kraken taxid.

  KmerToken* = object
    ## A single k-mer state in a Kraken2 signature.
    case kind*: KmerTokenKind
    of ktkTaxon:
      taxId*: uint32
    of ktkUnmapped, ktkAmbiguous:
      discard

  KmerRun* = object
    ## One run in Kraken2's run-length encoded signature.
    token*: KmerToken
    count*: int

  KmerMate* = seq[KmerRun]

  KrakenSignature* = object
    ## One Kraken2 signature. Single-end reads have one mate; paired-end reads
    ## normally have two mates, in read order.
    mates*: seq[KmerMate]

  TaxonomyIndex* = object
    ## Parent/depth view of taxonomy.nim's parent -> children TaxTree.
    parent*: Table[uint32, uint32]
    depth*: Table[uint32, int]
    roots*: seq[uint32]

  TokenSimilarityProc* = proc(a, b: KmerToken): float {.closure.}
    ## Scores two aligned k-mer tokens. Return values are expected in [0, 1].

  SignatureCountKey* = object
    ## A key for composition-based comparisons.
    mate*: int
    token*: KmerToken

proc raiseSignatureError(msg: string) {.noReturn.} =
  raise newException(KrakenSignatureError, msg)

proc unmappedToken*(): KmerToken {.inline.} =
  KmerToken(kind: ktkUnmapped)

proc ambiguousToken*(): KmerToken {.inline.} =
  KmerToken(kind: ktkAmbiguous)

proc taxonToken*(taxId: uint32): KmerToken {.inline.} =
  ## Construct a taxon token. Taxid 0 is normalized to Kraken's unmapped token.
  if taxId == 0'u32:
    unmappedToken()
  else:
    KmerToken(kind: ktkTaxon, taxId: taxId)

proc `==`*(a, b: KmerToken): bool =
  if a.kind != b.kind:
    return false

  case a.kind
  of ktkTaxon:
    a.taxId == b.taxId
  of ktkUnmapped, ktkAmbiguous:
    true

proc `$`*(token: KmerToken): string =
  case token.kind
  of ktkUnmapped:
    "0"
  of ktkAmbiguous:
    "A"
  of ktkTaxon:
    $token.taxId

proc hash*(token: KmerToken): Hash =
  var h = hash(ord(token.kind))
  if token.kind == ktkTaxon:
    h = h !& hash(token.taxId)
  !$h

proc hash*(key: SignatureCountKey): Hash =
  var h = hash(key.mate)
  h = h !& hash(key.token)
  !$h

proc isDigits(s: string): bool =
  if s.len == 0:
    return false
  for c in s:
    if c < '0' or c > '9':
      return false
  true

proc parseUint32Strict(s, what: string): uint32 =
  if not isDigits(s):
    raiseSignatureError("invalid " & what & " '" & s & "'")

  let value =
    try:
      parseUInt(s)
    except ValueError:
      raiseSignatureError("invalid " & what & " '" & s & "'")

  if value > uint(uint32.high):
    raiseSignatureError(what & " is too large for uint32: " & s)

  value.uint32

proc parsePositiveCount(s: string): int =
  if not isDigits(s):
    raiseSignatureError("invalid k-mer count '" & s & "'")

  let value =
    try:
      parseUInt(s)
    except ValueError:
      raiseSignatureError("invalid k-mer count '" & s & "'")

  if value == 0'u:
    raiseSignatureError("k-mer run count must be positive")
  if value > uint(int.high):
    raiseSignatureError("k-mer run count is too large: " & s)

  value.int

proc parseKmerToken*(s: string): KmerToken =
  ## Parse the taxid part of a Kraken2 run token.
  if s == "A":
    return ambiguousToken()

  let taxId = parseUint32Strict(s, "taxid")
  taxonToken(taxId)

proc parseKmerRun*(s: string): KmerRun =
  ## Parse one Kraken2 run token, for example "562:12" or "A:3".
  let colon = s.find(':')
  if colon <= 0 or colon == s.len - 1 or s.find(':', colon + 1) >= 0:
    raiseSignatureError("invalid k-mer run '" & s & "'")

  result.token = parseKmerToken(s[0 ..< colon])
  result.count = parsePositiveCount(s[colon + 1 .. ^1])

proc addRun*(mate: var KmerMate, run: KmerRun) =
  ## Add a run, merging adjacent equal tokens so signatures stay normalized.
  if run.count <= 0:
    raiseSignatureError("k-mer run count must be positive")

  if mate.len > 0 and mate[^1].token == run.token:
    if mate[^1].count > int.high - run.count:
      raiseSignatureError("merged k-mer run count overflow")
    mate[^1].count += run.count
  else:
    mate.add(run)

proc parseKrakenSignature*(field: string): KrakenSignature =
  ## Parse the final field of Kraken2 raw output.
  ##
  ## The "|:|" paired-end separator may be surrounded by spaces or not.
  let normalized = field.strip().replace("|:|", " |:| ")
  if normalized.len == 0:
    raiseSignatureError("empty Kraken2 signature")

  var mate: KmerMate = @[]
  var sawSeparator = false

  for token in normalized.splitWhitespace():
    if token == "|:|":
      if mate.len == 0:
        raiseSignatureError("empty mate before paired-end separator")
      result.mates.add(mate)
      mate = @[]
      sawSeparator = true
    else:
      mate.addRun(parseKmerRun(token))

  if mate.len == 0:
    if sawSeparator:
      raiseSignatureError("empty mate after paired-end separator")
    raiseSignatureError("empty Kraken2 signature")

  result.mates.add(mate)

proc parseKrakenLineSignature*(line: string): KrakenSignature =
  ## Extract and parse field 5 from a full Kraken2 raw-output line.
  var tabs = 0
  var start = -1

  for i, c in line:
    if c == '\t':
      inc tabs
      if tabs == 4:
        start = i + 1
        break

  if start < 0:
    raiseSignatureError("Kraken2 line has fewer than 5 tab-separated fields")

  parseKrakenSignature(line[start .. ^1])

proc tryParseKrakenSignature*(field: string, signature: var KrakenSignature): bool =
  ## Parse a signature without raising on invalid input.
  try:
    signature = parseKrakenSignature(field)
    true
  except ValueError:
    signature = KrakenSignature()
    false

proc toKrakenString*(run: KmerRun): string =
  $run.token & ":" & $run.count

proc toKrakenString*(mate: KmerMate): string =
  var parts: seq[string] = @[]
  for run in mate:
    parts.add(run.toKrakenString())
  parts.join(" ")

proc toKrakenString*(signature: KrakenSignature): string =
  var parts: seq[string] = @[]
  for mate in signature.mates:
    parts.add(mate.toKrakenString())
  parts.join(" |:| ")

proc `$`*(signature: KrakenSignature): string =
  signature.toKrakenString()

proc kmerCount*(mate: KmerMate): int =
  for run in mate:
    result += run.count

proc kmerCount*(signature: KrakenSignature): int =
  for mate in signature.mates:
    result += mate.kmerCount()

proc queriedKmerCount*(signature: KrakenSignature): int =
  ## Count k-mers queried against the database for Kraken2 confidence scoring.
  ## Ambiguous A runs are excluded; unmapped 0 runs are included.
  for mate in signature.mates:
    for run in mate:
      if run.token.kind != ktkAmbiguous:
        result += run.count

proc mateCount*(signature: KrakenSignature): int {.inline.} =
  signature.mates.len

proc isPaired*(signature: KrakenSignature): bool {.inline.} =
  signature.mates.len == 2

proc normalized*(signature: KrakenSignature): KrakenSignature =
  ## Return a copy with adjacent equal runs merged in each mate.
  for mate in signature.mates:
    var outMate: KmerMate = @[]
    for run in mate:
      outMate.addRun(run)
    result.mates.add(outMate)

proc buildTaxonomyIndex*(tree: TaxTree): TaxonomyIndex =
  ## Build parent and depth maps from taxonomy.nim's parent -> children tree.
  var allIds = initHashSet[uint32]()
  var childIds = initHashSet[uint32]()

  result.parent = initTable[uint32, uint32]()
  result.depth = initTable[uint32, int]()
  result.roots = @[]

  for parent, children in tree:
    allIds.incl(parent)
    for child in children:
      allIds.incl(child)
      childIds.incl(child)
      result.parent[child] = parent

  for id in allIds:
    if id notin childIds:
      result.roots.add(id)

  result.roots.sort()

  var queue: seq[tuple[id: uint32, depth: int]] = @[]
  for root in result.roots:
    queue.add((root, 0))

  var head = 0
  while head < queue.len:
    let (id, depth) = queue[head]
    inc head

    if id in result.depth and result.depth[id] <= depth:
      continue

    result.depth[id] = depth

    if id in tree:
      for child in tree[id]:
        queue.add((child, depth + 1))

  for id in allIds:
    if id notin result.depth:
      result.depth[id] = 0

proc hasTaxId*(taxonomy: TaxonomyIndex, taxId: uint32): bool {.inline.} =
  taxId in taxonomy.depth

proc parentOf*(taxonomy: TaxonomyIndex, taxId: uint32, parent: var uint32): bool =
  if taxId in taxonomy.parent:
    parent = taxonomy.parent[taxId]
    true
  else:
    false

proc depthOf*(taxonomy: TaxonomyIndex, taxId: uint32): int =
  taxonomy.depth.getOrDefault(taxId, -1)

proc validateConfidenceThreshold(threshold: float) =
  if threshold != threshold or threshold < 0.0 or threshold > 1.0:
    raiseSignatureError("confidence threshold must be in [0, 1]")

proc addAncestorCount(
    counts: var Table[uint32, int],
    taxonomy: TaxonomyIndex,
    taxId: uint32,
    count: int
  ) =
  var current = taxId

  while true:
    counts[current] = counts.getOrDefault(current, 0) + count

    var parent: uint32
    if not taxonomy.parentOf(current, parent) or parent == current:
      break

    current = parent

proc cladeKmerCounts(
    signature: KrakenSignature,
    taxonomy: TaxonomyIndex
  ): Table[uint32, int] =
  result = initTable[uint32, int]()

  for mate in signature.mates:
    for run in mate:
      if run.token.kind == ktkTaxon:
        result.addAncestorCount(taxonomy, run.token.taxId, run.count)

proc confidenceScore*(
    signature: KrakenSignature,
    taxId: uint32,
    taxonomy: TaxonomyIndex
  ): float =
  ## Return Kraken2's confidence score C/Q for a candidate taxid.
  ##
  ## C is the number of taxon k-mers assigned to the candidate's clade. Q is
  ## all queried k-mers, including unmapped 0 runs and excluding ambiguous A
  ## runs.
  if taxId == 0'u32:
    return 0.0

  let q = signature.queriedKmerCount()
  if q == 0:
    return 0.0

  let counts = signature.cladeKmerCounts(taxonomy)
  counts.getOrDefault(taxId, 0).float / q.float

proc deepestTaxId*(signature: KrakenSignature, taxonomy: TaxonomyIndex): uint32 =
  ## Infer the deepest taxid present in the signature.
  ##
  ## Prefer passing Kraken2's original column-3 taxid to taxIdAtConfidence when
  ## it is available; this helper is a fallback for signature-only callers.
  var bestDepth = -1

  for mate in signature.mates:
    for run in mate:
      if run.token.kind == ktkTaxon:
        let depth = taxonomy.depthOf(run.token.taxId)
        if result == 0'u32 or depth > bestDepth:
          result = run.token.taxId
          bestDepth = depth

proc taxIdAtConfidence*(
    signature: KrakenSignature,
    originalTaxId: uint32,
    taxonomy: TaxonomyIndex,
    threshold: float
  ): uint32 =
  ## Return the deepest ancestor of originalTaxId whose Kraken2-style
  ## confidence score meets threshold.
  ##
  ## If no ancestor, including the root, reaches threshold, returns 0
  ## (unclassified). Kraken 0 runs contribute to Q but not C; A runs contribute
  ## to neither.
  validateConfidenceThreshold(threshold)

  if originalTaxId == 0'u32:
    return 0'u32

  let q = signature.queriedKmerCount()
  if q == 0:
    return 0'u32

  let counts = signature.cladeKmerCounts(taxonomy)
  var current = originalTaxId

  while true:
    let score = counts.getOrDefault(current, 0).float / q.float
    if score >= threshold:
      return current

    var parent: uint32
    if not taxonomy.parentOf(current, parent) or parent == current:
      return 0'u32

    current = parent

proc taxIdAtConfidence*(
    signature: KrakenSignature,
    taxonomy: TaxonomyIndex,
    threshold: float
  ): uint32 =
  ## Infer the starting taxid from the signature, then apply taxIdAtConfidence.
  ## Prefer the overload that accepts Kraken2's original column-3 taxid when it
  ## is available.
  signature.taxIdAtConfidence(signature.deepestTaxId(taxonomy), taxonomy, threshold)

proc sameParent*(taxonomy: TaxonomyIndex, a, b: uint32): bool =
  ## Return true when both taxids have the same direct parent.
  var parentA, parentB: uint32
  taxonomy.parentOf(a, parentA) and
    taxonomy.parentOf(b, parentB) and
    parentA == parentB

proc isAncestorOf*(taxonomy: TaxonomyIndex, ancestor, child: uint32): bool =
  ## Return true when ancestor is child or an ancestor of child.
  if ancestor == child:
    return true

  var current = child
  while current in taxonomy.parent:
    current = taxonomy.parent[current]
    if current == ancestor:
      return true

  false

proc lowestCommonAncestor*(
    taxonomy: TaxonomyIndex,
    a, b: uint32,
    lca: var uint32
  ): bool =
  ## Find the lowest common ancestor of two taxids.
  if a == b:
    lca = a
    return true

  if a notin taxonomy.depth or b notin taxonomy.depth:
    return false

  var x = a
  var y = b
  var depthX = taxonomy.depth[x]
  var depthY = taxonomy.depth[y]

  while depthX > depthY:
    if x notin taxonomy.parent:
      return false
    x = taxonomy.parent[x]
    dec depthX

  while depthY > depthX:
    if y notin taxonomy.parent:
      return false
    y = taxonomy.parent[y]
    dec depthY

  while x != y:
    if x notin taxonomy.parent or y notin taxonomy.parent:
      return false
    x = taxonomy.parent[x]
    y = taxonomy.parent[y]

  lca = x
  true

proc exactTokenSimilarity*(a, b: KmerToken): float {.inline.} =
  if a == b: 1.0 else: 0.0

proc taxonomicTokenSimilarity*(
    a, b: KmerToken,
    taxonomy: TaxonomyIndex
  ): float =
  ## Token similarity using exact matching for 0/A and tree distance for taxids.
  ##
  ## For taxids, this uses the common normalized tree score:
  ##
  ##   2 * depth(LCA(a, b)) / (depth(a) + depth(b))
  ##
  ## Exact taxid equality still returns 1 even when the taxid is absent from
  ## the taxonomy index.
  if a == b:
    return 1.0

  if a.kind != ktkTaxon or b.kind != ktkTaxon:
    return 0.0

  var lca: uint32
  if not taxonomy.lowestCommonAncestor(a.taxId, b.taxId, lca):
    return 0.0

  let depthA = taxonomy.depthOf(a.taxId)
  let depthB = taxonomy.depthOf(b.taxId)
  let depthLca = taxonomy.depthOf(lca)
  let denominator = depthA + depthB

  if depthA < 0 or depthB < 0 or depthLca < 0 or denominator <= 0:
    return 0.0

  (2.0 * depthLca.float) / denominator.float

proc positionalSimilarity*(
    left, right: KrakenSignature,
    tokenSimilarity: TokenSimilarityProc,
    missingScore = 0.0
  ): float =
  ## Position-aware signature similarity without expanding the RLE stream.
  ##
  ## Mates are compared by index. If one signature has an extra mate or longer
  ## mate, those missing positions contribute missingScore, normally 0.
  var score = 0.0
  var total = 0
  let mateN = max(left.mates.len, right.mates.len)

  for mateIdx in 0 ..< mateN:
    if mateIdx >= left.mates.len:
      let n = right.mates[mateIdx].kmerCount()
      score += missingScore * n.float
      total += n
      continue

    if mateIdx >= right.mates.len:
      let n = left.mates[mateIdx].kmerCount()
      score += missingScore * n.float
      total += n
      continue

    let a = left.mates[mateIdx]
    let b = right.mates[mateIdx]
    var i = 0
    var j = 0
    var offsetA = 0
    var offsetB = 0

    while i < a.len and j < b.len:
      let remainingA = a[i].count - offsetA
      let remainingB = b[j].count - offsetB
      let step = min(remainingA, remainingB)

      score += step.float * tokenSimilarity(a[i].token, b[j].token)
      total += step

      offsetA += step
      offsetB += step

      if offsetA == a[i].count:
        inc i
        offsetA = 0

      if offsetB == b[j].count:
        inc j
        offsetB = 0

    while i < a.len:
      let remaining = a[i].count - offsetA
      score += missingScore * remaining.float
      total += remaining
      inc i
      offsetA = 0

    while j < b.len:
      let remaining = b[j].count - offsetB
      score += missingScore * remaining.float
      total += remaining
      inc j
      offsetB = 0

  if total == 0:
    0.0
  else:
    score / total.float

proc positionalSimilarity*(
    left, right: KrakenSignature,
    missingScore = 0.0
  ): float =
  ## Exact position-aware similarity. Taxids, 0, and A are all literal tokens.
  positionalSimilarity(
    left,
    right,
    proc(a, b: KmerToken): float = exactTokenSimilarity(a, b),
    missingScore
  )

proc taxonomicPositionalSimilarity*(
    left, right: KrakenSignature,
    taxonomy: TaxonomyIndex,
    missingScore = 0.0
  ): float =
  ## Position-aware similarity where real taxids get partial credit based on
  ## taxonomy. Kraken 0 and A remain exact literal signature tokens.
  positionalSimilarity(
    left,
    right,
    proc(a, b: KmerToken): float = taxonomicTokenSimilarity(a, b, taxonomy),
    missingScore
  )

proc tokenCounts*(
    signature: KrakenSignature,
    keepMates = true
  ): Table[SignatureCountKey, int] =
  ## Count token composition. By default R1 and R2 positions are kept separate.
  result = initTable[SignatureCountKey, int]()

  for mateIdx, mate in signature.mates:
    let keyMate = if keepMates: mateIdx else: 0
    for run in mate:
      let key = SignatureCountKey(mate: keyMate, token: run.token)
      result[key] = result.getOrDefault(key, 0) + run.count

proc weightedJaccardSimilarity*(
    left, right: KrakenSignature,
    keepMates = true
  ): float =
  ## Composition similarity over token counts, including 0 and A.
  ##
  ## This ignores exact positions but can tolerate small local shifts.
  let a = left.tokenCounts(keepMates)
  let b = right.tokenCounts(keepMates)
  var keys = initHashSet[SignatureCountKey]()

  for key in a.keys:
    keys.incl(key)
  for key in b.keys:
    keys.incl(key)

  if keys.len == 0:
    return 0.0

  var intersection = 0
  var union = 0

  for key in keys:
    let countA = a.getOrDefault(key, 0)
    let countB = b.getOrDefault(key, 0)
    intersection += min(countA, countB)
    union += max(countA, countB)

  if union == 0:
    0.0
  else:
    intersection.float / union.float

proc distanceFromSimilarity*(similarity: float): float {.inline.} =
  ## Convert a [0, 1] similarity to a [0, 1] distance.
  1.0 - similarity

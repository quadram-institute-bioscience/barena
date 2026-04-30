## Parse kraken2-inspect output and expand taxon IDs to include all descendants.

import std/[tables, sets, strutils]

type TaxTree* = Table[uint32, seq[uint32]]

proc parseInspectDb*(path: string, nodeCount: var int): TaxTree =
  ## Build a parent→children map from a kraken2-inspect output file.
  ## Lines starting with '#' are skipped. Depth is inferred from leading spaces
  ## in the name field (2 spaces per level). Sets nodeCount to the total number
  ## of taxon entries parsed.
  result = initTable[uint32, seq[uint32]]()
  var stack: seq[tuple[depth: int, taxId: uint32]] = @[]
  nodeCount = 0

  for line in lines(path):
    if line.len == 0 or line[0] == '#':
      continue

    # Fields: pct | clade_reads | direct_reads | rank | taxid | name
    let parts = line.split('\t')
    if parts.len < 6:
      continue

    let taxId = parseUInt(parts[4].strip()).uint32
    inc nodeCount

    # Indentation of name field encodes depth (2 spaces per level)
    var depth = 0
    for c in parts[5]:
      if c == ' ': inc depth
      else: break
    depth = depth div 2

    # Pop stack until we find the parent scope
    while stack.len > 0 and stack[^1].depth >= depth:
      discard stack.pop()

    if stack.len > 0:
      result.mgetOrPut(stack[^1].taxId, @[]).add(taxId)

    stack.add((depth: depth, taxId: taxId))

proc parseInspectDb*(path: string): TaxTree =
  var n: int
  result = parseInspectDb(path, n)

proc getDescendants*(tree: TaxTree, roots: seq[uint32]): HashSet[uint32] =
  ## BFS over tree to collect roots and all their descendants.
  result = initHashSet[uint32]()
  var queue = roots
  var head = 0
  while head < queue.len:
    let id = queue[head]; inc head
    if result.containsOrIncl(id): continue
    if id in tree:
      for child in tree[id]:
        queue.add(child)

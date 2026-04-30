## Fast parser for kraken2 raw output lines.
## Format: C/U \t read_name \t taxid \t length \t kmer_hits
## We only ever need fields 0 (classified flag) and 2 (taxid).

proc parseKrakenLine*(line: string, taxId: var uint32): bool =
  ## Returns true if the read is classified (C). Sets taxId.
  ## Skips name field entirely; finds taxid by counting two tabs.
  if line.len < 3:
    taxId = 0
    return false

  result = line[0] == 'C'

  var tabs = 0
  var i = 1
  while i < line.len:
    if line[i] == '\t':
      inc tabs
      if tabs == 2:
        taxId = 0
        inc i
        while i < line.len and line[i] != '\t' and line[i] != '\r' and line[i] != '\n':
          taxId = taxId * 10'u32 + uint32(ord(line[i]) - ord('0'))
          inc i
        return
    inc i

  taxId = 0

proc parseKrakenReadName*(line: string): string =
  ## Extract field 1 (read name) — only called in --strict mode.
  var i = 0
  while i < line.len and line[i] != '\t': inc i  # skip field 0
  inc i                                           # skip tab
  let start = i
  while i < line.len and line[i] != '\t': inc i  # span field 1
  return line[start ..< i]

proc stripPairSuffix*(name: string): string =
  ## Remove /1 or /2 read-pair suffix for name comparison.
  if name.len > 2 and name[^2] == '/' and (name[^1] == '1' or name[^1] == '2'):
    return name[0 ..< name.len - 2]
  return name

## Output file abstraction: handles plain and gzip output transparently.
## Uses a shared write buffer (passed by caller) to avoid per-record allocation.

import std/[os, strutils]
import readfx
import gzio

type
  OutKind* = enum okNone, okPlain, okGzip

  OutFile* = object
    kind*: OutKind
    plain: File
    gz: GzWriter

proc openOut*(path: string): OutFile =
  if path.len == 0:
    return OutFile(kind: okNone)
  let dir = path.parentDir
  if dir.len > 0:
    createDir(dir)
  if path.endsWith(".gz"):
    OutFile(kind: okGzip, gz: openGzWriter(path))
  else:
    OutFile(kind: okPlain, plain: open(path, fmWrite))

proc isActive*(o: OutFile): bool =
  o.kind != okNone

proc writeFq*(o: var OutFile, r: FQRecord, buf: var string) =
  ## Write a FASTQ record using a caller-owned buffer (avoids re-allocation).
  if o.kind == okNone: return
  buf.setLen(0)
  buf.add('@')
  buf.add(r.name)
  if r.comment.len > 0:
    buf.add(' ')
    buf.add(r.comment)
  buf.add('\n')
  buf.add(r.sequence)
  buf.add("\n+\n")
  buf.add(r.quality)
  buf.add('\n')
  case o.kind
  of okNone: discard
  of okPlain: o.plain.write(buf)
  of okGzip:  o.gz.write(buf)

proc close*(o: var OutFile) =
  case o.kind
  of okNone: discard
  of okPlain:
    if o.plain != nil: o.plain.close()
  of okGzip:
    o.gz.close()

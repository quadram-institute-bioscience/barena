## Thin zlib wrapper for writing gzip output files.

when defined(windows):
  const libz = "zlib1.dll"
elif defined(macosx):
  const libz = "libz.dylib"
else:
  const libz = "libz.so.1"

type GzWriteHandle = pointer

proc gz_open(path: cstring, mode: cstring): GzWriteHandle
    {.cdecl, dynlib: libz, importc: "gzopen".}
proc gz_write(f: GzWriteHandle, buf: pointer, len: cuint): cint
    {.cdecl, dynlib: libz, importc: "gzwrite".}
proc gz_close(f: GzWriteHandle): cint
    {.cdecl, dynlib: libz, importc: "gzclose".}

type
  GzWriter* = object
    handle: GzWriteHandle
    path: string

proc openGzWriter*(path: string): GzWriter =
  result.path = path
  result.handle = gz_open(path.cstring, "wb")
  if result.handle == nil:
    raise newException(IOError, "Cannot open gzip file for writing: " & path)

proc write*(w: var GzWriter, data: string) =
  if data.len == 0: return
  let n = gz_write(w.handle, unsafeAddr data[0], cuint(data.len))
  if n != cint(data.len):
    raise newException(IOError, "gzwrite error on: " & w.path)

proc close*(w: var GzWriter) =
  if w.handle != nil:
    discard gz_close(w.handle)
    w.handle = nil

#! /usr/bin/env python

import clang.cindex
import sys
import clang
from clang.cindex import CursorKind, TokenKind, StorageClass

clang.cindex.Config.set_library_file("/usr/lib/llvm-5.0/lib/libclang.so.1")
human_readable=False

def fully_qualified(c):
    if c is None:
        return ''
    elif c.kind == CursorKind.TRANSLATION_UNIT:
        return ''
    else:
        res = fully_qualified(c.semantic_parent)
        if res != '':
            return res + '::' + c.spelling
    return c.spelling

def output(name, start, end):
    if human_readable:
        print "{} {} {}".format(name, start, end)
    else:
        print "{},{}p".format(start,end)


if __name__ == "__main__":
    idx = clang.cindex.Index.create()
    tu = idx.parse(sys.argv[1], args=['-I../../src/','-Inewfiles/', '-DWITH_ZLIB'])

    headers = set()
    for idx, c in enumerate(tu.cursor.walk_preorder()):


        if c.kind == CursorKind.CALL_EXPR and c.referenced.get_definition() is None and c.referenced.storage_class is not StorageClass.STATIC:
            #print("\n".join(dir(c)))
            #sys.exit(0)
            #if str(c.location.file).endswith("c"):
            headers.add(str(c.referenced.location.file))

    print "\n".join(map(str, headers))
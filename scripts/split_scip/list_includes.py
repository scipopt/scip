#! /usr/bin/env python

import clang.cindex
import sys
import clang
import os.path
from clang.cindex import CursorKind, TokenKind, TypeKind, StorageClass

clang.cindex.Config.set_library_file("/usr/lib/llvm-5.0/lib/libclang.so.1")
includedirs = ['../../src/','newfiles/']

def remove_prefix(p,s):
    if s.startswith(p):
        return s[len(p):]
    return s

def get_headerset(file):
    idx = clang.cindex.Index.create()
    tu = idx.parse(file, args= map(lambda x: '-I' + x, includedirs) + ['-DSCIP_WITH_ZLIB'])

    headers = set()
    for idx, c in enumerate(tu.cursor.walk_preorder()):
        if str(c.location.file) != file:
            continue

        referenced = None

        if c.kind == CursorKind.PARM_DECL:
            thetype = c.type
            while thetype == TypeKind.POINTER:
                thetype = thetype.get_pointee()
            if thetype.kind == TypeKind.TYPEDEF:
                referenced = thetype.get_declaration()
        elif c.kind == CursorKind.FUNCTION_DECL:
            thetype = c.result_type
            while thetype == TypeKind.POINTER:
                thetype = thetype.get_pointee()
            if thetype.kind == TypeKind.TYPEDEF:
                referenced = thetype.get_declaration()
        elif c.kind == CursorKind.CALL_EXPR or c.kind == CursorKind.TYPE_REF or CursorKind.MEMBER_REF_EXPR or CursorKind.DECL_REF_EXPR:
            referenced = c.referenced

        if referenced is not None:
            #print("\n".join(dir(c)))
            #sys.exit(0)
            headerfile = str(referenced.location.file)

            if not headerfile.endswith("c"):
                for p in includedirs:
                    headerfile = remove_prefix(p, headerfile)

                headers.add(headerfile)
    return headers

if __name__ == "__main__":
    inclset = get_headerset(sys.argv[1])
    includelist = []
    if sys.argv[1].endswith(".c"):
        headerfile = sys.argv[1][:-1] + 'h'
        if os.path.isfile(headerfile):
            headerinclset = get_headerset(headerfile)
        else:
            headerinclset = set()

        for p in includedirs:
            headerfile = remove_prefix(p, headerfile)

        for f in headerinclset:
            if f in inclset:
                inclset.remove(f)

        inclset.add(headerfile)
    else:
        includelist.append('#include "scip/def.h"')

    for h in inclset:
        # for system headers only add them if they are not included by scip/def.h and also handle
        # that math functions are defined in mathcalls.h on our linux systems
        if h.startswith("/usr/include/"):
            if h.endswith("mathcalls.h"):
                h = "math.h"
            else:
                h = remove_prefix("/usr/include/", h)
            if h not in ["stdio.h", "stdint.h", "math.h", "limits.h", "float.h", "assert.h"]:
                includelist.append("#include <{}>".format(remove_prefix("/usr/include/", h)))
        else:
            includelist.append('#include "{}"'.format(h))

    print "\n".join(sorted(includelist))

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

def remove_prefix(p,s):
    if s.startswith(p):
        return s[len(p):]
    return s

if __name__ == "__main__":
    idx = clang.cindex.Index.create()
    includedirs = ['../../src/','newfiles/']
    tu = idx.parse(sys.argv[1], args= map(lambda x: '-I' + x, includedirs) + ['-DWITH_ZLIB'])

    headers = set()
    for idx, c in enumerate(tu.cursor.walk_preorder()):


        if c.kind == CursorKind.CALL_EXPR and c.referenced.get_definition() is None and c.referenced.storage_class is not StorageClass.STATIC:
            #print("\n".join(dir(c)))
            #sys.exit(0)
            headerfile = str(c.referenced.location.file)

            if not headerfile.endswith("c"):
                for p in includedirs:
                    headerfile = remove_prefix(p, headerfile)

                headers.add(headerfile)

    includelist = ['#include "scip/def.h"']
    for h in headers:
        # for system headers only add them if they are not included by scip/def.h and also handle
        # that math functions are defined in mathcalls.h on our linux systems
        if h.startswith("/usr/include/"):
            if h.endswith("mathcalls.h"):
                h = "math.h"
            else:
                h = remove_prefix("/usr/include/", h)
            if h == "string.h":
                includelist.append('#include <string.h>\n#if defined(_WIN32) || defined(_WIN64)\n#else\n#include <strings.h> /*lint --e{766}*/\n#endif')
            elif h not in ["stdio.h", "stdint.h", "math.h", "limits.h", "float.h", "assert.h"]:
                includelist.append("#include <{}>".format(remove_prefix("/usr/include/", h)))
        else:
            includelist.append('#include "{}"'.format(h))

    print "\n".join(sorted(includelist))
